"""
peptoid_decompose.py

Parse a full peptoid / beta-peptoid SMILES (linear or macrocyclic) and:
  1. Detect topology (linear / cyclic)
  2. Detect backbone type per residue (alpha: NCC(=O)  /  beta: NCCC(=O))
  3. Extract each residue, preserving all substituents on N and alpha-/beta-carbons
  4. Save a 2D structure image of the whole molecule
  5. Save a 2D grid image showing every individual residue
  6. Save a 3D conformer as an SDF mol-object file

Supported residue types
-----------------------
Alpha (peptoid / N-substituted glycine and derivatives):
    [N;H0,H1]CC(=O)          <- N may be tertiary (no H) OR secondary (1 H)
                                 alpha-C may carry substituents

Beta (N-substituted / alpha- or beta-substituted beta-amino acids):
    [N;H0,H1]CCC(=O)         <- same N generalisation; beta-C, alpha-C may carry substituents

Usage
-----
Linear alpha peptoid:
    python peptoid_decompose.py \\
      -s "NCC(=O)N(Cc1ccccc1)CC(=O)N(CC(C)C)CC(=O)N(CCCO)CC(=O)O"

Macrocyclic alpha peptoid:
    python peptoid_decompose.py \\
      -s "O=C1CN(Cc2ccccc2)CC(=O)N(Cc2ccccc2)CC(=O)N(CC(C)C)CC(=O)N(CCCO)C1"

Beta-peptoid (mixed N- and alpha-substituted):
    python peptoid_decompose.py \\
      -s "NCC(C)C(=O)N(Cc1ccccc1)CCC(=O)N(C)C(CC)C(=O)O"

Custom output prefix:
    python peptoid_decompose.py -s "<SMILES>" -o my_peptoid

Outputs (all prefixed by -o / --output, default "peptoid"):
    <prefix>_structure.png   - 2D image of the full molecule
    <prefix>_residues.png    - 2D grid image of each residue
    <prefix>_3d.sdf          - 3D conformer (SDF, with hydrogens)
"""

import argparse
import copy
import re
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, RWMol
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.Draw import MolToFile, MolsToGridImage


# ---------------------------------------------------------------------------
# Backbone SMARTS patterns
# ---------------------------------------------------------------------------

# Alpha residue: N (any substitution) - CH2_or_substituted - carbonyl C
# [N;H0,H1] covers both tertiary (peptoid) and secondary (alpha-subst.) nitrogens.
# The alpha-carbon is [CX4] so it matches CH2, CHR, CR2 alike.
_ALPHA_PATTERN = Chem.MolFromSmarts("[N;H0,H1:1][CX4:2]C(=O)")

# Beta residue: N - beta-C - alpha-C - carbonyl C
# Both carbons are [CX4] to allow any substitution.
_BETA_PATTERN  = Chem.MolFromSmarts("[N;H0,H1:1][CX4:2][CX4:3]C(=O)")

# Map pattern -> label
_PATTERNS = [
    (_BETA_PATTERN,  "beta"),   # beta checked first (longer chain, avoids subset match)
    (_ALPHA_PATTERN, "alpha"),
]


# ---------------------------------------------------------------------------
# Backbone collection
# ---------------------------------------------------------------------------

def _collect_all_matches(mol):
    """
    Return a list of dicts, one per unique backbone-N, containing:
        n_idx, backbone_type, atom_indices (all backbone atoms for this match)
    Beta matches take priority over alpha so overlapping patterns are not
    double-counted.
    """
    claimed_n = {}   # n_idx -> match dict

    for pattern, btype in _PATTERNS:
        for m in mol.GetSubstructMatches(pattern):
            n_idx = m[0]
            if n_idx in claimed_n:
                continue   # already claimed by longer (beta) pattern
            claimed_n[n_idx] = {
                "n_idx":          n_idx,
                "backbone_type":  btype,
                "atom_indices":   set(m),   # N, C(s), carbonyl-C implicitly in m
            }

    return list(claimed_n.values())


def _build_backbone_set(mol, match_list):
    """
    Expand atom_indices to include the carbonyl oxygen(s) for each match, then
    absorb any cap-group carbonyl attached directly to a residue N.

    Returns a frozenset of all backbone atom indices.
    """
    all_backbone = set()
    n_set = {e["n_idx"] for e in match_list}

    # Seed with all pattern atoms + their carbonyl oxygens
    for entry in match_list:
        all_backbone.update(entry["atom_indices"])
        # Find carbonyl C: last atom in atom_indices that has a =O neighbour
        for idx in entry["atom_indices"]:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(idx, nb.GetIdx())
                if (nb.GetAtomicNum() == 8
                        and bond.GetBondType() == Chem.BondType.DOUBLE):
                    all_backbone.add(nb.GetIdx())

    # Absorb cap-group carbonyls bonded directly to residue-N atoms
    for n_idx in n_set:
        for nb in mol.GetAtomWithIdx(n_idx).GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in all_backbone:
                continue
            if mol.GetAtomWithIdx(nb_idx).GetAtomicNum() != 6:
                continue
            # Is this C a carbonyl?
            for nb2 in mol.GetAtomWithIdx(nb_idx).GetNeighbors():
                bond = mol.GetBondBetweenAtoms(nb_idx, nb2.GetIdx())
                if (nb2.GetAtomicNum() == 8
                        and bond.GetBondType() == Chem.BondType.DOUBLE):
                    # BFS the full cap group (don't cross residue-N boundaries)
                    cap_queue = [nb_idx]
                    while cap_queue:
                        cur = cap_queue.pop(0)
                        if cur in all_backbone:
                            continue
                        all_backbone.add(cur)
                        for nb3 in mol.GetAtomWithIdx(cur).GetNeighbors():
                            ni3 = nb3.GetIdx()
                            if ni3 not in n_set and ni3 not in all_backbone:
                                cap_queue.append(ni3)
                    break

    return frozenset(all_backbone)


# ---------------------------------------------------------------------------
# Topology detection
# ---------------------------------------------------------------------------

def detect_topology(mol, match_list):
    """Return 'cyclic' if any backbone atom sits in a ring, else 'linear'."""
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    for entry in match_list:
        if entry["atom_indices"] & ring_atoms:
            return "cyclic"
    return "linear"


# ---------------------------------------------------------------------------
# Sidechain / substituent extraction
# ---------------------------------------------------------------------------

def _bfs_from(mol, start_idx, forbidden):
    """Return all atom indices reachable from start_idx without crossing forbidden."""
    visited = set()
    queue = [start_idx]
    while queue:
        cur = queue.pop(0)
        if cur in visited or cur in forbidden:
            continue
        visited.add(cur)
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in visited and ni not in forbidden:
                queue.append(ni)
    return visited


def _substituents_on(mol, atom_idx, all_backbone):
    """
    Return a list of substituent atom-index sets hanging off atom_idx into
    non-backbone territory.  Each set is one substituent branch.
    """
    subs = []
    for nb in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        ni = nb.GetIdx()
        if ni not in all_backbone:
            branch = _bfs_from(mol, ni, all_backbone)
            subs.append(branch)
    return subs


def _smiles_for_atoms(mol, keep_atoms, anchor_idx):
    """
    Build a SMILES for a fragment defined by keep_atoms, replacing anchor_idx
    with a [1*] dummy so that stereo context is preserved, then strip the dummy.
    """
    emol = RWMol(copy.deepcopy(mol))

    anchor_atom = emol.GetAtomWithIdx(anchor_idx)
    anchor_atom.SetAtomicNum(0)
    anchor_atom.SetIsotope(1)
    anchor_atom.SetNoImplicit(True)

    remove = sorted(set(range(mol.GetNumAtoms())) - keep_atoms, reverse=True)
    for idx in remove:
        emol.RemoveAtom(idx)

    try:
        Chem.SanitizeMol(emol)
    except Exception:
        return None

    smi = Chem.MolToSmiles(emol.GetMol(), isomericSmiles=True)
    smi = smi.replace("[1*]", "")
    smi = re.sub(r"^\(*-?", "", smi)
    smi = re.sub(r"\)*$", "", smi)
    smi = smi.strip("()")
    smi = re.sub(r"^[-=#]", "", smi)
    return smi if smi else None


def _get_substituent_smiles(mol, atom_idx, all_backbone):
    """
    Return a list of SMILES strings for all substituents on atom_idx that
    extend outside the backbone.
    """
    results = []
    for nb in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        ni = nb.GetIdx()
        if ni in all_backbone:
            continue
        branch_atoms = _bfs_from(mol, ni, all_backbone)
        keep = branch_atoms | {atom_idx}
        smi = _smiles_for_atoms(mol, keep, atom_idx)
        if smi:
            results.append(smi)
    return results


# ---------------------------------------------------------------------------
# Residue SMILES builder
# ---------------------------------------------------------------------------

def _build_residue_smiles(entry, n_sub, alpha_sub, beta_sub=None):
    """
    Construct a human-readable residue SMILES in the form:
      Alpha: *N(<n_sub>)C(<alpha_sub>...)C(=O)*
      Beta:  *N(<n_sub>)C(<beta_sub>...)C(<alpha_sub>...)C(=O)*

    Falls back gracefully when substitution lists are empty.
    """
    def _fmt_subs(subs):
        return "".join(f"({s})" for s in subs)

    n_part    = f"N{_fmt_subs(n_sub)}"
    al_part   = f"C{_fmt_subs(alpha_sub)}"
    tail      = "C(=O)*"

    if entry["backbone_type"] == "alpha":
        return f"*{n_part}{al_part}{tail}"
    else:
        be_part = f"C{_fmt_subs(beta_sub or [])}"
        return f"*{n_part}{be_part}{al_part}{tail}"


# ---------------------------------------------------------------------------
# Main extraction (importable)
# ---------------------------------------------------------------------------

def extract_residues(mol):
    """
    Decompose a peptoid / beta-peptoid molecule into its residues.

    Parameters
    ----------
    mol : rdkit.Chem.Mol

    Returns
    -------
    residues : list of dict
        Each dict contains:
          'index'           - 1-based residue number (int)
          'backbone_type'   - 'alpha' or 'beta' (str)
          'n_atom_idx'      - backbone N atom index in mol (int)
          'n_substituents'  - list of sidechain SMILES on N (list[str])
          'alpha_substituents' - list of substituent SMILES on alpha-C (list[str])
          'beta_substituents'  - list of substituent SMILES on beta-C (list[str])
                                 (empty list for alpha residues)
          'residue_smiles'  - full residue notation with * attachment points (str)
    topology : str
        'linear' or 'cyclic'

    Raises
    ------
    ValueError
        If no recognised residues are found.
    """
    match_list = _collect_all_matches(mol)
    if not match_list:
        raise ValueError(
            "No peptoid / beta-peptoid residues detected.\n"
            "Expected at least one [N;H0,H1]CC(=O)  (alpha) or\n"
            "                       [N;H0,H1]CCC(=O) (beta) unit."
        )

    topology      = detect_topology(mol, match_list)
    all_backbone  = _build_backbone_set(mol, match_list)

    residues = []
    for i, entry in enumerate(match_list):
        m_atoms = sorted(entry["atom_indices"])
        n_idx   = entry["n_idx"]

        # Identify alpha-C and (for beta) beta-C from the SMARTS match ordering.
        # atom_indices is a set; we reconstruct order by walking the backbone chain.
        # N -> next C -> (next C for beta) -> carbonyl C
        backbone_chain = _walk_backbone_chain(mol, entry, all_backbone)

        if entry["backbone_type"] == "alpha":
            # backbone_chain: [N, alpha_C, carbonyl_C]
            alpha_c_idx = backbone_chain[1] if len(backbone_chain) >= 3 else None
            beta_c_idx  = None
        else:
            # backbone_chain: [N, beta_C, alpha_C, carbonyl_C]
            beta_c_idx  = backbone_chain[1] if len(backbone_chain) >= 4 else None
            alpha_c_idx = backbone_chain[2] if len(backbone_chain) >= 4 else None

        n_sub    = _get_substituent_smiles(mol, n_idx,      all_backbone)
        al_sub   = _get_substituent_smiles(mol, alpha_c_idx, all_backbone) if alpha_c_idx is not None else []
        be_sub   = _get_substituent_smiles(mol, beta_c_idx,  all_backbone) if beta_c_idx  is not None else []

        res_smi  = _build_residue_smiles(entry, n_sub, al_sub, be_sub)

        residues.append(
            {
                "index":               i + 1,
                "backbone_type":       entry["backbone_type"],
                "n_atom_idx":          n_idx,
                "n_substituents":      n_sub,
                "alpha_substituents":  al_sub,
                "beta_substituents":   be_sub,
                "residue_smiles":      res_smi,
            }
        )

    return residues, topology


def _walk_backbone_chain(mol, entry, all_backbone):
    """
    Walk from the residue N through the backbone to return an ordered list:
      alpha: [N_idx, alpha_C_idx, carbonyl_C_idx]
      beta:  [N_idx, beta_C_idx,  alpha_C_idx,  carbonyl_C_idx]

    The carbonyl C is identified as the C in backbone with a =O neighbour.
    """
    n_idx = entry["n_idx"]
    local = entry["atom_indices"]

    # Collect only atoms in this residue's match (not the full backbone set)
    # and walk starting from N.
    chain = [n_idx]
    visited = {n_idx}
    current = n_idx

    for _ in range(3):   # at most 3 steps (alpha: 2, beta: 3 after N)
        found_next = False
        for nb in mol.GetAtomWithIdx(current).GetNeighbors():
            ni = nb.GetIdx()
            if ni in local and ni not in visited:
                chain.append(ni)
                visited.add(ni)
                current = ni
                found_next = True
                break
        if not found_next:
            break

    return chain


# ---------------------------------------------------------------------------
# Output helpers  (unchanged API, updated for new residue dict shape)
# ---------------------------------------------------------------------------

def save_structure_image(mol, path, size=(1200, 700)):
    """Save a 2D depiction of the complete molecule."""
    mol2d = copy.deepcopy(mol)
    Compute2DCoords(mol2d)
    MolToFile(mol2d, path, size=size)


def save_residue_grid(residues, path, cols=4, sub_size=(420, 320)):
    """
    Draw each residue in a labelled grid image.
    Uses residue_smiles with * -> [*] for RDKit parsing.
    """
    mols, legends = [], []
    for r in residues:
        rs = r["residue_smiles"].replace("*", "[*]")
        m = Chem.MolFromSmiles(rs)
        if m is None:
            # Fallback: try just the N-substituents joined
            fallback = ".".join(r["n_substituents"]) if r["n_substituents"] else None
            m = Chem.MolFromSmiles(fallback) if fallback else None
        if m:
            Compute2DCoords(m)
            mols.append(m)
            btype = r["backbone_type"]
            legends.append(
                f"Residue {r['index']} ({btype})\n{r['residue_smiles']}"
            )

    if not mols:
        print("  WARNING: no residue structures could be rendered for grid image.")
        return

    img = MolsToGridImage(
        mols,
        molsPerRow=min(cols, len(mols)),
        subImgSize=sub_size,
        legends=legends,
    )
    img.save(path)


def save_3d_sdf(mol, path):
    """
    Embed a 3D conformer (ETKDGv3 + MMFF minimisation) and write to SDF.
    Returns True on success, False if embedding fails.
    """
    mol3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol3d, params)

    if result == -1:
        result = AllChem.EmbedMolecule(mol3d, AllChem.EmbedParameters())

    if result == -1:
        print("  WARNING: 3D embedding failed - SDF not written.")
        return False

    ff_result = AllChem.MMFFOptimizeMolecule(mol3d, maxIters=2000)
    if ff_result == 1:
        print("  NOTE: MMFF did not fully converge (geometry may still be usable).")

    writer = Chem.SDWriter(path)
    writer.write(mol3d)
    writer.close()
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Decompose a peptoid / beta-peptoid SMILES (linear or macrocyclic)\n"
            "into its residues, generate 2D images, and a 3D SDF mol object.\n\n"
            "Supports:\n"
            "  - Alpha backbones  [N;H0,H1]CC(=O)  with any N- or alpha-C-substitution\n"
            "  - Beta  backbones  [N;H0,H1]CCC(=O) with any N-, beta-C-, or alpha-C-substitution\n"
            "  - Mixed alpha/beta sequences\n"
            "  - Linear and macrocyclic topologies"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-s", "--smiles",
        required=True,
        type=str,
        metavar="SMILES",
        help="Full-structure SMILES of the peptoid / beta-peptoid.",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="peptoid",
        metavar="PREFIX",
        help=(
            "Output filename prefix (default: 'peptoid').\n"
            "Three files are written:\n"
            "  <prefix>_structure.png  - 2D whole-molecule image\n"
            "  <prefix>_residues.png   - 2D residue grid image\n"
            "  <prefix>_3d.sdf         - 3D conformer (SDF with H atoms)"
        ),
    )
    args = parser.parse_args()

    # Parse molecule
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        print(f"ERROR: could not parse SMILES: {args.smiles}", file=sys.stderr)
        sys.exit(1)

    print(f"\nSMILES     : {args.smiles}")
    print(f"Atom count : {mol.GetNumAtoms()} heavy atoms")

    # Extract residues
    try:
        residues, topology = extract_residues(mol)
    except ValueError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Topology   : {topology}")
    print(f"\nResidues ({len(residues)} found):")
    for r in residues:
        btype = r["backbone_type"]
        print(f"  [{r['index']:2d}] ({btype:5s})  {r['residue_smiles']}")
        if r["n_substituents"]:
            print(f"        N-substituents    : {', '.join(r['n_substituents'])}")
        if r["alpha_substituents"]:
            print(f"        alpha-C-substituents: {', '.join(r['alpha_substituents'])}")
        if r["beta_substituents"]:
            print(f"        beta-C-substituents : {', '.join(r['beta_substituents'])}")

    # 2D whole-molecule image
    struct_path = f"{args.output}_structure.png"
    save_structure_image(mol, struct_path)
    print(f"\n2D structure image  ->  {struct_path}")

    # 2D residue grid
    grid_path = f"{args.output}_residues.png"
    save_residue_grid(residues, grid_path)
    print(f"Residue grid image  ->  {grid_path}")

    # 3D SDF
    sdf_path = f"{args.output}_3d.sdf"
    ok = save_3d_sdf(mol, sdf_path)
    if ok:
        print(f"3D mol object (SDF) ->  {sdf_path}")

    print()


if __name__ == "__main__":
    main()
