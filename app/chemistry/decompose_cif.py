"""
peptoid_decompose.py

Parse a full peptoid SMILES or a small-molecule CIF file (linear or macrocyclic)
and:
  1. Detect topology (linear / cyclic)
  2. Extract each residue as *N(sidechain)CC(=O)*  (stereochemistry preserved)
  3. Save a 2D structure image of the whole molecule
  4. Save a 2D grid image showing every individual residue
  5. Save a 3D conformer as an SDF mol-object file

CIF ingestion pipeline
----------------------
When a CIF file is supplied (-c / --cif), the following steps are applied
before decomposition:

  1. Parse the CIF with gemmi (small-molecule crystallographic reader).
  2. Remove alternate conformers / disordered rotamers: only the primary
     conformer (disorder group 1 / conformer A) is kept.  This prevents
     phantom bonds between duplicate atom positions from causing RDKit
     sanitization failures.
  3. Attempt to read explicit bond connectivity from the _geom_bond loop
     (present in CCDC-deposited structures).  Bond orders are inferred from
     interatomic distances when not stored in the CIF.
  4. Fall back to coordinate-only bond determination (RDKit
     rdDetermineBonds) when no _geom_bond data is found.
  5. Remove all counterions and solvent molecules by keeping only the
     largest connected fragment (by heavy-atom count).
  6. Strip explicit hydrogen atoms (RDKit re-adds them implicitly).
  7. Sanitize the resulting molecule and convert to SMILES, which is then
     fed into the standard decomposition pipeline.

Usage
-----
From SMILES (original behaviour):
    python peptoid_decompose.py \
      -s "NCC(=O)N(Cc1ccccc1)CC(=O)N(CC(C)C)CC(=O)N(CCCO)CC(=O)O"

From a CIF file:
    python peptoid_decompose.py -c my_peptoid.cif

Macrocyclic peptoid from CIF:
    python peptoid_decompose.py -c cyclic_peptoid.cif

Custom output prefix:
    python peptoid_decompose.py -c my_peptoid.cif -o my_output

Outputs (all prefixed by -o / --output, default "peptoid"):
    <prefix>_structure.png   - 2D image of the full molecule
    <prefix>_residues.png    - 2D grid image of each residue
    <prefix>_3d.sdf          - 3D conformer (SDF, with hydrogens)

Dependencies
------------
  pip install rdkit gemmi
"""

import argparse
import copy
import math
import re
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, RWMol
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.Draw import MolToFile, MolsToGridImage

try:
    import gemmi
    _GEMMI_AVAILABLE = True
except ImportError:
    _GEMMI_AVAILABLE = False


# ---------------------------------------------------------------------------
# Backbone identification
# ---------------------------------------------------------------------------

# Peptoid residue SMARTS: tertiary N (no H, 3 bonds) - alpha-C - carbonyl-C - O
_RESIDUE_PATTERN = Chem.MolFromSmarts("[NX3;H0]CC(=O)")


def _collect_backbone(mol, all_matches):
    """
    Return (all_backbone_atom_indices, n_backbone_dict).

    all_backbone includes:
      - Every residue N, alpha-C, carbonyl-C, and carbonyl-O
      - Any cap groups (C=O directly bonded to a residue N but not a residue
        carbonyl) - e.g. N-terminal free-amine predecessor NCC(=O)-, acetyl, etc.
    """
    n_backbone = {}   # n_idx -> {alpha_C, carbonyl_C, O, ...}
    for m in all_matches:
        n_idx = m[0]
        if n_idx not in n_backbone:
            n_backbone[n_idx] = set()
        n_backbone[n_idx].update(m[1:])

    all_backbone = set(n_backbone.keys())
    for conns in n_backbone.values():
        all_backbone.update(conns)

    # Absorb cap groups attached to residue N atoms
    for n_idx in list(n_backbone.keys()):
        for nb in mol.GetAtomWithIdx(n_idx).GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in all_backbone:
                continue
            nb_atom = mol.GetAtomWithIdx(nb_idx)
            if nb_atom.GetAtomicNum() != 6:
                continue
            # Check if this C is a carbonyl C (has a =O neighbor)
            for nb2 in nb_atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(nb_idx, nb2.GetIdx())
                if nb2.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    # BFS the entire cap group (do not cross residue-N boundary)
                    cap_queue = [nb_idx]
                    while cap_queue:
                        cur = cap_queue.pop(0)
                        if cur in all_backbone:
                            continue
                        all_backbone.add(cur)
                        for nb3 in mol.GetAtomWithIdx(cur).GetNeighbors():
                            ni = nb3.GetIdx()
                            if ni not in n_backbone and ni not in all_backbone:
                                cap_queue.append(ni)
                    break

    return all_backbone, n_backbone


def _get_sidechain_atoms(mol, n_idx, all_backbone):
    """BFS from N into non-backbone territory to collect sidechain atom indices."""
    sc = set()
    for nb in mol.GetAtomWithIdx(n_idx).GetNeighbors():
        ni = nb.GetIdx()
        if ni in all_backbone:
            continue
        queue = [ni]
        while queue:
            cur = queue.pop(0)
            if cur in sc:
                continue
            sc.add(cur)
            for nb2 in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni2 = nb2.GetIdx()
                if ni2 not in all_backbone and ni2 not in sc:
                    queue.append(ni2)
    return sc


def _sidechain_smiles(mol, n_idx, sc_atoms):
    """
    Build a SMILES string for the sidechain, preserving stereochemistry.

    Strategy: temporarily replace the residue N with an isotope-labelled dummy
    atom [1*] so that chiral atoms bonded to N retain their stereo context
    during SMILES generation, then strip the dummy token from the output string.
    """
    keep = sc_atoms | {n_idx}
    emol = RWMol(copy.deepcopy(mol))

    # Mark N as [1*] (dummy, isotope 1)
    n_edit = emol.GetAtomWithIdx(n_idx)
    n_edit.SetAtomicNum(0)
    n_edit.SetIsotope(1)
    n_edit.SetNoImplicit(True)

    # Remove all atoms not in (sidechain union {N})
    remove = sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True)
    for idx in remove:
        emol.RemoveAtom(idx)

    Chem.SanitizeMol(emol)
    smi = Chem.MolToSmiles(emol.GetMol(), isomericSmiles=True)

    # Strip the [1*] token and any orphaned punctuation left around it
    smi = smi.replace("[1*]", "")
    smi = re.sub(r"^\(*-?", "", smi)   # strip leading ( or (-
    smi = re.sub(r"\)*$", "", smi)     # strip trailing )
    smi = smi.strip("()")
    smi = re.sub(r"^[-=#]", "", smi)   # strip orphan bond char at start

    return smi if smi else "[H]"


# ---------------------------------------------------------------------------
# Topology detection
# ---------------------------------------------------------------------------

def detect_topology(mol, all_matches):
    """Return 'cyclic' if any backbone atom sits in a ring, else 'linear'."""
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    for m in all_matches:
        for idx in m:
            if idx in ring_atoms:
                return "cyclic"
    return "linear"


# ---------------------------------------------------------------------------
# Main extraction function (importable)
# ---------------------------------------------------------------------------

def extract_residues(mol):
    """
    Decompose a peptoid molecule into its residues.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Parsed molecule object (from Chem.MolFromSmiles or similar).

    Returns
    -------
    residues : list of dict
        Each dict contains:
          'index'            - 1-based residue number (int)
          'n_atom_idx'       - backbone N atom index in mol (int)
          'sidechain_smiles' - SMILES of sidechain only, stereo preserved (str)
          'residue_smiles'   - full *N(sidechain)CC(=O)* notation (str)
    topology : str
        'linear' or 'cyclic'

    Raises
    ------
    ValueError
        If no peptoid residues are found in the molecule.
    """
    all_matches = mol.GetSubstructMatches(_RESIDUE_PATTERN)
    if not all_matches:
        raise ValueError(
            "No peptoid residues detected.\n"
            "Expected at least one [NX3;H0]CC(=O) unit "
            "(tertiary N-substituted glycine backbone)."
        )

    topology = detect_topology(mol, all_matches)
    all_backbone, _ = _collect_backbone(mol, all_matches)

    # One entry per unique backbone N atom
    seen_n = {}
    for m in all_matches:
        if m[0] not in seen_n:
            seen_n[m[0]] = m
    unique_matches = list(seen_n.values())

    residues = []
    for i, m in enumerate(unique_matches):
        n_idx = m[0]
        sc_atoms = _get_sidechain_atoms(mol, n_idx, all_backbone)
        sc_smi = _sidechain_smiles(mol, n_idx, sc_atoms) if sc_atoms else "[H]"
        residues.append(
            {
                "index": i + 1,
                "n_atom_idx": n_idx,
                "sidechain_smiles": sc_smi,
                "residue_smiles": f"*N({sc_smi})CC(=O)*",
            }
        )

    return residues, topology


# ---------------------------------------------------------------------------
# CIF ingestion pipeline
# ---------------------------------------------------------------------------

# ---- Bond-order thresholds (Angstrom) ------------------------------------
# Distances below the cutoff are upgraded from SINGLE to DOUBLE (or TRIPLE).
_DOUBLE_BOND_MAX_DIST = {
    frozenset(["C", "O"]): 1.30,
    frozenset(["C", "N"]): 1.35,
    frozenset(["C", "C"]): 1.40,
    frozenset(["C", "S"]): 1.65,
    frozenset(["N", "O"]): 1.25,
    frozenset(["N", "N"]): 1.30,
    frozenset(["P", "O"]): 1.55,
}
_TRIPLE_BOND_MAX_DIST = {
    frozenset(["C", "C"]): 1.25,
    frozenset(["C", "N"]): 1.20,
}

# Common solvent / counterion formulae (by canonical SMILES) to strip if
# the user requests aggressive cleaning beyond largest-fragment selection.
_SOLVENT_SMILES = {
    "O",          # water
    "[Na+]", "[K+]", "[Li+]", "[Ca+2]", "[Mg+2]", "[Zn+2]",
    "[Cl-]", "[Br-]", "[F-]", "[I-]",
    "[NH4+]",
    "CC#N",       # acetonitrile
    "CO",         # methanol
    "CCO",        # ethanol
    "CCCO",       # 1-propanol
    "CC(C)=O",    # acetone
    "O=CO",       # formic acid
    "CC(O)=O",    # acetic acid
    "CS(C)=O",    # DMSO
    "ClCCl",      # DCM
    "ClC(Cl)Cl",  # chloroform
}


def _atom_dist(conf, i, j):
    """Euclidean distance between atom i and j in an RDKit conformer."""
    p1 = conf.GetAtomPosition(i)
    p2 = conf.GetAtomPosition(j)
    return math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)


def _infer_bond_orders_from_coords(mol):
    """
    Given a molecule whose bonds are all SINGLE (from DetermineConnectivity),
    upgrade bonds to DOUBLE or TRIPLE using interatomic distances.

    Returns a new Mol (the input is not modified).
    """
    rw = RWMol(mol)
    conf = mol.GetConformer()

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ei = mol.GetAtomWithIdx(i).GetSymbol()
        ej = mol.GetAtomWithIdx(j).GetSymbol()
        if ei == "H" or ej == "H":
            continue
        dist = _atom_dist(conf, i, j)
        pair = frozenset([ei, ej])
        if pair in _TRIPLE_BOND_MAX_DIST and dist < _TRIPLE_BOND_MAX_DIST[pair]:
            rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.TRIPLE)
        elif pair in _DOUBLE_BOND_MAX_DIST and dist < _DOUBLE_BOND_MAX_DIST[pair]:
            rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE)

    return rw.GetMol()


def _build_mol_from_geom_bonds(block, ss, sites=None):
    """
    Build an RDKit molecule using explicit bond data from the _geom_bond loop.

    Bond orders are inferred from the stored distances.  Only intra-asymmetric-
    unit bonds (symmetry code '.', '1_555', or blank) are included.

    Parameters
    ----------
    block : gemmi CIF block
    ss    : gemmi SmallStructure (for cell parameters)
    sites : list of gemmi SmallStructure.Site, optional
        If provided, use this site list instead of ss.sites (e.g. pre-filtered
        to remove alternate conformers).

    Returns None if the required CIF tags are absent.
    """
    # Check that the required tags exist
    required = [
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_distance",
    ]
    for tag in required:
        if block.find_value(tag) is None and block.find_loop(tag) is None:
            return None

    if sites is None:
        sites = ss.sites

    cell = ss.cell
    label_to_idx = {}
    rw = RWMol()
    conf_positions = []

    for i, site in enumerate(sites):
        elem = site.element.name
        if not elem or elem == "X":          # unknown element — skip
            continue
        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(elem)
        if atomic_num == 0:
            continue
        label_to_idx[site.label] = len(conf_positions)
        rw.AddAtom(Chem.Atom(atomic_num))
        pos = cell.orthogonalize(site.fract)
        conf_positions.append((pos.x, pos.y, pos.z))

    if not conf_positions:
        return None

    # Add 3D conformer so bond-order inference can use distances
    conf = Chem.Conformer(rw.GetNumAtoms())
    for i, (x, y, z) in enumerate(conf_positions):
        conf.SetAtomPosition(i, (x, y, z))
    rw.AddConformer(conf, assignId=True)

    # Read bond table — try with and without the symmetry column
    try:
        bond_table = block.find(
            ["_geom_bond_atom_site_label_1",
             "_geom_bond_atom_site_label_2",
             "_geom_bond_distance",
             "_geom_bond_site_symmetry_2"]
        )
        has_sym_col = True
    except Exception:
        bond_table = block.find(
            ["_geom_bond_atom_site_label_1",
             "_geom_bond_atom_site_label_2",
             "_geom_bond_distance"]
        )
        has_sym_col = False

    def _is_intra_sym(code):
        """Return True if a symmetry code represents the identity operation
        (i.e. the bond is within the asymmetric unit, not a crystal contact).
        Identity codes: '.', blank, '1_555', or any code whose numeric part
        is exactly 555 (P1 equivalent) with symmetry operator index 1."""
        code = code.strip()
        if code in (".", "", "1_555"):
            return True
        # Reject anything that looks like a symmetry operation, e.g. '2_655'
        # These are always non-identity (different operator or translation).
        import re as _re
        if _re.match(r"^\d+_\d+$", code):
            parts = code.split("_")
            return parts[0] == "1" and parts[1] == "555"
        return False

    for row in bond_table:
        a1, a2 = row[0].strip(), row[1].strip()
        raw_dist = row[2].strip().split("(")[0]   # strip esd, e.g. "1.234(5)"
        try:
            dist = float(raw_dist)
        except ValueError:
            continue
        if has_sym_col:
            sym = row[3].strip() if len(row) > 3 else "."
            if not _is_intra_sym(sym):
                continue                          # skip symmetry-generated contacts

        if a1 not in label_to_idx or a2 not in label_to_idx:
            continue
        i1, i2 = label_to_idx[a1], label_to_idx[a2]
        if rw.GetBondBetweenAtoms(i1, i2) is not None:
            continue                              # bond already added

        # Add bond as SINGLE for now; orders are assigned below.
        rw.AddBond(i1, i2, Chem.BondType.SINGLE)

    mol_connectivity = rw.GetMol()

    # Run ring perception so we can identify ring bonds.  Ring bonds are left
    # as SINGLE here — RDKit sanitization (called later) will perceive
    # aromaticity correctly from connectivity + atom types.  Naïvely applying
    # the distance heuristic to ring bonds tags all aromatic C–C (~1.39 Å) as
    # DOUBLE, which produces pentavalent ring carbons and sanitisation failures.
    try:
        Chem.SanitizeMol(
            mol_connectivity,
            Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
        )
    except Exception:
        pass

    ring_bond_ids = {
        bond.GetIdx()
        for bond in mol_connectivity.GetBonds()
        if bond.IsInRing()
    }

    # Upgrade only non-ring bonds to double/triple using interatomic distances.
    rw2 = RWMol(mol_connectivity)
    for bond in mol_connectivity.GetBonds():
        if bond.GetIdx() in ring_bond_ids:
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ei = mol_connectivity.GetAtomWithIdx(i).GetSymbol()
        ej = mol_connectivity.GetAtomWithIdx(j).GetSymbol()
        if ei == "H" or ej == "H":
            continue
        dist = _atom_dist(conf, i, j)
        pair = frozenset([ei, ej])
        if pair in _TRIPLE_BOND_MAX_DIST and dist < _TRIPLE_BOND_MAX_DIST[pair]:
            rw2.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.TRIPLE)
        elif pair in _DOUBLE_BOND_MAX_DIST and dist < _DOUBLE_BOND_MAX_DIST[pair]:
            rw2.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE)

    return rw2.GetMol()


def _build_mol_from_coords(ss, sites=None):
    """
    Fallback: build an RDKit molecule from 3-D coordinates only.

    Uses RDKit's rdDetermineBonds for connectivity, then upgrades bonds to
    double/triple via distance heuristics.  Hydrogen atoms in the CIF are
    included to improve bond-order determination; they are stripped later.

    Parameters
    ----------
    ss    : gemmi SmallStructure (for cell parameters)
    sites : list of gemmi SmallStructure.Site, optional
        If provided, use this site list instead of ss.sites.

    Returns None on complete failure.
    """
    if sites is None:
        sites = ss.sites

    cell = ss.cell
    atoms = []
    for site in sites:
        elem = site.element.name
        if not elem or elem == "X":
            continue
        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(elem)
        if atomic_num == 0:
            continue
        pos = cell.orthogonalize(site.fract)
        atoms.append((elem, pos.x, pos.y, pos.z))

    if not atoms:
        return None

    xyz_lines = [str(len(atoms)), "mol"]
    for elem, x, y, z in atoms:
        xyz_lines.append(f"{elem}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_block = "\n".join(xyz_lines)

    mol = Chem.MolFromXYZBlock(xyz_block)
    if mol is None:
        return None

    rdDetermineBonds.DetermineConnectivity(mol)

    # Try DetermineBondOrders first (works reliably when H atoms are present)
    has_h = any(a[0] == "H" for a in atoms)
    if has_h:
        for charge in [0, -1, 1, -2, 2]:
            try:
                mol_try = copy.deepcopy(mol)
                rdDetermineBonds.DetermineBondOrders(
                    mol_try, charge=charge, allowChargedFragments=True
                )
                Chem.SanitizeMol(mol_try)
                return mol_try
            except Exception:
                continue

    # Distance-based fallback
    mol = _infer_bond_orders_from_coords(mol)
    return mol


def _neutralize_overvalent_atoms(mol, verbose=True):
    """
    Fix atoms whose explicit valence exceeds the RDKit maximum by removing
    excess H bonds or, for charged atoms (e.g. protonated N), setting the
    appropriate formal charge.

    Common cases from crystallographic CIFs:
      - N with 4 bonds and no formal charge  → set charge +1  (ammonium)
      - O with 3 bonds and no formal charge  → set charge +1  (oxonium)
      - N-H bonds that were added as real bonds → remove the N–H bond

    Returns a new RWMol (or the original if no changes needed).
    """
    # Standard maximum valences for the elements we care about
    _MAX_VALENCE = {"N": 3, "O": 2, "S": 2, "C": 4, "P": 5}

    rw = RWMol(copy.deepcopy(mol))
    changed = False

    for atom in rw.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in _MAX_VALENCE:
            continue
        max_v = _MAX_VALENCE[sym]

        # Count explicit bonds (sum of bond orders, ignoring H neighbours
        # that are going to be removed anyway)
        explicit_valence = sum(
            int(b.GetBondTypeAsDouble())
            for b in atom.GetBonds()
        )
        if explicit_valence <= max_v:
            continue

        excess = explicit_valence - max_v

        # Strategy 1: if there are N–H or O–H bonds, remove them first
        # (orphaned H atoms from CIF that got bonded)
        h_bonds = [
            b for b in atom.GetBonds()
            if rw.GetAtomWithIdx(
                b.GetOtherAtomIdx(atom.GetIdx())
            ).GetAtomicNum() == 1
        ]
        removed = 0
        for bond in h_bonds:
            if removed >= excess:
                break
            h_idx = bond.GetOtherAtomIdx(atom.GetIdx())
            rw.RemoveBond(atom.GetIdx(), h_idx)
            removed += 1
            changed = True
        excess -= removed

        if excess <= 0:
            continue

        # Strategy 2: assign formal charge to account for remaining excess
        # N with 4 heavy bonds → quaternary ammonium [N+]
        # O with 3 bonds       → oxonium [O+]
        if sym in ("N", "O") and atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(excess)
            if verbose:
                print(
                    f"  NOTE: atom {atom.GetIdx()} ({sym}) assigned "
                    f"formal charge +{excess} to satisfy valence."
                )
            changed = True
            excess = 0

        if excess <= 0:
            continue

        # Strategy 3: for atoms where formal charge cannot fix excess valence
        # (e.g. C with 5 bonds), first try downgrading double/triple bonds to
        # single (preserves connectivity, corrects over-estimated bond orders),
        # then as a last resort remove the longest heavy-atom bond.
        heavy_bonds = [
            b for b in atom.GetBonds()
            if rw.GetAtomWithIdx(
                b.GetOtherAtomIdx(atom.GetIdx())
            ).GetAtomicNum() > 1
        ]

        # 3a: downgrade multi-order bonds (longest/highest-order first)
        try:
            conf = rw.GetConformer()
            heavy_bonds.sort(
                key=lambda b: _atom_dist(
                    conf, atom.GetIdx(), b.GetOtherAtomIdx(atom.GetIdx())
                ),
                reverse=True,
            )
        except Exception:
            pass  # no conformer — arbitrary order

        for bond in list(heavy_bonds):
            if excess <= 0:
                break
            bo = bond.GetBondTypeAsDouble()
            if bo >= 2.0:
                new_bo = Chem.BondType.SINGLE if bo >= 3.0 else Chem.BondType.SINGLE
                reduction = int(bo) - 1
                bond.SetBondType(new_bo)
                if verbose:
                    print(
                        f"  NOTE: downgraded bond "
                        f"{atom.GetIdx()}({sym})–"
                        f"{bond.GetOtherAtomIdx(atom.GetIdx())} "
                        f"from order {int(bo)} to 1 to fix excess valence."
                    )
                excess -= reduction
                changed = True

        # 3b: if still overvalent, remove the longest heavy-atom bond(s)
        for bond in heavy_bonds:
            if excess <= 0:
                break
            other_idx = bond.GetOtherAtomIdx(atom.GetIdx())
            rw.RemoveBond(atom.GetIdx(), other_idx)
            if verbose:
                print(
                    f"  NOTE: removed spurious bond "
                    f"{atom.GetIdx()}({sym})–{other_idx} to fix excess valence."
                )
            excess -= 1
            changed = True

    return rw.GetMol() if changed else mol


def _sanitize_and_clean(mol, min_heavy_atoms=5, verbose=True):
    """
    Sanitize an RDKit molecule produced from a CIF and apply cleaning steps:

      1. Partial sanitization to tolerate unusual valences.
      2. Split into disconnected fragments.
      3. Discard known solvents and small ions.
      4. Keep the largest remaining fragment (by heavy-atom count).
      5. Neutralize over-valent atoms (protonated N/O, stray H bonds).
      6. Remove explicit hydrogen atoms.
      7. Full re-sanitization and SMILES round-trip.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
    min_heavy_atoms : int
        Fragments with fewer heavy atoms than this threshold are discarded
        even if they are the only fragment (raises ValueError).
    verbose : bool

    Returns
    -------
    rdkit.Chem.Mol   (sanitized, no explicit H)

    Raises
    ------
    ValueError  on unrecoverable failures.
    """
    if mol is None:
        raise ValueError("Molecule object is None — CIF parsing produced no atoms.")

    # Step 1: partial sanitize (skip valence check so we can fix it ourselves)
    _PARTIAL_FLAGS = (
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS
        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
    )
    try:
        Chem.SanitizeMol(mol, _PARTIAL_FLAGS)
    except Exception:
        pass   # best-effort; errors surface at step 7

    # Step 2: split into fragments (unsanitized to avoid hard failures)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if verbose:
        print(f"  Fragments found : {len(frags)}")

    # Step 3: discard known solvents, small ions, and pure-H fragments
    kept = []
    removed_labels = []
    for frag in frags:
        # Drop fragments made entirely of H atoms
        if frag.GetNumHeavyAtoms() == 0:
            continue
        try:
            frag_noH = Chem.RemoveHs(frag, sanitize=False)
            Chem.SanitizeMol(frag_noH)
            smi = Chem.MolToSmiles(frag_noH)
        except Exception:
            smi = ""
        if smi in _SOLVENT_SMILES:
            removed_labels.append(smi or "?")
            continue
        kept.append(frag)

    if not kept:
        kept = frags      # nothing left — revert to all fragments

    if verbose and removed_labels:
        print(f"  Removed solvents/ions : {', '.join(removed_labels)}")

    # Step 4: keep largest by heavy-atom count
    largest = max(kept, key=lambda m: m.GetNumHeavyAtoms())
    n_heavy = largest.GetNumHeavyAtoms()
    if verbose:
        print(f"  Largest fragment    : {n_heavy} heavy atoms")
    if n_heavy < min_heavy_atoms:
        raise ValueError(
            f"Largest fragment has only {n_heavy} heavy atom(s) "
            f"(threshold: {min_heavy_atoms}).  "
            "Is this really a peptoid CIF?"
        )

    # Step 5: fix over-valent atoms (protonated N/O, stray H-bonds)
    largest = _neutralize_overvalent_atoms(largest, verbose=verbose)

    # Step 6: remove explicit H (suppress RDKit warnings about isolated H),
    # then re-fragment to drop any orphaned bare-H atoms that RemoveHs leaves
    # behind as disconnected components (e.g. when a CIF H–H "bond" exists).
    try:
        mol_noH = Chem.RemoveHs(largest, sanitize=False)
    except Exception:
        mol_noH = largest
    # Re-split and keep largest heavy fragment
    post_frags = Chem.GetMolFrags(mol_noH, asMols=True, sanitizeFrags=False)
    mol_clean = max(post_frags, key=lambda m: m.GetNumHeavyAtoms())

    # Step 7: full sanitize + SMILES round-trip to normalise the graph
    try:
        Chem.SanitizeMol(mol_clean)
    except Exception as e:
        raise ValueError(
            f"Post-cleaning sanitization failed: {e}\n"
            "Tip: the CIF may contain unusual bonding. Try opening it in "
            "Mercury/VESTA to confirm the connectivity is as expected."
        ) from e

    smi = Chem.MolToSmiles(mol_clean, isomericSmiles=True)
    mol_final = Chem.MolFromSmiles(smi)
    if mol_final is None:
        raise ValueError(
            f"SMILES round-trip failed for: {smi}\n"
            "The molecule may contain an unsupported valence or bond type."
        )
    return mol_final


def _filter_disorder(sites, verbose=True):
    """
    Given a list of gemmi SmallStructure.Site objects, return a filtered list
    that contains only one position per disordered atom group.

    CIF structures with multiple rotamers/conformers list the same atom at two
    (or more) positions with partial occupancies and a disorder-group tag
    (e.g. group 1 = conformer A, group 2 = conformer B).  Including all of
    them creates phantom bonds and causes RDKit sanitisation to fail.

    Strategy:
      - Sites with disorder_group == 0  →  always kept (well-ordered).
      - For disordered sites (disorder_group != 0), collect all groups present
        and keep only those belonging to the *lowest* group number (typically 1 /
        conformer A).  The lowest group usually has the highest occupancy in
        SHELX-refined structures.
    """
    disordered = [s for s in sites if s.disorder_group != 0]
    if not disordered:
        return sites          # nothing to do

    groups_present = sorted({s.disorder_group for s in disordered})
    keep_group = groups_present[0]   # e.g. 1 (conformer A)

    if verbose:
        print(
            f"  Disorder groups found : {groups_present}  "
            f"→ keeping group {keep_group} (conformer A)"
        )

    filtered = [s for s in sites if s.disorder_group == 0 or s.disorder_group == keep_group]
    if verbose:
        n_removed = len(sites) - len(filtered)
        print(f"  Removed {n_removed} alternate-conformer site(s)")
    return filtered


def cif_to_mol(cif_path, verbose=True):
    """
    Read a small-molecule CIF file and return a sanitized RDKit molecule
    ready for peptoid decomposition.

    Pipeline
    --------
    1. Parse CIF with gemmi (SmallStructure).
    2. Filter disordered atoms: keep only the primary conformer (group 1 /
       conformer A) to avoid phantom bonds from alternate rotamers.
    3. If _geom_bond data is present → use explicit bonds + distance-based
       bond-order assignment.
    4. Otherwise → use 3-D coordinates + RDKit rdDetermineBonds (with
       distance-heuristic upgrade for double/triple bonds).
    5. Clean up: remove counterions/solvents, keep largest fragment,
       strip explicit H, sanitize, SMILES round-trip.

    Parameters
    ----------
    cif_path : str
        Path to the CIF file.
    verbose : bool
        Print progress messages.

    Returns
    -------
    mol : rdkit.Chem.Mol
        Sanitized RDKit molecule (no explicit H).
    smiles : str
        Canonical SMILES of the cleaned molecule.

    Raises
    ------
    ImportError  if gemmi is not installed.
    ValueError   if the CIF cannot be parsed or yields no valid molecule.
    """
    if not _GEMMI_AVAILABLE:
        raise ImportError(
            "gemmi is required for CIF ingestion.\n"
            "Install it with:  pip install gemmi"
        )

    if verbose:
        print(f"\nReading CIF : {cif_path}")

    # ---- Parse with gemmi ------------------------------------------------
    try:
        doc = gemmi.cif.read(cif_path)
    except Exception as e:
        raise ValueError(f"gemmi could not read '{cif_path}': {e}") from e

    block = doc.sole_block()
    try:
        ss = gemmi.make_small_structure_from_block(block)
    except Exception as e:
        raise ValueError(
            f"gemmi could not interpret '{cif_path}' as a small-molecule "
            f"structure: {e}"
        ) from e

    n_sites = len(ss.sites)
    if n_sites == 0:
        raise ValueError(
            "No atomic sites found in the CIF.\n"
            "Ensure the file contains _atom_site_fract_* coordinates."
        )
    if verbose:
        print(f"  Asymmetric-unit sites : {n_sites}")

    # ---- Filter alternate conformers / disordered rotamers ---------------
    filtered_sites = _filter_disorder(list(ss.sites), verbose=verbose)

    # ---- Strategy 1: explicit _geom_bond data ----------------------------
    mol = None
    strategy_used = None

    mol_geom = _build_mol_from_geom_bonds(block, ss, sites=filtered_sites)
    if mol_geom is not None and mol_geom.GetNumAtoms() > 0:
        if verbose:
            print("  Bond strategy : _geom_bond loop (explicit bonds + distance orders)")
        mol = mol_geom
        strategy_used = "geom_bond"
    else:
        # ---- Strategy 2: coordinate-only ---------------------------------
        if verbose:
            print("  Bond strategy : coordinate-only (rdDetermineBonds + distance heuristic)")
        mol = _build_mol_from_coords(ss, sites=filtered_sites)
        strategy_used = "coords"

    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(
            "Could not construct a molecule from the CIF coordinates.\n"
            "Please verify the file contains valid fractional coordinates."
        )

    # ---- Clean and sanitize ----------------------------------------------
    if verbose:
        print("  Cleaning molecule (removing counterions / solvents) …")
    mol_clean = _sanitize_and_clean(mol, verbose=verbose)

    smiles = Chem.MolToSmiles(mol_clean, isomericSmiles=True)
    if verbose:
        print(f"  Extracted SMILES : {smiles}")
        print(f"  Heavy atoms      : {mol_clean.GetNumHeavyAtoms()}")

    return mol_clean, smiles


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_structure_image(mol, path, size=(1200, 700)):
    """Save a 2D depiction of the complete molecule."""
    mol2d = copy.deepcopy(mol)
    Compute2DCoords(mol2d)
    MolToFile(mol2d, path, size=size)


def save_residue_grid(residues, path, cols=4, sub_size=(420, 320)):
    """
    Draw each residue as *N(sidechain)CC(=O)* in a labelled grid image.
    Falls back to sidechain-only depiction if the residue SMILES cannot be parsed.
    """
    mols, legends = [], []
    for r in residues:
        rs = r["residue_smiles"].replace("*", "[*]")
        m = Chem.MolFromSmiles(rs)
        if m is None:
            m = Chem.MolFromSmiles(r["sidechain_smiles"])
        if m:
            Compute2DCoords(m)
            mols.append(m)
            legends.append(f"Residue {r['index']}\n{r['residue_smiles']}")

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
        # Fallback: distance-geometry with random coords
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
            "Decompose a peptoid (linear or macrocyclic) into its residues.\n"
            "Input can be a SMILES string (-s) or a small-molecule CIF file (-c).\n\n"
            "CIF files are cleaned automatically:\n"
            "  - counterions and solvent molecules are removed\n"
            "  - the largest fragment is kept\n"
            "  - bond orders are inferred from _geom_bond data or coordinates\n\n"
            "Three output files are written for every run:\n"
            "  <prefix>_structure.png  - 2D whole-molecule image\n"
            "  <prefix>_residues.png   - 2D residue grid image\n"
            "  <prefix>_3d.sdf         - 3D conformer (SDF with H atoms)"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-s", "--smiles",
        type=str,
        metavar="SMILES",
        help="Full-structure SMILES of the peptoid.",
    )
    input_group.add_argument(
        "-c", "--cif",
        type=str,
        metavar="CIF_FILE",
        help=(
            "Path to a small-molecule CIF file.\n"
            "The pipeline removes counterions/solvents, builds a sanitized\n"
            "molecule, and extracts the SMILES automatically."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="peptoid",
        metavar="PREFIX",
        help="Output filename prefix (default: 'peptoid').",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Obtain the molecule
    # ------------------------------------------------------------------
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"ERROR: could not parse SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        input_smiles = args.smiles
        print(f"\nSMILES     : {input_smiles}")

    else:  # args.cif
        try:
            mol, input_smiles = cif_to_mol(args.cif, verbose=True)
        except (ImportError, ValueError) as e:
            print(f"\nERROR: {e}", file=sys.stderr)
            sys.exit(1)

    print(f"Atom count : {mol.GetNumAtoms()} heavy atoms")

    # ------------------------------------------------------------------
    # Extract residues
    # ------------------------------------------------------------------
    try:
        residues, topology = extract_residues(mol)
    except ValueError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Topology   : {topology}")
    print(f"\nResidues ({len(residues)} found):")
    for r in residues:
        print(f"  [{r['index']:2d}]  {r['residue_smiles']}")
        print(f"        sidechain: {r['sidechain_smiles']}")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------

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
