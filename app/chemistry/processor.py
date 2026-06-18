import base64
import os
import tempfile

from rdkit import Chem

from app.chemistry import decompose
from app.chemistry.decompose_cif import cif_to_mol


def _encode_file(path):
    with open(path, 'rb') as file:
        return base64.b64encode(file.read()).decode('ascii')


def process_structure(smiles=None, cif_path=None):
    if bool(smiles) == bool(cif_path):
        raise ValueError('Provide exactly one SMILES string or CIF file.')

    if cif_path:
        mol, cleaned_smiles = cif_to_mol(cif_path, verbose=False)
        input_type = 'CIF'
    else:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            raise ValueError('RDKit could not parse the supplied SMILES.')
        cleaned_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        input_type = 'SMILES'

    residues, topology = decompose.extract_residues(mol)

    with tempfile.TemporaryDirectory() as temp_dir:
        structure_path = os.path.join(temp_dir, 'structure.png')
        residues_path = os.path.join(temp_dir, 'residues.png')

        decompose.save_structure_image(mol, structure_path)
        decompose.save_residue_grid(residues, residues_path)

        structure_image = _encode_file(structure_path)
        residue_image = _encode_file(residues_path)

    return {
        'input_type': input_type,
        'cleaned_smiles': cleaned_smiles,
        'topology': topology,
        'residue_count': len(residues),
        'residues': residues,
        'structure_image': structure_image,
        'residue_image': residue_image,
    }
