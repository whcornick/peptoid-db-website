import json
import os
import shutil
import sqlite3
from datetime import datetime

from rdkit import Chem

from app import app, db
from app.models import Author, Peptoid, Residue
from app.services.code_generation import CODE_PATTERN


class FinalizationError(Exception):
    pass


def _canonical_residue_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles.replace('*', '[*]'))
    if molecule is None:
        raise FinalizationError(
            'Could not interpret residue SMILES: {}'.format(smiles)
        )
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def _parse_authors(author_text):
    authors = []
    for line in author_text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [part.strip() for part in line.split(',', 1)]
        if len(parts) != 2 or not all(parts):
            raise FinalizationError(
                'Each author must use First name, Last name format.'
            )
        authors.append((parts[0], parts[1]))

    if not authors:
        raise FinalizationError('At least one author is required.')

    return authors


def _topology_code(topology):
    values = {
        'linear': 'A',
        'cyclic': 'C',
        'multicyclic': 'M',
        'A': 'A',
        'C': 'C',
        'M': 'M',
    }
    try:
        return values[topology]
    except KeyError:
        raise FinalizationError(
            'Unsupported topology: {}'.format(topology)
        )


def _database_path():
    prefix = 'sqlite:///'
    uri = app.config['SQLALCHEMY_DATABASE_URI']
    if not uri.startswith(prefix):
        raise FinalizationError(
            'Automatic backups currently require SQLite.'
        )
    return uri[len(prefix):]


def _create_backup():
    backup_dir = os.path.join(
        app.instance_path,
        'finalization_backups',
        datetime.utcnow().strftime('%Y%m%dT%H%M%S%f'),
    )
    os.makedirs(backup_dir, exist_ok=False)

    source = sqlite3.connect(_database_path())
    destination = sqlite3.connect(
        os.path.join(backup_dir, 'database.db')
    )
    try:
        source.backup(destination)
    finally:
        destination.close()
        source.close()

    image_source = os.path.join(app.static_folder, 'peptoids')
    if os.path.isdir(image_source):
        shutil.copytree(
            image_source,
            os.path.join(backup_dir, 'peptoid_images'),
        )

    return backup_dir


def finalize_submission(submission):
    if submission.status != 'pending':
        raise FinalizationError(
            'Only pending submissions may be approved.'
        )

    required = {
        'title': submission.title,
        'release': submission.release,
        'experiment': submission.experiment,
        'authors': submission.authors,
        'cleaned SMILES': submission.cleaned_smiles,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise FinalizationError(
            'Missing required metadata: {}'.format(', '.join(missing))
        )

    if len(submission.cleaned_smiles) > 128:
        raise FinalizationError(
            'Cleaned SMILES exceeds the finalized database limit.'
        )

    author_names = _parse_authors(submission.authors)
    residue_data = json.loads(submission.residue_data_json)
    if not residue_data:
        raise FinalizationError('No residues were detected.')

    code = (submission.proposed_code or '').strip()
    match = CODE_PATTERN.fullmatch(code)
    if not match:
        raise FinalizationError(
            'The submission does not have a valid reserved database code.'
        )

    code_parts = match.groupdict()
    if int(code_parts['length']) != len(residue_data):
        raise FinalizationError(
            'The reserved code does not match the detected residue count.'
        )
    if code_parts['topology'] != _topology_code(submission.topology):
        raise FinalizationError(
            'The reserved code does not match the detected topology.'
        )

    if Peptoid.query.filter_by(code=code).first():
        raise FinalizationError(
            'Generated peptoid code {} already exists.'.format(code)
        )

    existing_residues = {}
    for residue in Residue.query.all():
        key = _canonical_residue_smiles(residue.SMILES)
        existing_residues[key] = residue

    matched_residues = []
    sequence_names = []
    for generated in residue_data:
        key = _canonical_residue_smiles(
            generated['residue_smiles']
        )
        residue = existing_residues.get(key)
        if residue is None:
            raise FinalizationError(
                'No approved residue record matches residue {}: {}'
                .format(generated['index'], generated['residue_smiles'])
            )
        sequence_names.append(residue.long_name)
        if residue not in matched_residues:
            matched_residues.append(residue)

    if not submission.structure_image_path or not os.path.isfile(
        submission.structure_image_path
    ):
        raise FinalizationError('The generated structure image is missing.')

    backup_dir = _create_backup()
    output_path = os.path.join(
        app.static_folder, 'peptoids', '{}_2d.png'.format(code)
    )
    temp_output = output_path + '.tmp'
    image_created = False

    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        shutil.copy2(submission.structure_image_path, temp_output)

        peptoid = Peptoid(
            code=code,
            title=submission.title,
            release=submission.release,
            experiment=submission.experiment,
            pub_doi=submission.pub_doi,
            struct_doi=submission.struct_doi,
            citation=submission.citation,
            topology=_topology_code(submission.topology),
            sequence=', '.join(sequence_names),
            struct_smiles=submission.cleaned_smiles,
        )

        for first_name, last_name in author_names:
            author = Author.query.filter_by(
                first_name=first_name,
                last_name=last_name,
            ).order_by(Author.id.asc()).first()
            if author is None:
                author = Author(
                    first_name=first_name,
                    last_name=last_name,
                )
                db.session.add(author)
            peptoid.peptoid_author.append(author)

        for residue in matched_residues:
            peptoid.peptoid_residue.append(residue)

        db.session.add(peptoid)
        submission.status = 'approved'
        submission.finalized_code = code
        db.session.flush()

        os.replace(temp_output, output_path)
        image_created = True
        db.session.commit()
        return peptoid, backup_dir

    except Exception:
        db.session.rollback()
        if os.path.exists(temp_output):
            os.remove(temp_output)
        if image_created and os.path.exists(output_path):
            os.remove(output_path)
        raise
