# from the app module imports db instance of SQLAlchemy application
from enum import unique
from app import db
import datetime
from flask import url_for
from sqlalchemy import event
from flask_login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash

# peptoid-author helper table
peptoid_author = db.Table('peptoid-author',
                          db.Column('peptoid_id', db.String(16), db.ForeignKey(
                              'peptoid.code'), unique=False),
                          db.Column('author_id', db.Integer, db.ForeignKey(
                              'author.id'), unique=False)
                          )

# peptoid-residue helper table
peptoid_residue = db.Table('peptoid-residue',
                           db.Column('peptoid_id', db.String(16), db.ForeignKey(
                               'peptoid.code'), unique=False),
                           db.Column('residue_id', db.Integer, db.ForeignKey(
                               'residue.id'), unique=False)
                           )

# peptoid table: image file name, title to display on page, data base code, release date,
# experimental technique, doi of publication, relationship with the author and residue


class Peptoid(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    code = db.Column(db.String(16), unique=True)
    title = db.Column(db.Text, unique=False)
    release = db.Column(db.DateTime, unique=False)
    experiment = db.Column(db.Text, unique=False)
    pub_doi = db.Column(db.String(64), unique=False)
    struct_doi = db.Column(db.String(64), unique=False)
    citation = db.Column(db.String(1024),unique=False)
    topology = db.Column(db.String(1), unique=False)
    sequence = db.Column(db.String(1024), unique=False)
    n_term = db.Column(db.String(32), unique=False)
    c_term = db.Column(db.String(32), unique=False)
    cyclization_points = db.Column(db.String(8), unique=False)
    struct_smiles = db.Column(db.String(128),unique=False)

    peptoid_author = db.relationship('Author', secondary=peptoid_author, lazy='dynamic',
                                     backref=db.backref('peptoids', order_by='Peptoid.release.desc()'))
    peptoid_residue = db.relationship('Residue', secondary=peptoid_residue, lazy='dynamic',
                                      backref=db.backref('peptoids', order_by='Peptoid.release.desc()'))

    def to_dict(self):
        data = {
            'title': self.title,
            'release': self.release,
            'experiment': self.experiment,
            'pub_doi': self.pub_doi,
            'struct_doi': self.struct_doi,
            'citation': self.citation,
            'topology': self.topology,
            'sequence':self.sequence,
            'n_term':self.n_term,
            'c_term':self.c_term,
            'cyc_points':self.cyclization_points,
            'struct_smiles':self.struct_smiles,
            '_links': {
                'self': url_for('api.get_peptoid', code=self.code),
                'residues': url_for('api.get_pep_residues', code=self.code),
                'authors': url_for('api.get_pep_authors', code=self.code),
                '3d_image': url_for('static', filename='peptoids/'+self.code[:5]+'.png'),
                '2d_image': url_for('static', filename='peptoids/'+self.code+'_2d.png')
            }
        }
        return data

    def __repr__(self):
        return '<Peptoid {}>'.format(self.title)

# authors table: first name and last name
class Author(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    first_name = db.Column(db.Text, unique=False)
    last_name = db.Column(db.Text, unique=False)

    def to_dict(self):
        data = {
            'first_name': self.first_name,
            'last_name': self.last_name,
            'peptoids': {}
        }
        for p in self.peptoids:
            data['peptoids'][p.code] = url_for('api.get_peptoid', code=p.code)
        return data

    def __repr__(self):
        return '<Author {}>'.format(self.last_name + ", " + self.first_name)

# residues table: nomenclature


class Residue(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    long_name = db.Column(db.Text, unique=False)
    short_name = db.Column(db.Text, unique=False)
    pep_type = db.Column(db.Text, unique=False)
    monomer_structure = db.Column(db.Text, unique=False)
    SMILES = db.Column(db.Text, unique=False)

    def to_dict(self):
        data = {
            'full_nomenclature': self.long_name,
            'short_name': self.short_name,
            'type': self.pep_type,
            'monomer_structure': self.monomer_structure,
            'SMILES': self.SMILES,
            'peptoids': {}
        }
        for p in self.peptoids:
            data['peptoids'][p.code] = url_for('api.get_peptoid', code=p.code)
        return data

    def __repr__(self):
        return '<Residue {}>'.format(self.long_name)


class Contributor(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), unique=True, nullable=False, index=True)
    email = db.Column(db.String(255), unique=True, nullable=False, index=True)
    password_hash = db.Column(db.String(255), nullable=False)
    first_name = db.Column(db.String(100), nullable=False)
    last_name = db.Column(db.String(100), nullable=False)
    institution = db.Column(db.String(255), nullable=False)
    is_active = db.Column(db.Boolean, nullable=False, default=True)
    created_at = db.Column(
        db.DateTime, nullable=False, default=datetime.datetime.utcnow
    )
    last_login = db.Column(db.DateTime, nullable=True)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def __repr__(self):
        return '<Contributor {}>'.format(self.username)


class Submission(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    contributor_id = db.Column(
        db.Integer, db.ForeignKey('contributor.id'), nullable=False, index=True
    )
    status = db.Column(db.String(20), nullable=False, default='draft', index=True)
    created_at = db.Column(
        db.DateTime, nullable=False, default=datetime.datetime.utcnow
    )
    updated_at = db.Column(
        db.DateTime, nullable=False, default=datetime.datetime.utcnow,
        onupdate=datetime.datetime.utcnow
    )

    proposed_code = db.Column(db.String(16), nullable=True)
    finalized_code = db.Column(db.String(16), nullable=True)
    title = db.Column(db.Text, nullable=False)
    release = db.Column(db.DateTime, nullable=False)
    experiment = db.Column(db.Text, nullable=False)
    pub_doi = db.Column(db.String(64), nullable=True)
    struct_doi = db.Column(db.String(64), nullable=True)
    citation = db.Column(db.String(1024), nullable=True)
    authors = db.Column(db.Text, nullable=False)

    input_type = db.Column(db.String(16), nullable=False)
    original_smiles = db.Column(db.Text, nullable=True)
    cleaned_smiles = db.Column(db.Text, nullable=False)
    topology = db.Column(db.String(20), nullable=False)
    sequence = db.Column(db.String(1024), nullable=True)
    residue_data_json = db.Column(db.Text, nullable=False)
    warnings_json = db.Column(db.Text, nullable=True)

    structure_image_path = db.Column(db.Text, nullable=True)
    residue_image_path = db.Column(db.Text, nullable=True)
    staged_cif_path = db.Column(db.Text, nullable=True)
    original_cif_filename = db.Column(db.Text, nullable=True)
    rejection_reason = db.Column(db.Text, nullable=True)

    contributor = db.relationship(
        'Contributor',
        backref=db.backref('submissions', lazy='dynamic')
    )

    def __repr__(self):
        return '<Submission {} {}>'.format(self.id, self.status)
