from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, SubmitField, RadioField, TextAreaField, SelectField, DateField
from wtforms.validators import DataRequired, Optional


class SearchForm(FlaskForm):
    residue = StringField(
        'Residue',
        render_kw={"placeholder": "(S)-N-(1-cyclohexylethyl)glycine or Nsch"}
    )
    author = StringField(
        'Author',
        render_kw={"placeholder": "Kirshenbaum, Kent"}
    )
    doi = StringField(
        'DOI',
        render_kw={"placeholder": "10.1002/ejoc.202001401"}
    )
    experiment = RadioField(
        'Experimental Technique',
        choices=[
            ('', 'No experiment filter'),
            ('X-Ray Diffraction', 'X-Ray Diffraction'),
            ('Solution NMR', 'Solution NMR')
        ],
        default=''
    )
    topology = RadioField(
        'Topology',
        choices=[
            ('', 'No topology filter'),
            ('A', 'Linear'),
            ('C', 'Cyclic'),
            ('M', 'Multicyclic')
        ],
        default=''
    )
    submit = SubmitField('Submit')
class ImportPeptoidForm(FlaskForm):
    code = StringField('Peptoid code', validators=[DataRequired()])
    title = StringField('Title', validators=[DataRequired()])
    release = DateField('Release date', validators=[DataRequired()])

    experiment = SelectField(
        'Experimental technique',
        choices=[
            ('X-Ray Diffraction', 'X-Ray Diffraction'),
            ('Solution NMR', 'Solution NMR'),
            ('Other', 'Other')
        ],
        validators=[DataRequired()]
    )

    pub_doi = StringField('Publication DOI', validators=[Optional()])
    struct_doi = StringField('Structure DOI', validators=[Optional()])
    citation = TextAreaField('Citation', validators=[Optional()])

    authors = TextAreaField(
        'Authors',
        validators=[DataRequired()],
        description='Enter one author per line as First name, Last name.'
    )

    smiles = TextAreaField(
        'Full-structure SMILES',
        validators=[Optional()]
    )

    cif_file = FileField(
        'CIF file',
        validators=[
            Optional(),
            FileAllowed(['cif'], 'Please upload a .cif file.')
        ]
    )

    import_submit = SubmitField('Process and preview')

    def validate(self, extra_validators=None):
        valid = super().validate(extra_validators=extra_validators)
        has_smiles = bool(self.smiles.data and self.smiles.data.strip())
        has_cif = bool(self.cif_file.data and self.cif_file.data.filename)

        if has_smiles == has_cif:
            message = (
                'Provide either one CIF file or one SMILES string, but not both.'
            )
            self.smiles.errors.append(message)
            self.cif_file.errors.append(message)
            return False

        return valid
