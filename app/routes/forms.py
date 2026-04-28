from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, RadioField


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
