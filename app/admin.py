from flask_admin import Admin, AdminIndexView, expose
from flask_admin.contrib.sqla import ModelView
from app.models import Peptoid, Author, Residue, Contributor, Submission
from app import app, db, basic_auth
from flask_admin.contrib.fileadmin import FileAdmin
import os.path as op
from flask import redirect, Response, abort, send_file, flash, url_for
from werkzeug.exceptions import HTTPException
from wtforms import PasswordField, TextAreaField, SubmitField
from flask_wtf import FlaskForm
from wtforms.validators import DataRequired
from markupsafe import Markup
from app.services.finalization import finalize_submission, FinalizationError
import json
import os
from wtforms.validators import ValidationError

#class used to force admin to enter credentials
class AuthException(HTTPException):
    def __init__(self, message):
        super().__init__(message, Response(
            "You could not be authenticated. Please refresh the page.", 401,
            {'WWW-Authenticate': 'Basic realm="Login Required"'}
        ))


#inherits from AdminIndexView function that will make admin views only accesible if the user is authenticated with BasicAuth
class MyAdminIndexView(AdminIndexView):
    def is_accessible(self):
        if not basic_auth.authenticate():
            raise AuthException('Not authenticated.')
        else:
            return True

    def inaccessible_callback(self, name, **kwargs):
        return redirect(basic_auth.challenge())
        
#Creates Admin with custom index view and adds model views for database models
admin = Admin(app, name='PeptoidDB Admin', template_mode='bootstrap4',index_view=MyAdminIndexView(template='admin/home.html'))
admin.add_view(ModelView(Peptoid, db.session))
admin.add_view(ModelView(Residue, db.session))
admin.add_view(ModelView(Author, db.session))

class ContributorAdmin(ModelView):
    column_list = (
        'username', 'email', 'first_name', 'last_name',
        'institution', 'is_active', 'created_at', 'last_login'
    )
    form_columns = (
        'username', 'email', 'password', 'first_name',
        'last_name', 'institution', 'is_active'
    )
    form_extra_fields = {
        'password': PasswordField(
            'Password',
            description='Required for new accounts; leave blank to keep the current password.'
        )
    }

    def on_model_change(self, form, model, is_created):
        password = form.password.data
        if is_created and not password:
            raise ValidationError('A password is required for new contributors.')
        if password:
            model.set_password(password)

admin.add_view(ContributorAdmin(Contributor, db.session))


class ApproveSubmissionForm(FlaskForm):
    submit = SubmitField('Approve and publish')


class RejectSubmissionForm(FlaskForm):
    reason = TextAreaField('Rejection reason', validators=[DataRequired()])
    submit = SubmitField('Reject submission')


class SubmissionAdmin(ModelView):
    def is_accessible(self):
        if not basic_auth.authenticate():
            raise AuthException('Not authenticated.')
        return True

    def inaccessible_callback(self, name, **kwargs):
        return redirect(basic_auth.challenge())
    can_create = False
    can_delete = False
    can_view_details = True

    column_list = (
        'id', 'contributor', 'status', 'proposed_code',
        'title', 'experiment', 'topology', 'created_at', 'updated_at'
    )
    column_filters = ('status', 'experiment', 'topology', 'created_at')
    column_searchable_list = (
        'proposed_code', 'title', 'authors',
        'pub_doi', 'struct_doi'
    )
    column_default_sort = ('updated_at', True)
    column_labels = {'proposed_code': 'Database code'}
    column_formatters = {
        'id': lambda view, context, model, name: Markup(
            '<a href="{}">Review #{}</a>'.format(
                url_for('.review_view', submission_id=model.id), model.id
            )
        )
    }

    column_details_list = (
        'id', 'contributor', 'status', 'created_at', 'updated_at',
        'proposed_code', 'finalized_code', 'title', 'release',
        'experiment', 'pub_doi', 'struct_doi', 'citation', 'authors',
        'input_type', 'original_smiles', 'cleaned_smiles', 'topology',
        'sequence', 'residue_data_json', 'warnings_json',
        'original_cif_filename', 'rejection_reason'
    )

    form_columns = (
        'title', 'release', 'experiment',
        'pub_doi', 'struct_doi', 'citation', 'authors',
        'sequence', 'rejection_reason'
    )

    @expose('/review/<int:submission_id>')
    def review_view(self, submission_id):
        submission = db.session.get(Submission, submission_id)
        if submission is None:
            abort(404)
        return self.render(
            'admin/submission_review.html',
            submission=submission,
            residues=json.loads(submission.residue_data_json),
            reject_form=RejectSubmissionForm(),
            approve_form=ApproveSubmissionForm(),
        )

    @expose('/review/<int:submission_id>/image/<kind>')
    def review_image(self, submission_id, kind):
        submission = db.session.get(Submission, submission_id)
        if submission is None:
            abort(404)
        path = {
            'structure': submission.structure_image_path,
            'residues': submission.residue_image_path,
        }.get(kind)
        if not path or not os.path.isfile(path):
            abort(404)
        return send_file(path, mimetype='image/png')

    @expose('/review/<int:submission_id>/approve', methods=['POST'])
    def approve_view(self, submission_id):
        submission = db.session.get(Submission, submission_id)
        if submission is None:
            abort(404)
        form = ApproveSubmissionForm()
        if not form.validate_on_submit():
            abort(400)
        try:
            peptoid, backup_dir = finalize_submission(submission)
        except FinalizationError as error:
            flash(str(error), 'danger')
            return redirect(url_for('.review_view', submission_id=submission.id))
        except Exception as error:
            flash('Approval failed and was rolled back: {}'.format(error), 'danger')
            return redirect(url_for('.review_view', submission_id=submission.id))
        flash(
            'Submission approved as {}. Backup: {}'.format(
                peptoid.code, backup_dir
            ),
            'success',
        )
        return redirect(url_for('.review_view', submission_id=submission.id))

    @expose('/review/<int:submission_id>/reject', methods=['POST'])
    def reject_view(self, submission_id):
        submission = db.session.get(Submission, submission_id)
        if submission is None:
            abort(404)
        form = RejectSubmissionForm()
        if not form.validate_on_submit():
            abort(400)
        if submission.status != 'pending':
            flash('Only pending submissions may be rejected.', 'warning')
            return redirect(url_for('.review_view', submission_id=submission.id))
        submission.status = 'rejected'
        submission.proposed_code = None
        submission.rejection_reason = form.reason.data.strip()
        db.session.commit()
        flash('Submission rejected.', 'success')
        return redirect(url_for('.review_view', submission_id=submission.id))

admin.add_view(
    SubmissionAdmin(Submission, db.session, name='Submissions')
)

#Views for image uploads of peptoid structures and residues
class PeptoidImageAdmin(FileAdmin):
    can_upload=True
    can_delete=True
    can_delete_dirs=False
    can_mkdir=False
    can_rename=True
    allowed_extensions=['png']
    editable_extensions=[]

class ResidueImageAdmin(FileAdmin):
    can_upload=True
    can_delete=True
    can_delete_dirs=False
    can_mkdir=False
    can_rename=True
    allowed_extensions=['png']
    editable_extensions=[]

peptoid_image_path = op.join(op.dirname(__file__), 'static/peptoids')
admin.add_view(PeptoidImageAdmin(base_path=peptoid_image_path, name='Peptoid Images'))

residue_image_path = op.join(op.dirname(__file__), 'static/residues')
admin.add_view(ResidueImageAdmin(base_path=residue_image_path, name='Residue Images'))