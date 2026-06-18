from flask_admin import Admin, AdminIndexView
from flask_admin.contrib.sqla import ModelView
from app.models import Peptoid, Author, Residue, Contributor
from app import app, db, basic_auth
from flask_admin.contrib.fileadmin import FileAdmin
import os.path as op
from flask import redirect, Response
from werkzeug.exceptions import HTTPException
from wtforms import PasswordField
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
admin = Admin(app, name='PeptoidDB Admin', template_mode='bootstrap4',index_view=MyAdminIndexView())
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