#importing essential modules for instantiating application and extensions
from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_bootstrap import Bootstrap
from flask import Blueprint
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_graphql import GraphQLView
from flask_moment import Moment
from flask_admin import Admin, BaseView, expose
from flask_admin.contrib.sqla import ModelView
from flask_basicauth import BasicAuth
from flask_login import LoginManager

#instantiating app, database, config, limiter, bootstrap, moment, and basic_auth
app = Flask(__name__)
app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db, render_as_batch=True)
bootstrap = Bootstrap(app)
limiter = Limiter(key_func=get_remote_address, app=app, storage_uri=app.config.get('RATELIMIT_STORAGE_URI'))
moment = Moment(app)
basic_auth = BasicAuth(app)
login_manager = LoginManager(app)
login_manager.login_view = 'routes.contributor_login'
login_manager.login_message_category = 'info'

#setting up admin
from app.admin import admin

#registering blueprints
from app.api import bp as api_bp
app.register_blueprint(api_bp, url_prefix='/api')

from app.errors import bp as errors_bp
app.register_blueprint(errors_bp)

from app.routes import bp as routes_bp
app.register_blueprint(routes_bp)

#importing models from app module
from app import models

@login_manager.user_loader
def load_contributor(contributor_id):
    return db.session.get(models.Contributor, int(contributor_id))
from .schema import schema

class myGraphQLView(GraphQLView):
    decorators = [limiter.limit('1000 per minute')]
    def dispatch_request(self):
        return GraphQLView.dispatch_request(self)

app.add_url_rule('/graphql', view_func=myGraphQLView.as_view('graphql', schema=schema, graphiql=True))
