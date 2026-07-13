import os
basedir = os.path.abspath(os.path.dirname(__file__))

# SQL Alchemy configure from object


class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'app.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    PEPTOIDS_PER_PAGE = 6

    # Contributor/import hardening settings.
    # These keep public browsing online while limiting heavier SMILES/CIF work.
    CONTRIBUTIONS_ENABLED = os.environ.get('CONTRIBUTIONS_ENABLED', 'true').lower() in ('1', 'true', 'yes', 'on')
    MAX_CONTENT_LENGTH = int(os.environ.get('MAX_CONTENT_LENGTH', 5 * 1024 * 1024))
    MAX_SMILES_LENGTH = int(os.environ.get('MAX_SMILES_LENGTH', 10000))
    RATELIMIT_STORAGE_URI = os.environ.get('RATELIMIT_STORAGE_URI', 'memory://')

    BASIC_AUTH_USERNAME = 'kklab'
    BASIC_AUTH_PASSWORD = 'peptoids1234!'
