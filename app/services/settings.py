from app import app, db
from app.models import AppSetting

CONTRIBUTIONS_PAUSED_KEY = 'contributions_paused'


def get_setting(key, default=None):
    setting = db.session.get(AppSetting, key)
    if setting is None:
        return default
    return setting.value


def set_setting(key, value, description=None):
    setting = db.session.get(AppSetting, key)
    if setting is None:
        setting = AppSetting(key=key, value=str(value), description=description)
        db.session.add(setting)
    else:
        setting.value = str(value)
        if description is not None:
            setting.description = description
    db.session.commit()
    return setting


def _truthy(value):
    return str(value).strip().lower() in ('1', 'true', 'yes', 'on')


def contributions_are_paused():
    return _truthy(get_setting(CONTRIBUTIONS_PAUSED_KEY, 'false'))


def contributions_are_enabled():
    if not app.config.get('CONTRIBUTIONS_ENABLED', True):
        return False
    return not contributions_are_paused()


def set_contributions_paused(paused):
    return set_setting(
        CONTRIBUTIONS_PAUSED_KEY,
        'true' if paused else 'false',
        'Controls whether contributor submissions are temporarily paused.',
    )
