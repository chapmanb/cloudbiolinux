from fabric.api import sudo


def _read_boolean(env, name, default):
    property_str = env.get(name, str(default))
    return property_str.upper() in ["TRUE", "YES"]


def _chown_galaxy(env, path):
    chown_command = "chown --recursive %s:%s %s"
    galaxy_user = env.get("galaxy_user", "galaxy")
    sudo(chown_command % (galaxy_user, galaxy_user, path))
