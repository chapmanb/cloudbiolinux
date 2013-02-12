from fabric.api import sudo
from fabric.contrib.files import exists

def _read_boolean(env, name, default):
    property_str = env.get(name, str(default))
    return property_str.upper() in ["TRUE", "YES"]


def _chown_galaxy(env, path):
    """
    Recursively change ownership of ``path``, first checking if ``path`` exists.
    """
    chown_command = "chown --recursive %s:%s %s"
    galaxy_user = env.get("galaxy_user", "galaxy")
    if exists(path):
        sudo(chown_command % (galaxy_user, galaxy_user, path))

def _dir_is_empty(path):
    """
    Return ``True`` is ``path`` directory has no files or folders in it.
    Return ``False`` otherwise.
    """
    if "empty" in sudo('[ "$(ls -A {0})" ] || echo "empty"'.format(path)):
        return True
    return False
