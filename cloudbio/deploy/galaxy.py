"""Fabric (http://docs.fabfile.org) deployment file to set up galaxy.
"""

import os
import time

from cloudbio.custom.galaxy import install_galaxy_webapp

from fabric.api import sudo, run, env, cd
from fabric.contrib.files import append

from cloudbio.galaxy.tools import _setup_install_dir
from cloudbio.galaxy.utils import _chown_galaxy
from util import start_service

## TODO: Investigate whether it would make sense to move more of this
## into cloudbio.galaxy (or maybe cloudbio.custom.galaxy)


def wait_for_galaxy():

    while not "8080" in run("netstat -lant"):
        # Check if galaxy has started
        print "Waiting for galaxy to start."
        time.sleep(10)


def purge_galaxy():
    sudo("/etc/init.d/galaxy stop")
    sudo("rm -rf %s" % env.galaxy_home)
    init_script = "postgresql"
    # if env.postgres_version[0] < '9':
    #    # Postgres 8.4 had different name for script
    #    init_script = "postgresql-%s" % env.postgres_version
    sudo("/etc/init.d/%s restart" % init_script)
    sudo('psql  -c "drop database galaxy;"', user="postgres")
    sudo('psql  -c "create database galaxy;"', user="postgres")


def setup_galaxy(options, seed=True):
    """Deploy a Galaxy server along with some tools.
    """
    _setup_install_dir(env)  # Still needed? -John
    install_galaxy_webapp(env)
    #_fix_galaxy_permissions()
    _setup_shed_tools_dir()
    _setup_galaxy_log_dir()
    _migrate_galaxy_database()
    if seed:
        seed_database(options["galaxy"])
    _start_galaxy()


def _migrate_galaxy_database():
    with cd(env.galaxy_home):
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; python ./scripts/build_universe_config.py conf.d; python -ES ./scripts/fetch_eggs.py; ./create_db.sh'", user="galaxy")


def seed_database(galaxy_data):
    with cd(env.galaxy_home):
        sudo("rm -f seed.py")
        _setup_database_seed_file(galaxy_data)
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; python ./scripts/build_universe_config.py conf.d; python -ES ./scripts/fetch_eggs.py; python seed.py'", user="galaxy")


def seed_workflows(options):
    wait_for_galaxy()
    galaxy_data = options["galaxy"]
    with cd(os.path.join(env.galaxy_home, "workflows")):
        for user in galaxy_data["users"]:
            api_key = user["api_key"]
            workflows = None
            if "workflows" in user:
                workflows = user["workflows"]
            if not workflows:
                continue
            for workflow in workflows:
                sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; bash import_all.sh %s %s'" % (api_key, workflow), user=env.galaxy_user)


def _setup_database_seed_file(galaxy_data):
    _seed_append("""from scripts.db_shell import *
from galaxy.util.bunch import Bunch
from galaxy.security import GalaxyRBACAgent
bunch = Bunch( **globals() )
bunch.engine = engine
# model.flush() has been removed.
bunch.session = db_session
# For backward compatibility with "model.context.current"
bunch.context = db_session
security_agent = GalaxyRBACAgent( bunch )
security_agent.sa_session = sa_session

def add_user(email, password, key=None):
    query = sa_session.query( User ).filter_by( email=email )
    if query.count() > 0:
        return query.first()
    else:
        user = User(email)
        user.set_password_cleartext(password)
        sa_session.add(user)
        sa_session.flush()

        security_agent.create_private_user_role( user )
        if not user.default_permissions:
            security_agent.user_set_default_permissions( user, history=True, dataset=True )

        if key is not None:
            api_key = APIKeys()
            api_key.user_id = user.id
            api_key.key = key
            sa_session.add(api_key)
            sa_session.flush()
        return user

def add_history(user, name):
    query = sa_session.query( History ).filter_by( user=user ).filter_by( name=name )
    if query.count() == 0:
        history = History(user=user, name=name)
        sa_session.add(history)
        sa_session.flush()
        return history
    else:
        return query.first()

""")
    i = 0
    for user in galaxy_data["users"]:
        username = user["username"]
        password = user["password"]
        api_key = user["api_key"]
        histories = None
        if "histories" in user:
            histories = user["histories"]
        user_object = "user_%d" % i
        _seed_append("""%s = add_user("%s", "%s", "%s")""" % (user_object, username, password, api_key))
        _import_histories(user_object, histories)
        i = i + 1


def _import_histories(user_object, histories):
    if not histories:
        return
    for history_name in histories:
        _import_history(user_object, history_name)


def _import_history(user_object, history_name):
    history_name_stripped = history_name.strip()
    if history_name_stripped:
        _seed_append("""add_history(%s, "%s")""" % (user_object, history_name_stripped))


def _seed_append(text):
    append("%s/seed.py" % env.galaxy_home, text, use_sudo=True)


def _start_galaxy():
    # Create directory to store galaxy service's pid file.
    _make_dir_for_galaxy("/var/lib/galaxy")
    start_service("galaxy")


def refresh_galaxy(target_galaxy_repo):
    _update_galaxy(target_galaxy_repo)
    sudo("/etc/init.d/galaxy restart", pty=False)


def _setup_galaxy_log_dir():
    _make_dir_for_galaxy("/var/log/galaxy")


def _setup_shed_tools_dir():
    _make_dir_for_galaxy("%s/../shed_tools" % env.galaxy_home)


def _make_dir_for_galaxy(path):
    sudo("mkdir -p '%s'" % path)
    _chown_galaxy(env, path)


#def _fix_galaxy_permissions():
#    # Ensure that everything under install dir is owned by env.galaxy_user
#    _chown_galaxy(env, os.path.split(env.install_dir)[0])
#    sudo("chmod 755 %s" % os.path.split(env.install_dir)[0])


def _update_galaxy(target_galaxy_repo):
    # Need to merge? -John
    hg_command = "hg pull %s; hg update" % target_galaxy_repo
    with cd(env.galaxy_home):
        sudo(hg_command, user=env.galaxy_user)
