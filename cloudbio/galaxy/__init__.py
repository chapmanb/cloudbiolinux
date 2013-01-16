"""
Adapted from Enis Afgan's mi-deployment code:
https://bitbucket.org/afgane/mi-deployment
"""

import os
import contextlib

from fabric.api import sudo, run, cd, settings, hide
from fabric.contrib.files import exists, contains, sed, append
from fabric.colors import red

from cloudbio.custom.shared import _write_to_file, _setup_conf_file, _setup_simple_service,  _make_tmp_dir
from cloudbio.galaxy.tools import _install_tools
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean, _dir_is_empty


# -- Adjust this link if using content from another location
CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"
CM_REPO_ROOT_URL = "https://bitbucket.org/galaxy/cloudman/raw/tip/"

def _setup_users(env):
    def _add_user(username, uid=None):
        """ Add user with username to the system """
        if not contains('/etc/passwd', "%s:" % username):
            uid_str = "--uid %s" % uid if uid else ""
            sudo('useradd -d /home/%s --create-home --shell /bin/bash ' \
                 '-c"Galaxy-required user" %s --user-group %s' % \
                     (username, uid_str, username))
    # Must specify uid for 'galaxy' user because of the
    # configuration for proFTPd
    _add_user('galaxy', '1001')
    _add_user('sgeadmin')
    _add_user('postgres')
    _add_user('condor')
    env.logger.debug("Done setting up CloudMan users")


def _setup_galaxy_env_defaults(env):
    if "galaxy_user" not in env:
        env.galaxy_user = "galaxy"
    if "galaxy_home" not in env:
        env.galaxy_home = "/mnt/galaxyTools/galaxy-central"
    if "galaxy_tools_dir" not in env:
        # Was called install_dir in tools_fabfile.py
        env.galaxy_tools_dir = "/mnt/galaxyTools/tools"
    if "galaxy_loc_files" not in env:
        indicies_dir = env.get("data_files", "/mnt/galaxyIndices")
        env.galaxy_loc_files = os.path.join(indicies_dir, "galaxy", "galaxy-data")
    if "galaxy_jars_dir" not in env:
        env.galaxy_jars_dir = os.path.join(env.galaxy_home, "tool-data", "shared", "jars")
    if "galaxy_update_default" not in env:
        env.galaxy_update_default = True
    if "python_version" not in env:
        env.python_version = "2.7"  # Override in fabricrc if this is not the case.
    if "galaxy_indices_mount" not in env:
        indicies_dir = env.get("data_files", "/mnt/galaxyIndices")
        env.galaxy_indices_mount = indicies_dir
    if "galaxy_data_mount" not in env:
        env.galaxy_data_mount = "/mnt/galaxyData"
    if "galaxy_init_database" not in env:
        env.galaxy_init_database = False

def _install_galaxy(env):
    """
    Used to install Galaxy and setup its environment, including tools.
    This method is somewhat targeted for the cloud deployment so some
    tweaking of the code may be desirable. This method cannot be used
    to update an existing Galaxy.
    """
    _clone_galaxy_repo(env)
    _chown_galaxy(env, env.galaxy_home) # Make sure env.galaxy_user owns env.galaxy_home
    sudo("chmod 755 %s" % os.path.split(env.galaxy_home)[0])
    setup_db = _read_boolean(env, "galaxy_setup_database", False)
    if setup_db:
        _setup_galaxy_db(env)
    setup_service = _read_boolean(env, "galaxy_setup_service", False)
    if setup_service:
        _setup_service(env)
    _install_tools(env)
    setup_xvfb = _read_boolean(env, "galaxy_setup_xvfb", False)
    if setup_xvfb:
        _setup_xvfb(env)
    _chown_galaxy(env, env.galaxy_home) # Make sure env.galaxy_user owns env.galaxy_home
    _chown_galaxy(env, env.galaxy_loc_files) # Make sure env.galaxy_user owns env.galaxy_loc_files
    return True

def _clone_galaxy_repo(env):
    """
    Clone Galaxy source code repository from ``env.galaxy_repository`` to
    ``env.galaxy_home``, setting the directory ownership to ``env.galaxy_user``

    This method cannot be used to update an existing Galaxy installation.
    """
    # Make sure ``env.galaxy_home`` dir exists but without Galaxy in it
    if exists(env.galaxy_home):
        if exists(os.path.join(env.galaxy_home, '.hg')):
            env.logger.warning("Galaxy install dir {0} exists and seems to have " \
                "a Mercurial repository already there. Galaxy already installed?"\
                .format(env.galaxy_home))
            return False
    else:
        sudo("mkdir -p '%s'" % env.galaxy_home)
    with cd(env.galaxy_home):
        # Needs to be done as non galaxy user, otherwise we have a
        # permissions problem.
        galaxy_repository = env.get("galaxy_repository", 'https://bitbucket.org/galaxy/galaxy-central/')
        env.safe_sudo('hg clone %s .' % galaxy_repository)
    # Make sure ``env.galaxy_home`` is owned by ``env.galaxy_user``
    _chown_galaxy(env, env.galaxy_home)
    # If needed, custom-configure this freshly cloned Galaxy
    preconfigured = _read_boolean(env, "galaxy_preconfigured_repository", False)
    if not preconfigured:
        _configure_galaxy_repository(env)

def _setup_galaxy_db(env):
    """
    Create (if one already does not exist) and initialize a database for use
    by Galaxy.
    """
    if not _galaxy_db_exists(env):
        _create_galaxy_db(env)
    _init_galaxy_db(env)

def _get_galaxy_db_configs(env):
    """
    Extract configuration options for Galaxy database and return those as a dictionary
    """
    config = {}
    config['psql_data_dir'] = env.get('galaxy_database_path', '/mnt/galaxy/db')
    config['psql_conf_file'] = os.path.join(config['psql_data_dir'], 'postgresql.conf')
    config['psql_bin_dir'] = env.get('postgres_bin_dir', '/usr/lib/postgresql/9.1/bin')
    config['psql_user'] = env.get('postrges_user', 'postgres')
    config['psql_port'] = env.get('postgres_port', '5910')
    config['psql_log'] = '/tmp/pSQL.log'
    config['galaxy_db_name'] = env.get('galaxy_db_name', 'galaxy')
    config['galaxy_ftp_user_pwd'] = env.get('galaxy_ftp_user_password', 'fu5yOj2sn')
    # And a couple of useful command shortcuts
    config['pg_ctl_cmd'] = "{0} -D {2}".format(os.path.join(config['psql_bin_dir'], 'pg_ctl'),
        config['psql_port'], config['psql_data_dir'])
    config['pg_start_cmd'] = "{0} -w -l {1} start".format(config['pg_ctl_cmd'], config['psql_log'])
    config['pg_stop_cmd'] = "{0} -w stop".format(config['pg_ctl_cmd'])
    config['psql_cmd'] = "{0} -p {1}".format(os.path.join(config['psql_bin_dir'], 'psql'),
        config['psql_port'])
    return config

def _galaxy_db_exists(env):
    """
    Check if galaxy database already exists. Return ``True`` if it does,
    ``False`` otherwise.

    Note that this method does a best-effort attempt at starting the DB server
    if one is not already running to do a thorough test. It shuts the server
    down upon completion, but only it if also started it.
    """
    db_exists = False
    started = False
    c = _get_galaxy_db_configs(env)
    if exists(c['psql_data_dir']) and not _dir_is_empty(c['psql_data_dir']):
        sudo("chown --recursive {0}:{0} {1}".format(c['psql_user'], c['psql_data_dir']))
        env.logger.debug("Galaxy database directory {0} already exists.".format(c['psql_data_dir']))
        # Check if PostgreSQL is already running and try to start the DB if not
        if not _postgres_running(env):
            with settings(warn_only=True):
                env.logger.debug("Trying to start DB server in {0}".format(c['psql_data_dir']))
                sudo("{0}".format(c['pg_start_cmd']), user=c['psql_user'])
                started = True
        # Check if galaxy DB already exists
        if 'galaxy' in sudo("{0} -P pager --list | grep {1} || true".format(c['psql_cmd'],
            c['galaxy_db_name']), user=c['psql_user']):
            env.logger.warning("Galaxy database {0} already exists in {1}! Not creating it."\
                .format(c['galaxy_db_name'], c['psql_data_dir']))
            db_exists = True
        if started:
            with settings(warn_only=True):
                sudo("{0}".format(c['pg_stop_cmd']), user=c['psql_user'])
    return db_exists

def _create_galaxy_db(env):
    """
    Create a new PostgreSQL database for use by Galaxy
    """
    c = _get_galaxy_db_configs(env)
    if not exists(c['psql_data_dir']):
        sudo("mkdir -p {0}".format(c['psql_data_dir']))
    sudo("chown --recursive {0}:{0} {1}".format(c['psql_user'], c['psql_data_dir']))
    # Initialize a new database for Galaxy in ``psql_data_dir``
    if _dir_is_empty(c['psql_data_dir']):
        sudo("{0} -D {1}".format(os.path.join(c['psql_bin_dir'], 'initdb'), c['psql_data_dir']),
            user=c['psql_user'])
    # Set port for the database server
    sed(c['psql_conf_file'], '#port = 5432', 'port = {0}'.format(c['psql_port']), use_sudo=True)
    sudo("chown {0}:{0} {1}".format(c['psql_user'], c['psql_conf_file']))
    # Start PostgreSQL server so a role for Galaxy user can be created
    if not _postgres_running(env):
        sudo(c['pg_start_cmd'], user=c['psql_user'])
        started = True
    else:
        # Restart is required so port setting takes effect
        sudo("{0} -D {1} -w -l {2} restart".format(c['pg_ctl_cmd']), c['psql_data_dir'],
            c['psql_log'], user=c['psql_user'])
        started = False
    # Create a role for env.galaxy_user
    sudo('{0} -c"CREATE ROLE {1} LOGIN CREATEDB"'.format(c['psql_cmd'], env.galaxy_user),
        user=c['psql_user'])
    # Create a Galaxy database
    sudo('{0} -p {1} {2}'.format(os.path.join(c['psql_bin_dir'], 'createdb'),
        c['psql_port'], c['galaxy_db_name']), user=env.galaxy_user)
    # Create a role for 'galaxyftp' user
    sudo('{0} -c"CREATE ROLE galaxyftp LOGIN PASSWORD \'{1}\'"'.format(c['psql_cmd'],
        c['galaxy_ftp_user_pwd']), user=c['psql_user'])
    if started:
        with settings(warn_only=True):
            sudo("{0}".format(c['pg_stop_cmd']), user=c['psql_user'])
    exp = "export PATH={0}:$PATH".format(c['psql_bin_dir'])
    if not contains('/etc/bash.bashrc', exp):
        append('/etc/bash.bashrc', exp, use_sudo=True)

def _init_galaxy_db(env):
    """
    Initialize Galaxy's database with the tables and apply the migrations,
    fetching Galaxy eggs in the process.
    """
    with cd(env.galaxy_home):
        universe_wsgi_url = env.get('galaxy_universe_wsgi_url',
            os.path.join(CM_REPO_ROOT_URL, 'universe_wsgi.ini.cloud'))
        sudo("wget --output-document=universe_wsgi.ini {0}".format(universe_wsgi_url))
        started = False
        if not _postgres_running(env):
            c = _get_galaxy_db_configs(env)
            sudo(c['pg_start_cmd'], user=c['psql_user'])
            started = True
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; python -ES ./scripts/fetch_eggs.py; ./create_db.sh'",
            user=env.galaxy_user)
        sudo("rm universe_wsgi.ini")
        if started:
            with settings(warn_only=True):
                sudo("{0}".format(c['pg_stop_cmd']), user=c['psql_user'])

def _configure_galaxy_options(env, option_dict=None, prefix="galaxy_universe_"):
    """
    Read through fabric's environment and make sure any property of
    the form galaxy_universe_XXX=YYY lands up in Galaxy's universe_wsgi.ini
    options as XXX=YYY using Galaxy configuration directory:
    """
    galaxy_conf_directory = env.get("galaxy_conf_directory", None)
    if not galaxy_conf_directory:
        return False
    # By default just read the options from env (i.e. from fabricrc), but
    # allow override so the options can come from a YAML file (such as done
    # with galaxy-vm-launcher.)
    if option_dict is None:
        option_dict = env

    option_priority = env.get("galaxy_conf_priority", "200")
    for key, value in option_dict.iteritems():
        if key.startswith(prefix):
            key = key[len(prefix):]
            conf_file_name = "%s_override_%s.ini" % (option_priority, key)
            conf_file = os.path.join(galaxy_conf_directory, conf_file_name)
            contents = "[app:main]\n%s=%s" % (key, value)
            _write_to_file(contents, conf_file, 0700)
            _chown_galaxy(env, conf_file)

def _configure_galaxy_repository(env):
    """
    Custom-configure Galaxy repository. This is primarily targeted at a cloud
    deployment.

    mi-deployment would always edit the repository in this way, but
    galaxy-vm-launcher requires the configured Galaxy repository to pull
    in the changesets from https://bitbucket.org/jmchilton/cloud-galaxy-dist
    which prebakes these modifications in.
    """
    _chown_galaxy(env, env.galaxy_home) # Make sure env.galaxy_user owns env.galaxy_home
    with cd(env.galaxy_home):  # and settings(warn_only=True):
        # Make sure Galaxy runs in a new shell and does not
        # inherit the environment by adding the '-ES' flag
        # to all invocations of python within run.sh
        sudo("sed -i 's/python .\//python -ES .\//g' run.sh", user=env.galaxy_user)
        if _read_boolean(env, "galaxy_cloud", True):
            # Append DRMAA_LIBRARY_PATH in run.sh as well (this file will exist
            # once SGE is installed - which happens at instance contextualization)
            sudo("grep -q 'export DRMAA_LIBRARY_PATH=/opt/sge/lib/lx24-amd64/libdrmaa.so.1.0' run.sh; if [ $? -eq 1 ]; then sed -i '2 a export DRMAA_LIBRARY_PATH=/opt/sge/lib/lx24-amd64/libdrmaa.so.1.0' run.sh; fi", user=env.galaxy_user)
            # Upload the custom cloud welcome screen files
            if not exists("%s/static/images/cloud.gif" % env.galaxy_home):
                sudo("wget --output-document=%s/static/images/cloud.gif %s/cloud.gif" % (env.galaxy_home, CDN_ROOT_URL), user=env.galaxy_user)
            if not exists("%s/static/images/cloud_txt.png" % env.galaxy_home):
                sudo("wget --output-document=%s/static/images/cloud_text.png %s/cloud_text.png" % (env.galaxy_home, CDN_ROOT_URL), user=env.galaxy_user)
            sudo("wget --output-document=%s/static/welcome.html %s/welcome.html" % (env.galaxy_home, CDN_ROOT_URL), user=env.galaxy_user)
        # Set up the symlink for SAMTOOLS (remove this code once SAMTOOLS is converted to data tables)
        if exists("%s/tool-data/sam_fa_indices.loc" % env.galaxy_home):
            sudo("rm %s/tool-data/sam_fa_indices.loc" % env.galaxy_home, user=env.galaxy_user)
        tmp_loc = False
        if not exists(env.galaxy_loc_files):
            sudo("mkdir -p '%s'" % env.galaxy_loc_files)
            _chown_galaxy(env, env.galaxy_loc_files)
        if not exists("%s/sam_fa_indices.loc" % env.galaxy_loc_files):
            sudo("touch %s/sam_fa_indices.loc" % env.galaxy_loc_files, user=env.galaxy_user)
            tmp_loc = True
        sudo("ln -s %s/sam_fa_indices.loc %s/tool-data/sam_fa_indices.loc" % (env.galaxy_loc_files, env.galaxy_home), user=env.galaxy_user)
        if tmp_loc:
            sudo("rm %s/sam_fa_indices.loc" % env.galaxy_loc_files)
        # set up the special HYPHY link in tool-data/
        hyphy_dir = os.path.join(env.galaxy_tools_dir, 'hyphy', 'default')
        sudo('ln -s %s tool-data/HYPHY' % hyphy_dir, user=env.galaxy_user)
        # set up the jars directory for Java tools
        if not exists('tool-data/shared/jars'):
            sudo("mkdir -p tool-data/shared/jars", user=env.galaxy_user)
        srma_dir = os.path.join(env.galaxy_tools_dir, 'srma', 'default')
        haploview_dir = os.path.join(env.galaxy_tools_dir, 'haploview', 'default')
        picard_dir = os.path.join(env.galaxy_tools_dir, 'picard', 'default')
        sudo('ln -s %s/srma.jar tool-data/shared/jars/.' % srma_dir, user=env.galaxy_user)
        sudo('ln -s %s/haploview.jar tool-data/shared/jars/.' % haploview_dir, user=env.galaxy_user)
        sudo('ln -s %s/*.jar tool-data/shared/jars/.' % picard_dir, user=env.galaxy_user)
    return True


def _setup_service(env):
    _setup_conf_file(env, "/etc/init.d/galaxy", "galaxy_init", default_source="galaxy_init")
    _setup_conf_file(env, "/etc/default/galaxy", "galaxy_default")
    _setup_simple_service("galaxy")


def _setup_xvfb(env):
    _setup_conf_file(env, "/etc/init.d/xvfb", "xvfb_init", default_source="xvfb_init")
    _setup_conf_file(env, "/etc/default/xvfb", "xvfb_default", default_source="xvfb_default")
    _setup_simple_service("xvfb")
    env.safe_sudo("mkdir /var/lib/xvfb; chown root:root /var/lib/xvfb; chmod 0755 /var/lib/xvfb")


def _setup_nginx_service(env):
    # Setup system service for nginx, not needed for CloudMan but it is
    # useful if CloudMan is not being used (such as galaxy-vm-launcher work).
    _setup_conf_file(env, "/etc/init.d/nginx", "nginx_init", default_source="nginx_init")
    _setup_simple_service("nginx")


def _install_nginx_standalone(env):
    _install_nginx(env)
    _setup_nginx_service(env)


def _install_nginx(env):
    """Nginx open source web server.
    http://www.nginx.org/
    """
    version = "1.2.0"
    url = "http://nginx.org/download/nginx-%s.tar.gz" % version

    install_dir = os.path.join(env.install_dir, "nginx")
    remote_conf_dir = os.path.join(install_dir, "conf")

    # Skip install if already present
    if exists(remote_conf_dir) and contains(os.path.join(remote_conf_dir, "nginx.conf"), "/cloud"):
        env.logger.debug("Nginx already installed; not installing it again.")
        return

    with _make_tmp_dir() as work_dir:
        with contextlib.nested(cd(work_dir), settings(hide('stdout'))):
            modules = _get_nginx_modules(env)
            module_flags = " ".join(["--add-module=../%s" % x for x in modules])
            run("wget %s" % url)
            run("tar xvzf %s" % os.path.split(url)[1])
            with cd("nginx-%s" % version):
                run("./configure --prefix=%s --with-ipv6 %s "
                    "--user=galaxy --group=galaxy --with-debug "
                    "--with-http_ssl_module --with-http_gzip_static_module" %
                    (install_dir, module_flags))
                sed("objs/Makefile", "-Werror", "")
                run("make")
                sudo("make install")
                sudo("cd %s; stow nginx" % env.install_dir)

    defaults = {"galaxy_home": "/mnt/galaxyTools/galaxy-central"}
    _setup_conf_file(env, os.path.join(remote_conf_dir, "nginx.conf"), "nginx.conf", defaults=defaults)

    nginx_errdoc_file = 'nginx_errdoc.tar.gz'
    url = os.path.join(REPO_ROOT_URL, nginx_errdoc_file)
    remote_errdoc_dir = os.path.join(install_dir, "html")
    with cd(remote_errdoc_dir):
        sudo("wget --output-document=%s/%s %s" % (remote_errdoc_dir, nginx_errdoc_file, url))
        sudo('tar xvzf %s' % nginx_errdoc_file)

    sudo("mkdir -p %s" % env.install_dir)
    if not exists("%s/nginx" % env.install_dir):
        sudo("ln -s %s/sbin/nginx %s/nginx" % (install_dir, env.install_dir))
    # If the guessed symlinking did not work, force it now
    cloudman_default_dir = "/opt/galaxy/sbin"
    if not exists(cloudman_default_dir):
        sudo("mkdir -p %s" % cloudman_default_dir)
    if not exists(os.path.join(cloudman_default_dir, "nginx")):
        sudo("ln -s %s/sbin/nginx %s/nginx" % (install_dir, cloudman_default_dir))
    env.logger.debug("Nginx {0} installed to {1}".format(version, install_dir))


def _get_nginx_modules(env):
    """Retrieve add-on modules compiled along with nginx.
    """
    modules = {
        "upload": True,
        "chunk": True,
        "ldap": False
    }

    module_dirs = []

    for module, enabled_by_default in modules.iteritems():
        enabled = _read_boolean(env, "nginx_enable_module_%s" % module, enabled_by_default)
        if enabled:
            module_dirs.append(eval("_get_nginx_module_%s" % module)(env))

    return module_dirs


def _get_nginx_module_upload(env):
    upload_module_version = "2.2.0"
    upload_url = "http://www.grid.net.ru/nginx/download/" \
                 "nginx_upload_module-%s.tar.gz" % upload_module_version
    run("wget %s" % upload_url)
    upload_fname = os.path.split(upload_url)[1]
    run("tar -xvzpf %s" % upload_fname)
    return upload_fname.rsplit(".", 2)[0]


def _get_nginx_module_chunk(env):
    chunk_module_version = "0.22"
    chunk_git_version = "b46dd27"

    chunk_url = "https://github.com/agentzh/chunkin-nginx-module/tarball/v%s" % chunk_module_version
    chunk_fname = "agentzh-chunkin-nginx-module-%s.tar.gz" % (chunk_git_version)
    run("wget -O %s %s" % (chunk_fname, chunk_url))
    run("tar -xvzpf %s" % chunk_fname)
    return chunk_fname.rsplit(".", 2)[0]


def _get_nginx_module_ldap(env):
    run("rm -rf nginx-auth-ldap")  # Delete it if its there or git won't clone
    run("git clone https://github.com/kvspb/nginx-auth-ldap")
    return "nginx-auth-ldap"


def _setup_postgresql(env):
    # Handled by CloudMan, but if configuring standalone galaxy, this
    # will need to be executed to create a postgres user for Galaxy.
    _configure_postgresql(env)
    _init_postgresql_data(env)


def _configure_postgresql(env, delete_main_dbcluster=False):
    """ This method is intended for cleaning up the installation when
    PostgreSQL is installed from a package. Basically, when PostgreSQL
    is installed from a package, it creates a default database cluster
    and splits the config file away from the data.
    This method can delete the default database cluster that was automatically
    created when the package is installed. Deleting the main database cluster
    also has the effect of stopping the auto-start of the postmaster server at
    machine boot. The method adds all of the PostgreSQL commands to the PATH.
    """
    pg_ver = sudo("dpkg -s postgresql | grep Version | cut -f2 -d':'")
    pg_ver = pg_ver.strip()[:3]  # Get first 3 chars of the version since that's all that's used for dir name
    got_ver = False
    while(not got_ver):
        try:
            pg_ver = float(pg_ver)
            got_ver = True
        except Exception:
            print(red("Problems trying to figure out PostgreSQL version."))
            pg_ver = raw_input(red("Enter the correct one (eg, 9.1; not 9.1.3): "))
    if delete_main_dbcluster:
        env.safe_sudo('pg_dropcluster --stop %s main' % pg_ver, user='postgres')
    # Not sure why I ever added this to gvl, doesn't seem needed. -John
    #_put_installed_file_as_user("postgresql-%s.conf" % env.postgres_version, "/etc/postgresql/%s/main/postgresql.conf" % env.postgres_version, user='root')
    exp = "export PATH=/usr/lib/postgresql/%s/bin:$PATH" % pg_ver
    if not contains('/etc/bash.bashrc', exp):
        append('/etc/bash.bashrc', exp, use_sudo=True)


def _init_postgresql_data(env):
    if "galaxy" not in env.safe_sudo("psql -P pager --list | grep galaxy || true", user="postgres"):
        env.safe_sudo("createdb galaxy", user="postgres")
        env.safe_sudo("psql -c 'create user galaxy; grant all privileges on database galaxy to galaxy;'", user="postgres")

def _postgres_running(env):
    """
    Return ``True`` if PostgreSQL is running, ``False`` otherwise.
    """
    c = _get_galaxy_db_configs(env)
    if 'no server running' in sudo("{0} status || true".format(c['pg_ctl_cmd']), user=c['psql_user']):
        return False
    return True
