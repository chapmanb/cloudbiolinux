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
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean


# -- Adjust this link if using content from another location
CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"


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
        indicies_dir = env.get("data_files", "/mnt/galaxyIndcies")
        env.galaxy_loc_files = os.path.join(indicies_dir, "galaxy", "galaxy-data")
    if "galaxy_jars_dir" not in env:
        env.galaxy_jars_dir = os.path.join(env.galaxy_home, "tool-data", "shared", "jars")
    if "galaxy_update_default" not in env:
        env.galaxy_update_default = True
    if "python_version" not in env:
        env.python_version = "2.7"  # Override in fabricrc if this is not the case.


def _install_galaxy(env):
    """ Used to install Galaxy and setup its environment.
    This method cannot be used to update an existing instance of Galaxy code;
    see volume_manipulations_fab.py script for that functionality.
    Also, this method is somewhat targeted for the EC2 deployment so some
    tweaking of the code may be desirable."""
    _clone_galaxy_repo(env)
    # MP: Ensure that everything under install dir is owned by env.galaxy_user
    sudo("chown --recursive %s:%s %s"
       % (env.galaxy_user, env.galaxy_user, os.path.split(env.galaxy_tools_dir)[0]))
    sudo("chmod 755 %s" % os.path.split(env.galaxy_tools_dir)[0])
    setup_service = _read_boolean(env, "galaxy_setup_service", False)
    if setup_service:
        _setup_service(env)
    _install_tools(env)
    setup_xvfb = _read_boolean(env, "galaxy_setup_xvfb", False)
    if setup_xvfb:
        _setup_xvfb(env)
    return True

def _clone_galaxy_repo(env):
    # MP: we need to have a tmp directory available if files already exist
    # in the galaxy install directory
    install_cmd = sudo if env.get("use_sudo", True) else run

    base_tmp_dir = env.get("galaxy_tmp_dir", "/mnt")
    tmp_dir = os.path.join(base_tmp_dir, "fab_tmp")
    if exists(tmp_dir):
        install_cmd("rm -rf %s" % tmp_dir)
    if exists(env.galaxy_home):
        if exists(os.path.join(env.galaxy_home, '.hg')):
            env.logger.info("Galaxy install dir '%s' exists and seems to have a Mercurial repository already there. Galaxy already installed.")
            return False
        else:
            # MP: need to move any files already in galaxy home so that hg
            # can checkout files.
            if not exists(tmp_dir):
                install_cmd("mkdir %s" % tmp_dir)
                install_cmd("chown %s %s" % (env.user, tmp_dir))
            install_cmd("mv %s/* %s" % (env.galaxy_home, tmp_dir))
    ## This is slightly different than mi-deployment to handle the
    ## case when the bucket url doesn't match the desired path.
    if not exists(env.galaxy_home):
        sudo("mkdir -p '%s'" % env.galaxy_home)
        _chown_galaxy(env, env.galaxy_home)
    with cd(env.galaxy_home):
        # MP needs to be done as non galaxy user, otherwise we have a
        # permissions problem.
        galaxy_repository = env.get("galaxy_repository", 'https://bitbucket.org/galaxy/galaxy-central/')
        sudo('hg clone %s .' % galaxy_repository)
    # MP: now we need to move the files back into the galaxy directory.
    if exists(tmp_dir):
        install_cmd("cp -R %s/* %s" % (tmp_dir, env.galaxy_home))
        install_cmd("rm -rf %s" % tmp_dir)
    preconfigured = _read_boolean(env, "galaxy_preconfigured_repository", False)
    if not preconfigured:
        _configure_galaxy_repository(env)


def _configure_galaxy_options(env, option_dict=None, prefix="galaxy_universe_"):
    """
    Read through fabric's environment and make sure any property of
    the form galaxy_universe_XXX=YYY, lands up in Galaxy's universe_wsgi.ini
    options as XXX=YYY using John Chilton's configuration directory work:

    https://bitbucket.org/galaxy/galaxy-central/pull-request/44/allow-usage-of-directory-of-configuration

    Until, the above pull request is accepted, its changeset should be pulled
    into the configured Galaxy repository.
    """
    galaxy_conf_directory = env.get("galaxy_conf_directory", False)
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
    """ mi-deployment would always edit the repository in this way, but
    galaxy-vm-launcher requires the configured Galaxy repository to pull
    in the changesets from https://bitbucket.org/jmchilton/cloud-galaxy-dist
    which prebakes these modifications in.
    """
    sudo("chown -R %s %s" % (env.galaxy_user, env.galaxy_home))
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
            sudo("rm %s/sam_fa_indices.loc" % env.galaxy_loc_files, user=env.galaxy_user)
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
    run("git clone https://code.google.com/p/nginx-auth-ldap/")
    return "nginx-auth-ldap"


def _setup_postgresql(env):
    # Handled by CloudMan, but if configuring standalone galaxy, this
    # will need to be executed to create a postgres user for Galaxy.
    _configure_postgresql(env)
    _init_postgresql_data()


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
        sudo('pg_dropcluster --stop %s main' % pg_ver, user='postgres')
    # Not sure why I ever added this to gvl, doesn't seem needed. -John
    #_put_installed_file_as_user("postgresql-%s.conf" % env.postgres_version, "/etc/postgresql/%s/main/postgresql.conf" % env.postgres_version, user='root')
    exp = "export PATH=/usr/lib/postgresql/%s/bin:$PATH" % pg_ver
    if not contains('/etc/bash.bashrc', exp):
        append('/etc/bash.bashrc', exp, use_sudo=True)


def _init_postgresql_data():
    if "galaxy" not in sudo("psql -P pager --list | grep galaxy || true", user="postgres"):
        sudo("createdb galaxy", user="postgres")
        sudo("psql -c 'create user galaxy; grant all privileges on database galaxy to galaxy;'", user="postgres")
