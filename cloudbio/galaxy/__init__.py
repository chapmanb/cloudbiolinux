"""
Adapted from Enis Afgan's mi-deployment code:
https://bitbucket.org/afgane/mi-deployment
"""

import os

from fabric.api import sudo, run, cd
from fabric.contrib.files import exists, contains

from cloudbio.custom.shared import _write_to_file, _setup_conf_file, _setup_simple_service
from cloudbio.galaxy.tools import _install_tools
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean

# -- Adjust this link if using content from another location
CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"


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
