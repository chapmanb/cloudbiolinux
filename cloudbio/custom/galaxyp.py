"""
"""

from cloudbio.galaxy.utils import _chown_galaxy

from fabric.contrib.files import *

from shared import _write_to_file


def install_protvis(env):
    """ Installs Andrew Brock's proteomics visualize tool.
    https://bitbucket.org/Andrew_Brock/proteomics-visualise/
    """
    _setup_protvis_env(env)
    protvis_home = env["protvis_home"]
    env.safe_sudo("sudo apt-get -y --force-yes install libxml2-dev libxslt-dev")

    run("rm -rf protvis")
    run("git clone -b lorikeet https://github.com/jmchilton/protvis.git")
    with cd("protvis"):
        run("git submodule init")
        run("git submodule update")
        env.safe_sudo("rsync -avur --delete-after . %s" % (protvis_home))
        _chown_galaxy(env, protvis_home)
        with cd(protvis_home):
            env.safe_sudo("./setup.sh", user=env.get("galaxy_user", "galaxy"))

    #default_revision = "8cc6af1c492c"
    #
    #revision = env.get("protvis_revision", default_revision)
    #url = _get_bitbucket_download_url(revision, "https://bitbucket.org/Andrew_Brock/proteomics-visualise")
    #def _make(env):
    #_get_install(url, env, _make)

    galaxy_data_dir = env.get('galaxy_data_dir', "/mnt/galaxyData/")
    protvis_converted_files_dir = env.get('protvis_converted_files_dir')
    _write_to_file('''GALAXY_ROOT = "%s"
PATH_WHITELIST = ["%s/files/", "%s"]
CONVERTED_FILES = "%s"
''' % (env.galaxy_home, galaxy_data_dir, protvis_converted_files_dir, protvis_converted_files_dir), "%s/conf.py" % protvis_home, 0755)
    _setup_protvis_service(env)


def _setup_protvis_env(env):
    if not "protvis_home" in env:
        env["protvis_home"] = "%s/%s" % (env.galaxy_tools_dir, "protvis")
    if not "protvis_user" in env:
        env["protvis_user"] = "galaxy"
    if not "protvis_port" in env:
        env["protvis_port"] = "8500"
    if not "protvis_converted_files_dir" in env:
        galaxy_data_dir = env.get('galaxy_data_dir', "/mnt/galaxyData/")
        env['protvis_converted_files_dir'] = "%s/tmp/protvis" % galaxy_data_dir


def _setup_protvis_service(env):
    _setup_conf_file(env, os.path.join("/etc/init.d/protvis"), "protvis_init", default_source="protvis_init")
    _setup_conf_file(env, os.path.join("/etc/default/protvis"), "protvis_default")
    _setup_simple_service("protvis")
