"""
"""

from cloudbio.custom.galaxy import _prep_galaxy
from cloudbio.custom.shared import _setup_conf_file, _setup_simple_service
from cloudbio.galaxy.utils import _chown_galaxy
from cloudbio.galaxy.tools import _install_tools

from fabric.context_managers import prefix
from fabric.contrib.files import *

from shared import _get_install, _write_to_file

import yaml


def install_proteomics_tools(env):
    _prep_galaxy(env)
    galaxyp_tools_conf = yaml.load(open("contrib/flavor/cloudman/galaxyp_tools.yaml", "r"))
    _install_tools(env, galaxyp_tools_conf)


def install_protk(env):
    """Installs Ira Cooke's ProtK framework for the
    galaxy user"""
    _prep_galaxy(env)
    default_version = "0.95@e81050c1c658"
    version_and_revision = env.get("protk_version", default_version)
    (version, revision) = version_and_revision.split("@")
    url = _get_bitbucket_download_url(revision, "https://bitbucket.org/iracooke/protk")
    with prefix("HOME=~%s" % env.galaxy_user):
        if not exists("$HOME/.rvm"):
            env.safe_sudo("curl -L get.rvm.io | bash -s stable; source ~%s/.rvm/scripts/rvm" % (env.galaxy_user), user=env.galaxy_user)
            env.safe_sudo(". $HOME/.rvm/scripts/rvm; rvm install 1.8.7 | cat", user=env.galaxy_user)
            env.safe_sudo(". $HOME/.rvm/scripts/rvm; rvm 1.8.7 do gem install rake --no-rdoc --no-ri", user=env.galaxy_user)

    install_dir = os.path.join(env.galaxy_tools_dir, "protk", version)

    def _make(env):
        env.safe_sudo("PROTK_DIR=%s; rm -rf $PROTK_DIR; mkdir -p $PROTK_DIR; mv * $PROTK_DIR" % (install_dir))

    _get_install(url, env, _make)
    _chown_galaxy(env, install_dir)
    with prefix("HOME=~%s" % env.galaxy_user):
        import StringIO
        output = StringIO.StringIO()
        get("%s/config.yml.sample" % install_dir, output)
        sample_config = yaml.load(output.getvalue())
        sample_config['test']['tpp_bin'] = os.path.join(env.galaxy_tools_dir, "transproteomic_pipeline", "default", "bin")
        sample_config['test']['omssa_bin'] = os.path.join(env.galaxy_tools_dir, "omssa", "default", "bin")
        sample_config['test']['ncbi_tools_bin'] = os.path.join(env.galaxy_tools_dir, "blast", "default", "bin")
        sample_config['test']['openms_bin'] = os.path.join(env.galaxy_tools_dir, "openms", "default", "bin")
        #sample_config['test']['galaxy_root']=env.galaxy_home
        _write_to_file(yaml.dump(sample_config), "%s/config.yml" % install_dir, 0755)
        with cd(install_dir):
            #env.safe_sudo("mv /tmp/config.yaml .")
            env.safe_sudo(". $HOME/.rvm/scripts/rvm; echo No | ./setup.sh", user=env.galaxy_user)
            env.safe_sudo("gcc -o make_random make_random.c -lm; ln -s ../make_random bin", user=env.galaxy_user)
    with cd("%s/.." % install_dir):
        env.safe_sudo("ln -s -f %s default" % version)


def install_protvis(env):
    """ Installs Andrew Brock's proteomics visualize tool.
    https://bitbucket.org/Andrew_Brock/proteomics-visualise/
    """
    _setup_protvis_env(env)
    default_revision = "8cc6af1c492c"
    protvis_home = env["protvis_home"]
    revision = env.get("protvis_revision", default_revision)
    url = _get_bitbucket_download_url(revision, "https://bitbucket.org/Andrew_Brock/proteomics-visualise")
    galaxy_data_dir = env.get('galaxy_data_dir', "/mnt/galaxyData/")
    protvis_converted_files_dir = env.get('protvis_converted_files_dir')

    env.safe_sudo("sudo apt-get -y --force-yes install libxml2-dev libxslt-dev")

    def _make(env):
        env.safe_sudo("rsync -avur --delete-after . %s" % (protvis_home))
        _chown_galaxy(env, protvis_home)
        with cd(protvis_home):
            env.safe_sudo("./setup.sh", user=env.get("galaxy_user", "galaxy"))
    _get_install(url, env, _make)
    _write_to_file('''GALAXY_ROOT = "%s"
PATH_WHITELIST = ["%s/files/", "%s"]
CONVERTED_FILES = "%s"
''' % (env.galaxy_home, galaxy_data_dir, protvis_converted_files_dir, protvis_converted_files_dir), "%s/conf.py" % protvis_home, 0755)
    _setup_protvis_service(env)


def _get_bitbucket_download_url(revision, default_repo):
    if revision.startswith("http"):
        url = revision
    else:
        url = "%s/get/%s.tar.gz" % (default_repo, revision)
    return url


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
