"""
"""

from cloudbio.custom.galaxy import _prep_galaxy
from cloudbio.custom.shared import _setup_conf_file, _setup_simple_service
from cloudbio.galaxy.utils import _chown_galaxy

from fabric.context_managers import prefix
from fabric.contrib.files import *

from shared import _write_to_file

import yaml


def install_galaxy_protk(env):
    """This method installs Ira Cooke's ProtK framework for the galaxy user.

    By default this will install ProtK from rubygems server, but if
    env.protk_version is set to <version>@<url> (e.g.
    1.1.5@https://bitbucket.org/iracooke/protk-working) the
    gem will be cloned with hg and installed from source.
    """
    _prep_galaxy(env)
    default_version = "1.2.2"
    version_and_revision = env.get("protk_version", default_version)
    install_from_source = version_and_revision.find("@") > 0
    # e.g. protk_version = 1.1.5@https://bitbucket.org/iracooke/protk-working
    if install_from_source:
        (version, revision) = version_and_revision.split("@")
        url = _get_bitbucket_download_url(revision, "https://bitbucket.org/iracooke/protk")
    else:
        version = version_and_revision

    ruby_version = "1.9.3"
    force_rvm_install = False
    with prefix("HOME=~%s" % env.galaxy_user):
        def rvm_exec(env, cmd="", rvm_cmd="use", with_gemset=False):
            target = ruby_version if not with_gemset else "%s@%s" % (ruby_version, "protk-%s" % version)
            prefix = ". $HOME/.rvm/scripts/rvm; rvm %s %s; " % (rvm_cmd, target)
            env.safe_sudo("%s %s" % (prefix, cmd), user=env.galaxy_user)
        if not exists("$HOME/.rvm") or force_rvm_install:
            env.safe_sudo("curl -L get.rvm.io | bash -s stable; source ~%s/.rvm/scripts/rvm" % (env.galaxy_user), user=env.galaxy_user)
            rvm_exec(env, rvm_cmd="install")
            rvm_exec(env, cmd="rvm gemset create protk-%s" % version)
        if not install_from_source:
            # Typical rubygem install
            rvm_exec(env, "gem install  --no-ri --no-rdoc protk -v %s" % version, with_gemset=True)
        else:
            with cd("~%s" % env.galaxy_user):
                env.safe_sudo("rm -rf protk_source; hg clone '%s' protk_source" % url, user=env.galaxy_user)
                rvm_exec(env, "cd protk_source; gem build protk.gemspec; gem install protk", with_gemset=True)

        protk_properties = {}
        ## ProtK can set these up itself, should make that an option.
        protk_properties["tpp_root"] = os.path.join(env.galaxy_tools_dir, "transproteomic_pipeline", "default")
        protk_properties['openms_root'] = "/usr"  # os.path.join(env.galaxy_tools_dir, "openms", "default", "bin")
        ### Assumes omssa, blast, and transproteomic_pipeline CBL galaxy installs.
        protk_properties['omssa_root'] = os.path.join(env.galaxy_tools_dir, "omssa", "default", "bin")
        protk_properties['blast_root'] = os.path.join(env.galaxy_tools_dir, "blast+", "default")
        protk_properties['pwiz_root'] = os.path.join(env.galaxy_tools_dir, "transproteomic_pipeline", "default", "bin")
        # Other properties: log_file, blast_root
        env.safe_sudo("mkdir -p \"$HOME/.protk\"", user=env.galaxy_user)
        env.safe_sudo("mkdir -p \"$HOME/.protk/Databases\"", user=env.galaxy_user)

        _write_to_file(yaml.dump(protk_properties), "/home/%s/.protk/config.yml" % env.galaxy_user, 0755)

        rvm_exec(env, "protk_setup.rb galaxyenv", with_gemset=True)

        install_dir = os.path.join(env.galaxy_tools_dir, "galaxy_protk", version)
        env.safe_sudo("mkdir -p '%s'" % install_dir)
        _chown_galaxy(env, install_dir)
        env.safe_sudo('ln -s -f "$HOME/.protk/galaxy/env.sh" "%s/env.sh"' % install_dir, user=env.galaxy_user)
        with cd(install_dir):
            with cd(".."):
                env.safe_sudo("ln -s -f '%s' default" % version)


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
