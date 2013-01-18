"""Build instructions associated with CloudMan.

http://wiki.g2.bx.psu.edu/Admin/Cloud

Adapted from Enis Afgan's code: https://bitbucket.org/afgane/mi-deployment
"""

cm_upstart = """
description     "Start CloudMan contextualization script"

start on runlevel [2345]

task
exec python %s 2> %s.log
"""
import os

from fabric.api import sudo, run, put, cd
from fabric.contrib.files import exists, settings, append

from cloudbio.galaxy import _setup_users
from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.custom.shared import (_make_tmp_dir, _write_to_file, _get_install,
                                    _configure_make,_if_not_installed)
from cloudbio.package.deb import (_apt_packages, _setup_apt_automation)

MI_REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"
CM_REPO_ROOT_URL = "https://bitbucket.org/galaxy/cloudman/raw/tip"

def _configure_cloudman(env, use_repo_autorun=False):
    _setup_users(env)
    _setup_env(env)
    _configure_logrotate(env)
    _configure_ec2_autorun(env, use_repo_autorun)
    _configure_sge(env)
    _configure_hadoop(env)
    _configure_condor(env)
    _configure_nfs(env)
    install_s3fs(env)

def _setup_env(env):
    """
    Setup the system environment required to run CloudMan. This means
    installing required system-level packages (as defined in CBL's
    ``packages.yaml``, or a flavor thereof) and Python dependencies
    (i.e., libraries) as defined in CloudMan's ``requirements.txt`` file.
    """
    # Get and install required system packages
    if env.distribution in ["debian", "ubuntu"]:
        config_file = get_config_file(env, "packages.yaml")
        (packages, _) = _yaml_to_packages(config_file.base, 'cloudman')
        # Allow editions and flavors to modify the package list
        packages = env.edition.rewrite_config_items("packages", packages)
        packages = env.flavor.rewrite_config_items("packages", packages)
        _setup_apt_automation()
        _apt_packages(pkg_list=packages)
    elif env.distribution in ["centos", "scientificlinux"]:
        env.logger.warn("No CloudMan system package dependencies for CentOS")
        pass
    # Get and install required Python libraries
    reqs_file = 'requirements.txt'
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            url = os.path.join(CM_REPO_ROOT_URL, reqs_file)
            run("wget --output-document=%s %s" % (reqs_file, url))
            sudo("pip install --upgrade --requirement={0}".format(reqs_file))
    # Add a custom vimrc
    vimrc_url = os.path.join(MI_REPO_ROOT_URL, 'conf_files', 'vimrc')
    remote_file = '/etc/vim/vimrc'
    sudo("wget --output-document=%s %s" % (remote_file, vimrc_url))
    env.logger.debug("Added a custom vimrc to {0}".format(remote_file))
    env.logger.debug("Done setting up CloudMan's environment")

def _configure_logrotate(env):
    """
    Add logrotate config file, which will automatically rotate CloudMan's log
    """
    conf_file = "cloudman.logrotate"
    remote = '/etc/logrotate.d/cloudman'
    url = os.path.join(MI_REPO_ROOT_URL, 'conf_files', conf_file)
    sudo("wget --output-document=%s %s" % (remote, url))
    env.logger.debug("----- Added logrotate file to {0} -----".format(remote))

def _configure_ec2_autorun(env, use_repo_autorun=False):
    script = "ec2autorun.py"
    remote = os.path.join(env.install_dir, "bin", script)
    if not exists(os.path.dirname(remote)):
        sudo('mkdir -p {0}'.format(os.path.dirname(remote)))
    if use_repo_autorun:
        url = os.path.join(MI_REPO_ROOT_URL, script)
        sudo("wget --output-document=%s %s" % (remote, url))
    else:
        install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
        put(os.path.join(install_file_dir, script), remote, mode=0777, use_sudo=True)
    # Create upstart configuration file for boot-time script
    cloudman_boot_file = 'cloudman.conf'
    remote_file = '/etc/init/%s' % cloudman_boot_file
    _write_to_file(cm_upstart % (remote, os.path.splitext(remote)[0]), remote_file, mode=0644)
    env.logger.debug("Done configuring CloudMan ec2_autorun")


def _configure_sge(env):
    """This method only sets up the environment for SGE w/o actually setting up SGE"""
    sge_root = '/opt/sge'
    if not exists(sge_root):
        sudo("mkdir -p %s" % sge_root)
        sudo("chown sgeadmin:sgeadmin %s" % sge_root)
    # Link our installed SGE to CloudMan's expected directory
    sge_package_dir = "/opt/galaxy/pkg"
    sge_dir = "ge6.2u5"
    if not exists(os.path.join(sge_package_dir, sge_dir)):
        sudo("mkdir -p %s" % sge_package_dir)
    if not exists(os.path.join(sge_package_dir, sge_dir)):
        sudo("ln --force -s %s/%s %s/%s" % (env.install_dir, sge_dir, sge_package_dir, sge_dir))
    env.logger.debug("Done configuring CloudMan SGE")

def _configure_hadoop(env):
    """
    Grab files required by CloudMan to setup a Hadoop cluster atop SGE.
    """
    hadoop_root = '/opt/hadoop'
    if not exists(hadoop_root):
        sudo("mkdir -p %s" % hadoop_root)
    with cd(hadoop_root):
        if not exists('hadoop.tar.gz'):
            sudo("wget --output-document=hadoop.tar.gz https://s3.amazonaws.com/cloudman/hadoop.tar.gz")
        if not exists('sge_integration.tar.gz'):
            sudo("wget --output-document=sge_integration.tar.gz https://s3.amazonaws.com/cloudman/sge_integration.tar.gz")
    sudo("chown -R ubuntu:ubuntu {0}".format(hadoop_root))
    env.logger.debug("Done configuring Hadoop for CloudMan")

def _configure_condor(env):
    """
    Grab files required by CloudMan to setup HTCondor
    """
    condor_root = '/opt/condor'
    filename = 'condor-7.9.1-x86_64_ubuntu_10.04.4-stripped.tar.gz'
    if not exists(condor_root):
        sudo("mkdir -p %s" % condor_root)
    with cd(condor_root):
        if not exists(filename):
            sudo('wget --output-document={0} https://s3.amazonaws.com/cloudman/{0}'\
                .format(filename))
    sudo("chown -R condor:condor {0}".format(condor_root))
    env.logger.debug("Done configuring HTCondor for CloudMan")

def _configure_nfs(env):
    nfs_dir = "/export/data"
    cloudman_dir = "/mnt/galaxyData/export"
    if not exists(nfs_dir):
        # For the case of rerunning this script, ensure the nfs_dir does
        # not exist (exists() method does not recognize it as a file because
        # by default it points to a non-existing dir/file).
        with settings(warn_only=True):
            sudo('rm -rf {0}'.format(nfs_dir))
        sudo("mkdir -p %s" % os.path.dirname(nfs_dir))
        sudo("ln -s %s %s" % (cloudman_dir, nfs_dir))
    sudo("chown -R %s %s" % (env.user, os.path.dirname(nfs_dir)))
    # Setup /etc/exports paths, to be used as NFS mount points
    galaxy_data_mount = env.get("galaxy_data_mount", "/mnt/galaxyData")
    galaxy_indices_mount = env.get("galaxy_indices_mount", "/mnt/galaxyIndices")
    galaxy_tools_mount = env.get("galaxy_tools_mount", "/mnt/galaxyTools")
    exports = [ '/opt/sge           *(rw,sync,no_root_squash,no_subtree_check)',
                '/opt/condor           *(rw,sync,no_root_squash,no_subtree_check)',
                '%s    *(rw,sync,no_root_squash,subtree_check,no_wdelay)' % galaxy_data_mount,
                '%s *(rw,sync,no_root_squash,no_subtree_check)' % galaxy_indices_mount,
                '%s   *(rw,sync,no_root_squash,no_subtree_check)' % galaxy_tools_mount,
                '%s       *(rw,sync,no_root_squash,no_subtree_check)' % nfs_dir,
                '%s/openmpi         *(rw,sync,no_root_squash,no_subtree_check)' % env.install_dir]
    extra_nfs_exports = env.get("extra_nfs_exports", "")
    for extra_nfs_export in extra_nfs_exports.split(","):
        exports.append('%s   *(rw,sync,no_root_squash,no_subtree_check)' % extra_nfs_export)
    append('/etc/exports', exports, use_sudo=True)
    # Create a symlink for backward compatibility where all of CloudMan's
    # stuff is expected to be in /opt/galaxy
    old_dir = '/opt/galaxy'
    # Because stow is used, the equivalent to CloudMan's expected path
    # is actually the parent of the install_dir so use it for the symlink
    new_dir = os.path.dirname(env.install_dir)
    if not exists(old_dir) and exists(new_dir):
        sudo('ln -s {0} {1}'.format(new_dir, old_dir))
    env.logger.debug("Done configuring CloudMan NFS")

@_if_not_installed("s3fs")
def install_s3fs(env):
    """
    Install s3fs, allowing S3 buckets to be mounted as ~POSIX file systems
    """
    default_version = "1.61"
    version = env.get("tool_version", default_version)
    url = "http://s3fs.googlecode.com/files/s3fs-%s.tar.gz" % version
    _get_install(url, env, _configure_make)

def _cleanup_ec2(env):
    """Clean up any extra files after building.
    """
    env.logger.info("Cleaning up for EC2 AMI creation")
    fnames = [".bash_history", "/var/log/firstboot.done", ".nx_setup_done",
              "/var/crash/*", "%s/ec2autorun.py.log" % env.install_dir,
              "%s/ec2autorun.err"  % env.install_dir, "%s/ec2autorun.log" % env.install_dir,
              "%s/bin/ec2autorun.log" % env.install_dir]
    for fname in fnames:
        sudo("rm -f %s" % fname)
    rmdirs = ["/mnt/galaxyData", "/mnt/cm", "/tmp/cm"]
    for rmdir in rmdirs:
        sudo("rm -rf %s" % rmdir)
    # Stop Apache from starting automatically at boot (it conflicts with Galaxy's nginx)
    sudo('/usr/sbin/update-rc.d -f apache2 remove')
    with settings(warn_only=True):
        # RabbitMQ fails to start if its database is embedded into the image
        # because it saves the current IP address or host name so delete it now.
        # When starting up, RabbitMQ will recreate that directory.
        sudo('/etc/init.d/rabbitmq-server stop')
        sudo('service rabbitmq-server stop')
        # Clean up packages that are causing issues or are unnecessary
        pkgs_to_remove = ['tntnet', 'tntnet-runtime', 'libtntnet9']
        for ptr in pkgs_to_remove:
            sudo('apt-get -y --force-yes remove --purge {0}'.format(ptr))
    sudo('initctl reload-configuration')
    for db_location in ['/var/lib/rabbitmq/mnesia', '/mnesia']:
        if exists(db_location):
            sudo('rm -rf %s' % db_location)
    # remove existing ssh host key pairs
    # http://docs.amazonwebservices.com/AWSEC2/latest/UserGuide/AESDG-chapter-sharingamis.html
    sudo("rm -f /etc/ssh/ssh_host_*")
    sudo("rm -f ~/.ssh/authorized_keys*")
    sudo("rm -f /root/.ssh/authorized_keys*")
