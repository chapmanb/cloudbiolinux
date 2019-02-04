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

from fabric.api import sudo, cd, run, put
from fabric.contrib.files import exists, settings

from cloudbio.galaxy import _setup_users
from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.custom.shared import (_make_tmp_dir, _write_to_file, _get_install,
                                    _configure_make, _if_not_installed,
                                    _setup_conf_file, _add_to_profiles,
                                    _create_python_virtualenv,
                                    _setup_simple_service,
                                    _read_boolean)
from cloudbio.package.deb import (_apt_packages, _setup_apt_automation)

MI_REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"
CM_REPO_ROOT_URL = "https://bitbucket.org/galaxy/cloudman/raw/tip"


def _configure_cloudman(env, use_repo_autorun=False):
    """
    Configure the machine to be capable of running CloudMan.

    ..Also see: ``custom/cloudman.py``
    """
    env.logger.debug("Configuring CloudMan")
    _setup_users(env)
    _setup_env(env)
    _configure_logrotate(env)
    _configure_ec2_autorun(env, use_repo_autorun)
    _configure_sge(env)
    _configure_hadoop(env)
    _configure_nfs(env)
    _configure_novnc(env)
    _configure_desktop(env)
    install_s3fs(env)


def _configure_desktop(env):
    """
    Configure a desktop manager to work with VNC. Note that `xfce4` (or `jwm`)
    and `vnc4server` packages need to be installed for this to have effect.
    """
    if not _read_boolean(env, "configure_desktop", False):
        return
    # Set nginx PAM module to allow logins for any system user
    if env.safe_exists("/etc/pam.d"):
        env.safe_sudo('echo "@include common-auth" > /etc/pam.d/nginx')
    env.safe_sudo('usermod -a -G shadow galaxy')
    # Create a start script for X
    _setup_conf_file(env, "/home/ubuntu/.vnc/xstartup", "xstartup", default_source="xstartup")
    # Create jwmrc config file (uncomment this if using jwm window manager)
    # _setup_conf_file(env, "/home/ubuntu/.jwmrc", "jwmrc.xml",
    #     default_source="jwmrc.xml", mode="0644")
    env.logger.info("----- Done configuring desktop -----")


def _configure_novnc(env):
    if not _read_boolean(env, "configure_novnc", False):
        # Longer term would like this enabled by default. -John
        return
    if not "novnc_install_dir" in env:
        env.novnc_install_dir = "/opt/novnc"
    if not "vnc_password" in env:
        env.vnc_password = "cl0udbi0l1nux"
    if not "vnc_user" in env:
        env.vnc_user = env.user
    if not "vnc_display" in env:
        env.vnc_display = "1"
    if not "vnc_depth" in env:
        env.vnc_depth = "16"
    if not "vnc_geometry" in env:
        env.vnc_geometry = "1024x768"

    _configure_vncpasswd(env)

    novnc_dir = env.novnc_install_dir
    env.safe_sudo("mkdir -p '%s'" % novnc_dir)
    env.safe_sudo("chown %s '%s'" % (env.user, novnc_dir))
    clone_cmd = "NOVNC_DIR='%s'; rm -rf $NOVNC_DIR; git clone https://github.com/kanaka/noVNC.git $NOVNC_DIR" % novnc_dir
    run(clone_cmd)
    ## Move vnc_auto.html which takes vnc_password as query argument
    ## to index.html and rewrite it so that password is autoset, no
    ## need to specify via query parameter.
    run("sed s/password\\ =\\ /password\\ =\\ \\\'%s\\\'\\;\\\\\\\\/\\\\\\\\// '%s/vnc_auto.html' > '%s/index.html'" % (env.vnc_password, novnc_dir, novnc_dir))

    _setup_conf_file(env, "/etc/init.d/novnc", "novnc_init", default_source="novnc_init")
    _setup_conf_file(env, "/etc/default/novnc", "novnc_default", default_source="novnc_default.template")
    _setup_conf_file(env, "/etc/init.d/vncserver", "vncserver_init", default_source="vncserver_init")
    _setup_conf_file(env, "/etc/default/vncserver", "vncserver_default", default_source="vncserver_default.template")
    _setup_simple_service("novnc")
    _setup_simple_service("vncserver")


def _configure_vncpasswd(env):
    with cd("~"):
        run("mkdir -p ~/.vnc")
        run("rm -rf vncpasswd")
        run("git clone https://github.com/trinitronx/vncpasswd.py vncpasswd")
        run("python vncpasswd/vncpasswd.py '%s' -f ~/.vnc/passwd" % env.vnc_password)
        run("chmod 600 ~/.vnc/passwd")
        run("rm -rf vncpasswd")


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
        # Allow flavors to modify the package list
        packages = env.flavor.rewrite_config_items("packages", packages)
        _setup_apt_automation()
        _apt_packages(pkg_list=packages)
    elif env.distribution in ["centos", "scientificlinux"]:
        env.logger.warn("No CloudMan system package dependencies for CentOS")
        pass
    # Get and install required Python libraries
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            url = os.path.join(CM_REPO_ROOT_URL, 'requirements.txt')
            _create_python_virtualenv(env, 'CM', reqs_url=url)
    # Add a custom vimrc
    vimrc_url = os.path.join(MI_REPO_ROOT_URL, 'conf_files', 'vimrc')
    remote_file = '/etc/vim/vimrc'
    if env.safe_exists("/etc/vim"):
        env.safe_sudo("wget --output-document=%s %s" % (remote_file, vimrc_url))
        env.logger.debug("Added a custom vimrc to {0}".format(remote_file))
    # Setup profile
    aliases = ['alias lt="ls -ltr"', 'alias ll="ls -l"']
    for alias in aliases:
        _add_to_profiles(alias, ['/etc/bash.bashrc'])
    env.logger.info("Done setting up CloudMan's environment")


def _configure_logrotate(env):
    """
    Add logrotate config file, which will automatically rotate CloudMan's log
    """
    conf_file = "cloudman.logrotate"
    remote = '/etc/logrotate.d/cloudman'
    url = os.path.join(MI_REPO_ROOT_URL, 'conf_files', conf_file)
    env.safe_sudo("wget --output-document=%s %s" % (remote, url))
    env.logger.info("----- Added logrotate file to {0} -----".format(remote))


def _configure_ec2_autorun(env, use_repo_autorun=False):
    """
    ec2autorun.py is a script that launches CloudMan on instance boot
    and is thus required on an instance. See the script itself for the
    details of what it does.

    This script also adds a cloudman service to ``/etc/init``, which
    actually runs ec2autorun.py as a system-level service at system boot.
    """
    script = "ec2autorun.py"
    remote = os.path.join(env.install_dir, "bin", script)
    if not env.safe_exists(os.path.dirname(remote)):
        env.safe_sudo('mkdir -p {0}'.format(os.path.dirname(remote)))
    if use_repo_autorun:
        # Is this used, can we eliminate use_repo_autorun?
        url = os.path.join(MI_REPO_ROOT_URL, script)
        env.safe_sudo("wget --output-document=%s %s" % (remote, url))
    else:
        install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
        tmp_remote = os.path.join("/tmp", os.path.basename(remote))
        env.safe_put(os.path.join(install_file_dir, script), tmp_remote)
        env.safe_sudo("mv %s %s" % (tmp_remote, remote))
        env.safe_sudo("chmod 0777 %s" % remote)
    # Create upstart configuration file for boot-time script
    cloudman_boot_file = 'cloudman.conf'
    remote_file = '/etc/init/%s' % cloudman_boot_file
    _write_to_file(cm_upstart % (remote, os.path.splitext(remote)[0]), remote_file, mode="0644")
    # Setup default image user data (if configured by image_user_data_path or
    # image_user_data_template_path). This specifies defaults for CloudMan when
    # used with resulting image, normal userdata supplied by user will override
    # these defaults.
    image_user_data_path = os.path.join(env.install_dir, "bin", "IMAGE_USER_DATA")
    if "image_user_data_dict" in env:
        # Explicit YAML contents defined in env, just dump them as is.
        import yaml
        _write_to_file(yaml.dump(env.get("image_user_data_dict")), image_user_data_path, mode="0644")
    else:
        # Else use file or template file.
        _setup_conf_file(env, image_user_data_path, "image_user_data", default_source="image_user_data")
    env.logger.info("Done configuring CloudMan's ec2_autorun")


def _configure_sge(env):
    """
    This method sets up the environment for SGE w/o
    actually setting up SGE; it basically makes sure system paths expected
    by CloudMan exist on the system.

    TODO: Merge this with ``install_sge`` method in ``custom/cloudman.py``.
    """
    sge_root = '/opt/sge'
    if not env.safe_exists(sge_root):
        env.safe_sudo("mkdir -p %s" % sge_root)
        env.safe_sudo("chown sgeadmin:sgeadmin %s" % sge_root)
    # Link our installed SGE to CloudMan's expected directory
    sge_package_dir = "/opt/galaxy/pkg"
    sge_dir = "ge6.2u5"
    if not env.safe_exists(os.path.join(sge_package_dir, sge_dir)):
        env.safe_sudo("mkdir -p %s" % sge_package_dir)
    if not env.safe_exists(os.path.join(sge_package_dir, sge_dir)):
        env.safe_sudo("ln --force -s %s/%s %s/%s" % (env.install_dir, sge_dir, sge_package_dir, sge_dir))
    env.logger.info("Done configuring SGE for CloudMan")


def _configure_hadoop(env):
    """
    Grab files required by CloudMan to setup a Hadoop cluster atop SGE.
    """
    hadoop_root = '/opt/hadoop'
    url_root = 'https://s3.amazonaws.com/cloudman'
    hcm_file = 'hadoop.1.0.4__1.0.tar.gz'
    si_file = 'sge_integration.1.0.tar.gz'
    # Make sure we're working with a clean hadoop_home dir to avoid any version conflicts
    env.safe_sudo("rm -rf {0}".format(hadoop_root))
    env.safe_sudo("mkdir -p %s" % hadoop_root)
    with cd(hadoop_root):
        env.safe_sudo("wget --output-document={0} {1}/{0}".format(hcm_file, url_root))
        env.safe_sudo("wget --output-document={0} {1}/{0}".format(si_file, url_root))
    env.safe_sudo("chown -R {0} {1}".format(env.user, hadoop_root))
    env.logger.info("Done configuring Hadoop for CloudMan")


def _configure_nfs(env):
    """
    Edit ``/etc/exports`` to append paths that are shared over NFS by CloudMan.

    In addition to the hard coded paths listed here, additional paths
    can be included by setting ``extra_nfs_exports`` in ``fabricrc.txt`` as
    a comma-separated list of directories.
    """
    nfs_dir = "/export/data"
    cloudman_dir = "/mnt/galaxy/export"
    if not env.safe_exists(nfs_dir):
        # For the case of rerunning this script, ensure the nfs_dir does
        # not exist (exists() method does not recognize it as a file because
        # by default it points to a non-existing dir/file).
        with settings(warn_only=True):
            env.safe_sudo('rm -rf {0}'.format(nfs_dir))
        env.safe_sudo("mkdir -p %s" % os.path.dirname(nfs_dir))
        env.safe_sudo("ln -s %s %s" % (cloudman_dir, nfs_dir))
    env.safe_sudo("chown -R %s %s" % (env.user, os.path.dirname(nfs_dir)))
    # Setup /etc/exports paths, to be used as NFS mount points
    # galaxy_data_mount = env.get("galaxy_data_mount", "/mnt/galaxyData")
    # galaxy_indices_mount = env.get("galaxy_indices_mount", "/mnt/galaxyIndices")
    # galaxy_tools_mount = env.get("galaxy_tools_mount", "/mnt/galaxyTools")
    exports = ['/opt/sge           *(rw,sync,no_root_squash,no_subtree_check)',
               '/opt/hadoop           *(rw,sync,no_root_squash,no_subtree_check)',
               # '%s    *(rw,sync,no_root_squash,subtree_check,no_wdelay)' % galaxy_data_mount,
               # '%s *(rw,sync,no_root_squash,no_subtree_check)' % galaxy_indices_mount,
               # '%s   *(rw,sync,no_root_squash,no_subtree_check)' % galaxy_tools_mount,
               # '%s       *(rw,sync,no_root_squash,no_subtree_check)' % nfs_dir,
               # '%s/openmpi         *(rw,sync,no_root_squash,no_subtree_check)' % env.install_dir
               ]
    extra_nfs_exports = env.get("extra_nfs_exports", "")
    if extra_nfs_exports:
        for extra_nfs_export in extra_nfs_exports.split(","):
            exports.append('%s   *(rw,sync,no_root_squash,no_subtree_check)' % extra_nfs_export)
    env.safe_append('/etc/exports', exports, use_sudo=True)
    # Create a symlink for backward compatibility where all of CloudMan's
    # stuff is expected to be in /opt/galaxy
    old_dir = '/opt/galaxy'
    # Because stow is used, the equivalent to CloudMan's expected path
    # is actually the parent of the install_dir so use it for the symlink
    new_dir = os.path.dirname(env.install_dir)
    if not env.safe_exists(old_dir) and exists(new_dir):
        env.safe_sudo('ln -s {0} {1}'.format(new_dir, old_dir))
    env.logger.info("Done configuring NFS for CloudMan")


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
    """
    Clean up any extra files after building. This method must be called
    on an instance after being built and before creating a new machine
    image. *Note* that after this method has run, key-based ssh access
    to the machine is no longer possible.
    """
    env.logger.info("Cleaning up for EC2 AMI creation")
    # Clean up log files and such
    fnames = [".bash_history", "/var/log/firstboot.done", ".nx_setup_done",
              "/var/crash/*", "%s/ec2autorun.py.log" % env.install_dir,
              "%s/ec2autorun.err" % env.install_dir, "%s/ec2autorun.log" % env.install_dir,
              "%s/bin/ec2autorun.log" % env.install_dir]
    for fname in fnames:
        sudo("rm -f %s" % fname)
    rmdirs = ["/mnt/galaxyData", "/mnt/cm", "/tmp/cm"]
    for rmdir in rmdirs:
        sudo("rm -rf %s" % rmdir)
    # Seed the history with frequently used commands
    env.logger.debug("Setting bash history")
    local = os.path.join(env.config_dir, os.pardir, "installed_files", "bash_history")
    remote = os.path.join('/home', 'ubuntu', '.bash_history')
    put(local, remote, mode="0660", use_sudo=True)
    # Make sure the default config dir is owned by ubuntu
    sudo("chown ubuntu:ubuntu ~/.config")
    # Stop Apache from starting automatically at boot (it conflicts with Galaxy's nginx)
    sudo('/usr/sbin/update-rc.d -f apache2 remove')
    with settings(warn_only=True):
        # RabbitMQ fails to start if its database is embedded into the image
        # because it saves the current IP address or host name so delete it now.
        # When starting up, RabbitMQ will recreate that directory.
        sudo('/etc/init.d/rabbitmq-server stop')
        sudo('service rabbitmq-server stop')
        # Clean up packages that are causing issues or are unnecessary
        pkgs_to_remove = ['tntnet', 'tntnet-runtime', 'libtntnet9', 'vsftpd']
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
