"""Custom install scripts for CloudMan environment.

From Enis Afgan: https://bitbucket.org/afgane/mi-deployment
"""
import os
import contextlib

from fabric.api import cd
from fabric.contrib.files import settings, hide

from cloudbio.custom.shared import (_make_tmp_dir, _setup_conf_file)
from cloudbio.cloudman import (_configure_cloudman, _configure_novnc,
    _configure_desktop, _configure_ec2_autorun)
from cloudbio.galaxy import _install_nginx

CDN_ROOT_URL = "http://linuxcourse.rutgers.edu/rate/Clusters/download"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"


def install_cloudman(env):
    """ A meta method for installing all of CloudMan components.
        Allows CloudMan and all of its dependencies to be installed via:
        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman
    """
    env.logger.debug("Installing CloudMan")
    _configure_cloudman(env, use_repo_autorun=False)
    install_nginx(env)
    install_proftpd(env)
    install_sge(env)
    install_novnc(env)


def install_ec2_autorun(env):
    _configure_ec2_autorun(env)


def install_novnc(env):
    _configure_novnc(env)
    _configure_desktop(env)


def install_nginx(env):
    _install_nginx(env)


def install_proftpd(env):
    """Highly configurable GPL-licensed FTP server software.
    http://proftpd.org/
    """
    version = "1.3.4c"
    postgres_ver = "9.1"
    url = "ftp://ftp.tpnet.pl/pub/linux/proftpd/distrib/source/proftpd-%s.tar.gz" % version
    modules = "mod_sql:mod_sql_postgres:mod_sql_passwd"
    extra_modules = env.get("extra_proftp_modules", "")  # Comma separated list of extra modules
    if extra_modules:
        modules = "%s:%s" % (modules, extra_modules.replace(",", ":"))
    install_dir = os.path.join(env.install_dir, 'proftpd')
    remote_conf_dir = os.path.join(install_dir, "etc")
    # Skip install if already available
    if env.safe_exists(remote_conf_dir):
        env.logger.debug("ProFTPd seems to already be installed in {0}".format(install_dir))
        return
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            env.safe_run("wget %s" % url)
            with settings(hide('stdout')):
                env.safe_run("tar xvzf %s" % os.path.split(url)[1])
            with cd("proftpd-%s" % version):
                env.safe_run("CFLAGS='-I/usr/include/postgresql' ./configure --prefix=%s "
                    "--disable-auth-file --disable-ncurses --disable-ident --disable-shadow "
                    "--enable-openssl --with-modules=%s "
                    "--with-libraries=/usr/lib/postgresql/%s/lib" % (install_dir, modules, postgres_ver))
                env.safe_sudo("make")
                env.safe_sudo("make install")
                env.safe_sudo("make clean")
    # Get the init.d startup script
    initd_script = 'proftpd.initd'
    initd_url = os.path.join(REPO_ROOT_URL, 'conf_files', initd_script)
    remote_file = "/etc/init.d/proftpd"
    env.safe_sudo("wget --output-document=%s %s" % (remote_file, initd_url))
    env.safe_sed(remote_file, 'REPLACE_THIS_WITH_CUSTOM_INSTALL_DIR', install_dir, use_sudo=True)
    env.safe_sudo("chmod 755 %s" % remote_file)
    # Set the configuration file
    conf_file = 'proftpd.conf'
    remote_file = os.path.join(remote_conf_dir, conf_file)
    if "postgres_port" not in env:
        env.postgres_port = '5910'
    if "galaxy_ftp_user_password" not in env:
        env.galaxy_ftp_user_password = 'fu5yOj2sn'
    proftpd_conf = {'galaxy_uid': env.safe_run('id -u galaxy'),
                    'galaxy_fs': '/mnt/galaxy',  # Should be a var but uncertain how to get it
                    'install_dir': install_dir}
    _setup_conf_file(env, remote_file, conf_file, overrides=proftpd_conf,
        default_source="proftpd.conf.template")
    # Get the custom welcome msg file
    welcome_msg_file = 'welcome_msg.txt'
    welcome_url = os.path.join(REPO_ROOT_URL, 'conf_files', welcome_msg_file)
    env.safe_sudo("wget --output-document=%s %s" %
       (os.path.join(remote_conf_dir, welcome_msg_file), welcome_url))
    # Stow
    env.safe_sudo("cd %s; stow proftpd" % env.install_dir)
    env.logger.debug("----- ProFTPd %s installed to %s -----" % (version, install_dir))


def install_sge(env):
    """Sun Grid Engine.
    """
    out_dir = "ge6.2u5"
    url = "%s/ge62u5_lx24-amd64.tar.gz" % CDN_ROOT_URL
    install_dir = env.install_dir
    if env.safe_exists(os.path.join(install_dir, out_dir)):
        return
    with _make_tmp_dir() as work_dir:
        with contextlib.nested(cd(work_dir), settings(hide('stdout'))):
            env.safe_run("wget %s" % url)
            env.safe_sudo("chown %s %s" % (env.user, install_dir))
            env.safe_run("tar -C %s -xvzf %s" % (install_dir, os.path.split(url)[1]))
    env.logger.debug("SGE setup")
