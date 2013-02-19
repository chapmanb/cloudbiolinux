"""Custom install scripts for CloudMan environment.

From Enis Afgan: https://bitbucket.org/afgane/mi-deployment
"""
import os
import contextlib

from fabric.api import sudo, run, cd
from fabric.contrib.files import exists, settings, hide, sed

from cloudbio.custom.shared import _make_tmp_dir
from cloudbio.cloudman import _configure_cloudman, _configure_novnc
from cloudbio.galaxy import _install_nginx

CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"


def install_cloudman(env):
    """ A meta method for installing all of CloudMan components.
        Allows CloudMan and all of its dependencies to be installed via:
        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman
    """
    _configure_cloudman(env, use_repo_autorun=False)
    install_nginx(env)
    install_proftpd(env)
    install_sge(env)

def install_novnc(env):
    _configure_novnc(env)

def install_nginx(env):
    _install_nginx(env)

def install_proftpd(env):
    """Highly configurable GPL-licensed FTP server software.
    http://proftpd.org/
    """
    version = "1.3.4a"
    postgres_ver = "9.1"
    url = "ftp://ftp.tpnet.pl/pub/linux/proftpd/distrib/source/proftpd-%s.tar.gz" % version
    install_dir = os.path.join(env.install_dir, 'proftpd')
    remote_conf_dir = os.path.join(install_dir, "etc")
    # skip install if already present
    if exists(remote_conf_dir):
        env.logger.debug("ProFTPd seems to already be installed in {0}".format(install_dir))
        return
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            with settings(hide('stdout')):
                run("tar xvzf %s" % os.path.split(url)[1])
            with cd("proftpd-%s" % version):
                run("CFLAGS='-I/usr/include/postgresql' ./configure --prefix=%s " \
                    "--disable-auth-file --disable-ncurses --disable-ident --disable-shadow " \
                    "--enable-openssl --with-modules=mod_sql:mod_sql_postgres:mod_sql_passwd " \
                    "--with-libraries=/usr/lib/postgresql/%s/lib" % (install_dir, postgres_ver))
                sudo("make")
                sudo("make install")
                sudo("make clean")
    # Get the init.d startup script
    initd_script = 'proftpd.initd'
    initd_url = os.path.join(REPO_ROOT_URL, 'conf_files', initd_script)
    remote_file = "/etc/init.d/proftpd"
    sudo("wget --output-document=%s %s" % (remote_file, initd_url))
    sed(remote_file, 'REPLACE_THIS_WITH_CUSTOM_INSTALL_DIR', install_dir, use_sudo=True)
    sudo("chmod 755 %s" % remote_file)
    # Set the configuration file
    conf_file = 'proftpd.conf'
    conf_url = os.path.join(REPO_ROOT_URL, 'conf_files', conf_file)
    remote_file = os.path.join(remote_conf_dir, conf_file)
    sudo("wget --output-document=%s %s" % (remote_file, conf_url))
    sed(remote_file, 'REPLACE_THIS_WITH_CUSTOM_INSTALL_DIR', install_dir, use_sudo=True)
    # Get the custom welcome msg file
    welcome_msg_file = 'welcome_msg.txt'
    welcome_url = os.path.join(REPO_ROOT_URL, 'conf_files', welcome_msg_file)
    sudo("wget --output-document=%s %s" % (os.path.join(remote_conf_dir, welcome_msg_file), welcome_url))
    # Stow
    sudo("cd %s; stow proftpd" % env.install_dir)
    env.logger.debug("----- ProFTPd %s installed to %s -----" % (version, install_dir))

def install_sge(env):
    """Sun Grid Engine.
    """
    out_dir = "ge6.2u5"
    url = "%s/ge62u5_lx24-amd64.tar.gz" % CDN_ROOT_URL
    install_dir = env.install_dir
    if exists(os.path.join(install_dir, out_dir)):
        return
    with _make_tmp_dir() as work_dir:
        with contextlib.nested(cd(work_dir), settings(hide('stdout'))):
            run("wget %s" % url)
            sudo("chown %s %s" % (env.user, install_dir))
            run("tar -C %s -xvzf %s" % (install_dir, os.path.split(url)[1]))
    env.logger.debug("SGE setup")
