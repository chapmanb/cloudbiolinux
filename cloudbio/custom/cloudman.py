"""Custom install scripts for CloudMan environment.

From Enis Afgan: https://bitbucket.org/afgane/mi-deployment
"""
import os, contextlib

from fabric.api import sudo, run, cd 
from fabric.contrib.files import exists, settings, hide, contains, sed

from cloudbio.custom.shared import _make_tmp_dir
from cloudbio.cloudman import _configure_cloudman

CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"

def install_cloudman(env):
    """ A meta method for installing all of CloudMan components.
        Allows CloudMan and all of its dependencies to be installed via:
        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman
    """
    _configure_cloudman(env, use_repo_autorun=True)
    install_nginx(env)
    install_proftpd(env)
    install_sge(env)

def install_nginx(env):
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

    nginx_conf_file = 'nginx.conf'
    url = os.path.join(REPO_ROOT_URL, nginx_conf_file)
    with cd(remote_conf_dir):
        sudo("wget --output-document=%s/%s %s" % (remote_conf_dir, nginx_conf_file, url))

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
    upload_module_version = "2.2.0"
    chunk_module_version = "0.22"
    chunk_git_version = "b46dd27"
    modules = []
    upload_url = "http://www.grid.net.ru/nginx/download/" \
                 "nginx_upload_module-%s.tar.gz" % upload_module_version
    run("wget %s" % upload_url)
    upload_fname = os.path.split(upload_url)[1]
    modules.append(upload_fname.rsplit(".", 2)[0])
    run("tar -xvzpf %s" % upload_fname)
    chunk_url = "https://github.com/agentzh/chunkin-nginx-module/tarball/v%s" % chunk_module_version
    chunk_fname = "agentzh-chunkin-nginx-module-%s.tar.gz" % (chunk_git_version)
    run("wget -O %s %s" % (chunk_fname, chunk_url))
    run("tar -xvzpf %s" % chunk_fname)
    modules.append(chunk_fname.rsplit(".", 2)[0])
    return modules

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
