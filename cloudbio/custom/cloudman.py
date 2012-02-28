"""Custom install scripts for CloudMan environment.

From Enis Afgan: https://bitbucket.org/afgane/mi-deployment
"""
import os, contextlib

from fabric.api import sudo, run, env, cd, put, local
from fabric.contrib.console import confirm
from fabric.contrib.files import exists, settings, hide, contains, append, sed

from cloudbio.custom.shared import _make_tmp_dir

CDN_ROOT_URL = "http://userwww.service.emory.edu/~eafgan/content"
REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"

def install_nginx(env):
    version = "0.7.67"
    url = "http://nginx.org/download/nginx-%s.tar.gz" % version

    install_dir = os.path.join(env.install_dir, "nginx")
    remote_conf_dir = os.path.join(install_dir, "conf")

    # skip install if already present
    if exists(remote_conf_dir) and contains(os.path.join(remote_conf_dir, "nginx.conf"), "/cloud"):
        return

    with _make_tmp_dir() as work_dir:
        with contextlib.nested(cd(work_dir), settings(hide('stdout'))):
            modules = _get_nginx_modules(env)
            module_flags = " ".join(["--add-module=../%s" % x for x in modules])
            run("wget %s" % url)
            run("tar xvzf %s" % os.path.split(url)[1])
            with cd("nginx-%s" % version):
                run("./configure --prefix=%s --with-ipv6 %s "
                    "--user=galaxy --group=galaxy "
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

def _get_nginx_modules(env):
    """Retrieve add-on modules compiled along with nginx.
    """
    upload_module_version = "2.0.12"
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
    version = "1.3.3d"
    postgres_ver = "8.4"
    url = "ftp://mirrors.ibiblio.org/proftpd/distrib/source/proftpd-%s.tar.gz" % version
    install_dir = os.path.join(env.install_dir, 'proftpd')
    remote_conf_dir = os.path.join(install_dir, "etc")
    # skip install if already present
    if exists(remote_conf_dir):
        return
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            with settings(hide('stdout')):
                run("tar xvzf %s" % os.path.split(url)[1])
            with cd("proftpd-%s" % version):
                run("CFLAGS='-I/usr/include/postgresql' ./configure --prefix=%s --disable-auth-file --disable-ncurses --disable-ident --disable-shadow --enable-openssl --with-modules=mod_sql:mod_sql_postgres:mod_sql_passwd --with-libraries=/usr/lib/postgres/%s/lib" % (install_dir, postgres_ver))
                sudo("make")
                sudo("make install")
                sudo("make clean")
                # Get init.d startup script
                initd_script = 'proftpd'
                initd_url = os.path.join(REPO_ROOT_URL, 'conf_files', initd_script)
                sudo("wget --output-document=%s %s" % (os.path.join('/etc/init.d', initd_script), initd_url))
                sudo("chmod 755 %s" % os.path.join('/etc/init.d', initd_script))
                # Get configuration files
                proftpd_conf_file = 'proftpd.conf'
                welcome_msg_file = 'welcome_msg.txt'
                conf_url = os.path.join(REPO_ROOT_URL, 'conf_files', proftpd_conf_file)
                welcome_url = os.path.join(REPO_ROOT_URL, 'conf_files', welcome_msg_file)
                sudo("wget --output-document=%s %s" % (os.path.join(remote_conf_dir, proftpd_conf_file), conf_url))
                sudo("wget --output-document=%s %s" % (os.path.join(remote_conf_dir, welcome_msg_file), welcome_url))
                sudo("cd %s; stow proftpd" % env.install_dir)

def install_sge(env):
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
