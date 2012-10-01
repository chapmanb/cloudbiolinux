from fabric.api import local, sudo, env, put, get
from fabric.contrib.files import exists, append

import os


def setup_install_dir():
    """Sets up install dir and ensures its owned by Galaxy"""
    if not exists(env.install_dir):
        sudo("mkdir -p %s" % env.install_dir)
    if not exists(env.jars_dir):
        sudo("mkdir -p %s" % env.jars_dir)
    chown_galaxy(os.path.split(env.install_dir)[0])


def ensure_can_sudo_into(user):
    sudoers_append("%admin  ALL = (" + user + ") NOPASSWD: ALL")


def sudoers_append(line):
    append("/etc/sudoers", line, use_sudo=True)


def start_service(service_name):
    # For reasons I don't understand this doesn't work for galaxy init
    # script unless pty=False
    sudo("/etc/init.d/%s start" % service_name, pty=False)


def wget(url, install_command=sudo, file_name=None):
    if not file_name:
        file_name = os.path.split(url)[-1]
        if '?' in file_name:
            file_name = file_name[0:file_name.index('?')]
    if ("cache_source_downloads" in env) and (not env.cache_source_downloads):
        install_command("wget %s -O %s" % (url, file_name))
    else:
        cache_dir = env.source_cache_dir
        if not cache_dir:
            cache_dir = ".downloads"
        cached_file = os.path.join(cache_dir, file_name)
        if os.path.exists(cached_file):
            put(cached_file, file_name)
        else:
            install_command("wget %s -O %s" % (url, file_name))
            local("mkdir -p '%s'" % cache_dir)
            get(file_name, cached_file)
