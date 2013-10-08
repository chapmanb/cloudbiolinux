"""
Install system programs not available from packages.
"""
import os

from fabric.api import cd

from cloudbio.custom import shared
from cloudbio.custom.shared import _if_not_installed, _get_install, _configure_make

@_if_not_installed("brew")
def install_homebrew(env):
    """Homebrew package manager for OSX and Linuxbrew for linux systems.

    https://github.com/mxcl/homebrew
    https://github.com/Homebrew/linuxbrew
    """
    if env.distribution == "macosx":
        # XXX Test homebrew install on mac
        env.safe_run('ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go)"')
    else:
        brew_cmd = os.path.join(env.system_install, "bin", "brew")
        if not env.safe_exists(brew_cmd):
            with shared._make_tmp_dir() as tmp_dir:
                with cd(tmp_dir):
                    env.safe_run("git clone https://github.com/Homebrew/linuxbrew.git" )
                    with cd("linuxbrew"):
                        env.safe_sudo("chown %s %s" % (env.user, env.system_install))
                        paths = ["bin", "etc", "include", "lib", "lib/pkgconfig", "Library",
                                 "sbin", "share", "var", "var/log", "share/locale",
                                 "share/man", "share/man/man1", "share/man/man2",
                                 "share/man/man3", "share/man/man4", "share/man/man5",
                                 "share/man/man6", "share/man/man7", "share/man/man8",
                                 "share/info", "share/doc", "share/aclocal"]
                        for path in paths:
                            if env.safe_exists("%s/%s" % (env.system_install, path)):
                                env.safe_sudo("chown %s %s/%s" % (env.user, env.system_install, path))
                        env.safe_run("mv bin/brew %s/bin" % env.system_install)
                        env.safe_run("mv Library %s" % env.system_install)
                        env.safe_run("mv .git %s" % env.system_install)

@_if_not_installed("s3fs")
def install_s3fs(env):
    """FUSE-based file system backed by Amazon S3.
    https://code.google.com/p/s3fs/
    """
    version = "1.61"
    url = "http://s3fs.googlecode.com/files/s3fs-{0}.tar.gz".format(version)
    _get_install(url, env, _configure_make)
