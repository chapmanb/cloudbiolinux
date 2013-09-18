"""
Install system programs not available from packages.
"""
import os

from shared import _if_not_installed, _get_install, _configure_make

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
        outdir = os.path.join(env.local_install, "homebrew")
        if env.safe_exists(os.path.join(outdir, "bin", "brew")):
            return
        env.safe_sudo("chown %s %s" % (env.user, env.local_install))
        env.safe_run("git clone https://github.com/Homebrew/linuxbrew.git %s" % outdir)
        comment_line = "# Linuxbrew -- added by CloudBioLinux"
        if not env.safe_contains(env.shell_config, comment_line):
            env.safe_append(env.shell_config, comment_line)
            env.safe_append(env.shell_config, "export PATH=%s/bin:$PATH" % outdir)
            env.safe_append(env.shell_config, "export LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH" % outdir)

@_if_not_installed("s3fs")
def install_s3fs(env):
    """FUSE-based file system backed by Amazon S3.
    https://code.google.com/p/s3fs/
    """
    version = "1.61"
    url = "http://s3fs.googlecode.com/files/s3fs-{0}.tar.gz".format(version)
    _get_install(url, env, _configure_make)
