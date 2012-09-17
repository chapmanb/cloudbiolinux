"""Reusable decorators and functions for custom installations.
"""
import tempfile
import os
import functools
import urllib
from tempfile import NamedTemporaryFile
from string import Template
from contextlib import contextmanager

from fabric.api import *
from fabric.contrib.files import *

CBL_REPO_ROOT_URL = "https://raw.github.com/chapmanb/cloudbiolinux/master/"

# -- decorators and context managers

def _if_not_installed(pname):
    """Decorator that checks if a callable program is installed.
    """
    def argcatcher(func):
        functools.wraps(func)

        def decorator(*args, **kwargs):
            if _galaxy_tool_install(args):
                run_function = not _galaxy_tool_present(args)
            else:
                run_function = _executable_not_on_path(pname)

            if run_function:
                return func(*args, **kwargs)
        return decorator
    return argcatcher


def _executable_not_on_path(pname):
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):
        result = run(pname)
    return result.return_code == 127


def _galaxy_tool_install(args):
    try:
        return args[0]["galaxy_tool_install"]
    except:
        return False


def _galaxy_tool_present(args):
    return exists(os.path.join(args[0]["system_install"], "env.sh"))


def _if_not_python_lib(library):
    """Decorator that checks if a python library is installed.
    """
    def argcatcher(func):
        functools.wraps(func)
        def decorator(*args, **kwargs):
            with settings(warn_only=True):
                pyver = env.python_version_ext if env.has_key("python_version_ext") else ""
                result = run("python%s -c 'import %s'" % (pyver, library))
            if result.failed:
                return func(*args, **kwargs)
        return decorator
    return argcatcher

@contextmanager
def _make_tmp_dir():
    tmp_dir = run("echo $TMPDIR").strip()
    if not tmp_dir:
        home_dir = run("echo $HOME")
        tmp_dir = os.path.join(home_dir, "tmp")
    work_dir = os.path.join(tmp_dir, "cloudbiolinux")
    if not exists(work_dir):
        run("mkdir -p %s" % work_dir)
    yield work_dir
    if exists(work_dir):
        run("rm -rf %s" % work_dir)

# -- Standard build utility simplifiers

def _get_expected_file(url):
    tar_file = os.path.split(url.split("?")[0])[-1]
    safe_tar = "--pax-option='delete=SCHILY.*,delete=LIBARCHIVE.*'"
    exts = {(".tar.gz", ".tgz") : "tar %s -xzpf" % safe_tar,
            (".tar.bz2",): "tar %s -xjpf" % safe_tar,
            (".zip",) : "unzip"}
    for ext_choices, tar_cmd in exts.iteritems():
        for ext in ext_choices:
            if tar_file.endswith(ext):
                return tar_file, tar_file[:-len(ext)], tar_cmd
    raise ValueError("Did not find extract command for %s" % url)

def _safe_dir_name(dir_name, need_dir=True):
    replace_try = ["", "-src", "_core"]
    for replace in replace_try:
        check = dir_name.replace(replace, "")
        if exists(check):
            return check
    # still couldn't find it, it's a nasty one
    for check_part in (dir_name.split("-")[0].split("_")[0],
                       dir_name.split("-")[-1].split("_")[-1],
                       dir_name.split(".")[0]):
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            dirs = run("ls -d1 *%s*/" % check_part).split("\n")
            dirs = [x for x in dirs if "cannot access" not in x and "No such" not in x]
        if len(dirs) == 1:
            return dirs[0]
    if need_dir:
        raise ValueError("Could not find directory %s" % dir_name)

def _fetch_and_unpack(url, need_dir=True):
    if url.startswith(("git", "svn", "hg", "cvs")):
        run(url)
        base = os.path.basename(url.split()[-1])
        return os.path.splitext(base)[0]
    else:
        tar_file, dir_name, tar_cmd = _get_expected_file(url)
        if not exists(tar_file):
            run("wget --no-check-certificate -O %s '%s'" % (tar_file, url))
        run("%s %s" % (tar_cmd, tar_file))
        return _safe_dir_name(dir_name, need_dir)

def _configure_make(env):
    run("./configure --disable-werror --prefix=%s " % env.system_install)
    run("make")
    env.safe_sudo("make install")

def _make_copy(find_cmd=None, premake_cmd=None, do_make=True):
    def _do_work(env):
        if premake_cmd:
            premake_cmd()
        if do_make:
            run("make")
        if find_cmd:
            install_dir = _get_bin_dir(env)
            for fname in run(find_cmd).split("\n"):
                env.safe_sudo("mv -f %s %s" % (fname.rstrip("\r"), install_dir))
    return _do_work

def _get_install(url, env, make_command, post_unpack_fn=None):
    """Retrieve source from a URL and install in our system directory.
    """
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            dir_name = _fetch_and_unpack(url)
            with cd(dir_name):
                if post_unpack_fn:
                    post_unpack_fn(env)
                make_command(env)

def _get_install_local(url, env, make_command, dir_name=None):
    """Build and install in a local directory.
    """
    (_, test_name, _) = _get_expected_file(url)
    test1 = os.path.join(env.local_install, test_name)
    if dir_name is not None:
        test2 = os.path.join(env.local_install, dir_name)
    elif "-" in test1:
        test2, _ = test1.rsplit("-", 1)
    else:
        test2 = os.path.join(env.local_install, test_name.split("_")[0])
    if not exists(test1) and not exists(test2):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                dir_name = _fetch_and_unpack(url)
                if not exists(os.path.join(env.local_install, dir_name)):
                    with cd(dir_name):
                        make_command(env)
                    run("mv %s %s" % (dir_name, env.local_install))

# --- Language specific utilities

def _symlinked_java_version_dir(pname, version, env):
    base_dir = os.path.join(env.system_install, "share", "java", pname)
    install_dir = "%s-%s" % (base_dir, version)
    if not exists(install_dir):
        env.safe_sudo("mkdir -p %s" % install_dir)
        if exists(base_dir):
            env.safe_sudo("rm -f %s" % base_dir)
        env.safe_sudo("ln -s %s %s" % (install_dir, base_dir))
        return install_dir
    return None

def _java_install(pname, version, url, env, install_fn=None):
    install_dir = _symlinked_java_version_dir(pname, version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                dir_name = _fetch_and_unpack(url)
                with cd(dir_name):
                    if install_fn is not None:
                        install_fn(env, install_dir)
                    else:
                        env.safe_sudo("mv *.jar %s" % install_dir)

def _python_make(env):
    run("python%s setup.py build" % env.python_version_ext)
    # handle standard Ubuntu case: only specify a prefix if we're not
    # in the standard system directory
    if env.system_install in ["/usr"]:
        prefix = ""
    else:
        prefix = "--prefix '%s'" % env.system_install
    env.safe_sudo("python%s setup.py install --skip-build %s" % (env.python_version_ext, prefix))
    for clean in ["dist", "build", "lib/*.egg-info"]:
        env.safe_sudo("rm -rf %s" % clean)


def _get_installed_file(env, local_file):
    path = os.path.join('installed_files', local_file)
    if not os.path.exists(path):
        # If using cloudbiolinux as a library, this won't be available,
        # download the file from github instead
        f = NamedTemporaryFile(delete=False)
        cloudbiolinx_repo_url = env.get("cloudbiolinux_repo_url", CBL_REPO_ROOT_URL)
        url = os.path.join(cloudbiolinx_repo_url, 'installed_files', local_file)
        urllib.urlretrieve(url, f.name)
        path = f.name
    return path


def _get_installed_file_contents(env, local_file):
    return open(_get_installed_file(env, local_file), "r").read()


def _write_to_file(contents, path, mode):
    """
    Use fabric to write string contents to remote file specified by path.
    """
    fd, local_path = tempfile.mkstemp()
    try:
        os.write(fd, contents)
        put(local_path, path, use_sudo=True, mode=mode)
        os.close(fd)
    finally:
        os.unlink(local_path)


def _get_bin_dir(env):
    """ 
    When env.system_install is /usr this exists, but in the Galaxy
    it may not already exist.
    """
    return _get_install_subdir(env, "bin")


def _get_include_dir(env):
    return _get_install_subdir(env, "include")


def _get_lib_dir(env):
    return _get_install_subdir(env, "lib")


def _get_install_subdir(env, subdir):
    path = os.path.join(env.system_install, subdir)
    if not exists(path):
        env.safe_sudo("mkdir -p '%s'" % path)
    return path


def _set_default_config(env, install_dir):
    """
    Sets up default galaxy config directory symbolic link (if needed). Needed
    when it doesn't exists or when installing a new version of software.
    """
    version = env["tool_version"]
    if exists(install_dir):
        install_dir_root = os.path.join(install_dir, "..")
        replace_default = False
        if not exists(os.path.join(install_dir_root, "default")):
            replace_default = True
        if not replace_default:
            default_version = env.safe_sudo("basename `readlink -f %s/default`" % install_dir_root)
            if version > default_version:  # Bug: Wouldn't work for 1.9 < 1.10
                print "default version %s is older than version %s just installed, replacing..." % (default_version, version)
                replace_default = True
        if replace_default:
            env.safe_sudo("rm -rf %s/default; ln -f -s %s %s/default" % (install_dir_root, install_dir, install_dir_root))


def _setup_simple_service(service_name):
    """
    Very Ubuntu/Debian specific, will need to be modified if used on other
    archs.
    """
    sudo("ln -f -s /etc/init.d/%s /etc/rc0.d/K01%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc1.d/K01%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc2.d/S99%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc3.d/S99%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc4.d/S99%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc5.d/S99%s" % (service_name, service_name))
    sudo("ln -f -s /etc/init.d/%s /etc/rc6.d/K01%s" % (service_name, service_name))


def _render_config_file_temlate(env, name, defaults={}, overrides={}, default_source=None):
    param_prefix = name.replace(".", "_")
    # Deployer can specify absolute path for config file, check this first
    path_key_name = "%s_path" % param_prefix
    template_key_name = "%s_template_path" % param_prefix
    if env.get(path_key_name, None):
        source_path = env[path_key_name]
        source_template = False
    elif env.get(template_key_name, None):
        source_path = env[template_key_name]
        source_template = True
    elif default_source:
        source_path = _get_installed_file(env, default_source)
        source_template = source_path.endswith(".template")
    else:
        default_template_name = "%s.template" % name
        source_path = _get_installed_file(env, default_template_name)
        source_template = True

    if source_template:
        template = Template(open(source_path, "r").read())
        template_params = _extend_env(env, defaults=defaults, overrides=overrides)
        contents = template.substitute(template_params)
    else:
        contents = open(source_path, "r").read()
    return contents


def _extend_env(env, defaults={}, overrides={}):
    new_env = {}
    for key, value in defaults.iteritems():
        new_env[key] = value
    for key, value in env.iteritems():
        new_env[key] = value
    for key, value in overrides.iteritems():
        new_env[key] = value
    return new_env


def _setup_conf_file(env, dest, name, defaults={}, overrides={}, default_source=None):
    conf_file_contents = _render_config_file_temlate(env, name, defaults, overrides, default_source)
    _write_to_file(conf_file_contents, dest, mode=0755)
