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
from cloudbio.fabutils import quiet

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
        result = env.safe_run("export PATH=$PATH:%s/bin && "
                              "export LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH && %s" %
                              (env.system_install, env.system_install, pname))
    return result.return_code == 127


def _galaxy_tool_install(args):
    try:
        return args[0]["galaxy_tool_install"]
    except:
        return False


def _galaxy_tool_present(args):
    return env.safe_exists(os.path.join(args[0]["system_install"], "env.sh"))


def _if_not_python_lib(library):
    """Decorator that checks if a python library is installed.
    """
    def argcatcher(func):
        functools.wraps(func)

        def decorator(*args, **kwargs):
            with settings(warn_only=True):
                result = env.safe_run("%s -c 'import %s'" % (_python_cmd(env), library))
            if result.failed:
                return func(*args, **kwargs)
        return decorator
    return argcatcher


@contextmanager
def _make_tmp_dir():
    """
    Setup a temporary working directory for building custom software. First checks
    fabric environment for a `work_dir` path, if that is not set it will use the
    remote path $TMPDIR/cloudbiolinux if $TMPDIR is defined remotely, finally falling
    back on remote $HOME/cloudbiolinux otherwise.
    """
    work_dir = __work_dir()
    use_sudo = False
    if not env.safe_exists(work_dir):
        with settings(warn_only=True):
            # Try to create this directory without using sudo, but
            # if needed fallback.
            result = env.safe_run("mkdir -p '%s'" % work_dir)
            if result.return_code != 0:
                use_sudo = True
        if use_sudo:
            env.safe_sudo("mkdir -p '%s'" % work_dir)
            env.safe_sudo("chown -R %s '%s'" % (env.user, work_dir))
    yield work_dir
    if env.safe_exists(work_dir):
        run_func = env.safe_sudo if use_sudo else env.safe_run
        run_func("rm -rf %s" % work_dir)


def __work_dir():
    work_dir = env.get("work_dir", None)
    if not work_dir:
        with quiet():
            tmp_dir = env.safe_run_output("echo $TMPDIR")
        if tmp_dir.failed or not tmp_dir.strip():
            home_dir = env.safe_run_output("echo $HOME")
            tmp_dir = os.path.join(home_dir, "tmp")
        work_dir = os.path.join(tmp_dir.strip(), "cloudbiolinux")
    return work_dir


# -- Standard build utility simplifiers


def _get_expected_file(url, dir_name=None, safe_tar=False, tar_file_name=None):
    if tar_file_name:
        tar_file = tar_file_name
    else:
        tar_file = os.path.split(url.split("?")[0])[-1]
    safe_tar = "--pax-option='delete=SCHILY.*,delete=LIBARCHIVE.*'" if safe_tar else ""
    exts = {(".tar.gz", ".tgz"): "tar %s -xzpf" % safe_tar,
            (".tar",): "tar %s -xpf" % safe_tar,
            (".tar.bz2",): "tar %s -xjpf" % safe_tar,
            (".zip",): "unzip"}
    for ext_choices, tar_cmd in exts.iteritems():
        for ext in ext_choices:
            if tar_file.endswith(ext):
                if dir_name is None:
                    dir_name = tar_file[:-len(ext)]
                return tar_file, dir_name, tar_cmd
    raise ValueError("Did not find extract command for %s" % url)


def _safe_dir_name(dir_name, need_dir=True):
    replace_try = ["", "-src", "_core"]
    for replace in replace_try:
        check = dir_name.replace(replace, "")
        if env.safe_exists(check):
            return check
    # still couldn't find it, it's a nasty one
    for check_part in (dir_name.split("-")[0].split("_")[0],
                       dir_name.split("-")[-1].split("_")[-1],
                       dir_name.split(".")[0],
                       dir_name.lower().split(".")[0]):
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            dirs = env.safe_run_output("ls -d1 *%s*/" % check_part).split("\n")
            dirs = [x for x in dirs if "cannot access" not in x and "No such" not in x]
        if len(dirs) == 1 and dirs[0]:
            return dirs[0]
    if need_dir:
        raise ValueError("Could not find directory %s" % dir_name)


def _fetch_and_unpack(url, need_dir=True, dir_name=None, revision=None,
                      safe_tar=False, tar_file_name=None):
    if url.startswith(("git", "svn", "hg", "cvs")):
        base = os.path.splitext(os.path.basename(url.split()[-1]))[0]
        if env.safe_exists(base):
            env.safe_sudo("rm -rf {0}".format(base))
        env.safe_run(url)
        if revision:
            if url.startswith("git"):
                env.safe_run("cd %s && git checkout %s" % (base, revision))
            else:
                raise ValueError("Need to implement revision retrieval for %s" % url.split()[0])
        return base
    else:
        # If tar_file_name is provided, use it instead of the inferred one
        tar_file, dir_name, tar_cmd = _get_expected_file(url, dir_name, safe_tar, tar_file_name=tar_file_name)
        if not env.safe_exists(tar_file):
            env.safe_run("wget --no-check-certificate -O %s '%s'" % (tar_file, url))
        env.safe_run("%s %s" % (tar_cmd, tar_file))
        return _safe_dir_name(dir_name, need_dir)


def _configure_make(env):
    env.safe_run("export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:%s/lib/pkgconfig && " \
                 "./configure --disable-werror --prefix=%s " %
                 (env.system_install, env.system_install))
    env.safe_run("make")
    env.safe_sudo("make install")


def _make_copy(find_cmd=None, premake_cmd=None, do_make=True):
    def _do_work(env):
        if premake_cmd:
            premake_cmd()
        if do_make:
            env.safe_run("make")
        if find_cmd:
            install_dir = _get_bin_dir(env)
            for fname in env.safe_run_output(find_cmd).split("\n"):
                env.safe_sudo("cp -rf %s %s" % (fname.rstrip("\r"), install_dir))
    return _do_work


def _get_install(url, env, make_command, post_unpack_fn=None, revision=None, dir_name=None,
                 safe_tar=False, tar_file_name=None):
    """Retrieve source from a URL and install in our system directory.
    """
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            dir_name = _fetch_and_unpack(url, revision=revision, dir_name=dir_name,
                                         safe_tar=safe_tar, tar_file_name=tar_file_name)
        with cd(os.path.join(work_dir, dir_name)):
            if post_unpack_fn:
                post_unpack_fn(env)
            make_command(env)


def _get_install_local(url, env, make_command, dir_name=None,
                       post_unpack_fn=None, safe_tar=False, tar_file_name=None):
    """Build and install in a local directory.
    """
    (_, test_name, _) = _get_expected_file(url, safe_tar=safe_tar, tar_file_name=tar_file_name)
    test1 = os.path.join(env.local_install, test_name)
    if dir_name is not None:
        test2 = os.path.join(env.local_install, dir_name)
    elif "-" in test1:
        test2, _ = test1.rsplit("-", 1)
    else:
        test2 = os.path.join(env.local_install, test_name.split("_")[0])
    if not env.safe_exists(test1) and not env.safe_exists(test2):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                dir_name = _fetch_and_unpack(url, dir_name=dir_name, safe_tar=safe_tar,
                    tar_file_name=tar_file_name)
                if not env.safe_exists(os.path.join(env.local_install, dir_name)):
                    with cd(dir_name):
                        if post_unpack_fn:
                            post_unpack_fn(env)
                        make_command(env)
                    # Copy instead of move because GNU mv does not have --parents flag.
                    # The source dir will get cleaned up anyhow so just leave it.
                    destination_dir = env.local_install
                    env.safe_sudo("mkdir -p '%s'" % destination_dir)
                    env.safe_sudo("cp --recursive %s %s" % (dir_name, destination_dir))

# --- Language specific utilities


def _symlinked_shared_dir(pname, version, env, extra_dir=None):
    """Create a symlinked directory of files inside the shared environment.
    """
    if extra_dir:
        base_dir = os.path.join(env.system_install, "share", extra_dir, pname)
    else:
        base_dir = os.path.join(env.system_install, "share", pname)
    install_dir = "%s-%s" % (base_dir, version)
    if not env.safe_exists(install_dir):
        env.safe_sudo("mkdir -p %s" % install_dir)
        if env.safe_exists(base_dir):
            env.safe_sudo("rm -f %s" % base_dir)
        env.safe_sudo("ln -s %s %s" % (install_dir, base_dir))
        return install_dir
    return None


def _symlinked_java_version_dir(pname, version, env):
    return _symlinked_shared_dir(pname, version, env, extra_dir="java")


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


def _python_cmd(env):
    """Retrieve python command, handling tricky situations on CentOS.
    """
    if "python_version_ext" in env and env.python_version_ext:
        major, minor = env.safe_run("python --version").split()[-1].split(".")[:2]
        check_major, check_minor = env.python_version_ext.split(".")[:2]
        if major != check_major or int(check_minor) > int(minor):
            return "python%s" % env.python_version_ext
        else:
            return "python"
    else:
        return "python"


def _pip_cmd(env):
    """Retrieve pip command for installing python packages, allowing configuration.
    """
    to_check = ["pip"]
    if "pip_cmd" in env and env.pip_cmd:
        to_check.append(env.pip_cmd)
    if not env.use_sudo:
        to_check.append(os.path.join(env.system_install, "bin", "pip"))
    if "python_version_ext" in env and env.python_version_ext:
        to_check.append("pip-{0}".format(env.python_version_ext))
    for cmd in to_check:
        with quiet():
            pip_version = env.safe_run("%s --version" % cmd)
        if pip_version.succeeded:
            return cmd
    raise ValueError("Could not find pip installer from: %s" % to_check)


def _python_make(env):
    env.safe_sudo("%s install --upgrade ." % _pip_cmd(env))
    for clean in ["dist", "build", "lib/*.egg-info"]:
        env.safe_sudo("rm -rf %s" % clean)


def _get_installed_file(env, local_file):
    installed_files_dir = \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "installed_files")
    path = os.path.join(installed_files_dir, local_file)
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
    if not env.safe_exists(path):
        env.safe_sudo("mkdir -p '%s'" % path)
    return path


def _set_default_config(env, install_dir, sym_dir_name="default"):
    """
    Sets up default galaxy config directory symbolic link (if needed). Needed
    when it doesn't exists or when installing a new version of software.
    """
    version = env["tool_version"]
    if env.safe_exists(install_dir):
        install_dir_root = "%s/.." % install_dir
        sym_dir = "%s/%s" % (install_dir_root, sym_dir_name)
        replace_default = False
        if not env.safe_exists(sym_dir):
            replace_default = True
        if not replace_default:
            default_version = env.safe_sudo("basename `readlink -f %s`" % sym_dir)
            if version > default_version:  # Bug: Wouldn't work for 1.9 < 1.10
                print "default version %s is older than version %s just installed, replacing..." % (default_version, version)
                replace_default = True
        if replace_default:
            env.safe_sudo("rm -rf %s; ln -f -s %s %s" % (sym_dir, install_dir, sym_dir))


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


def _render_config_file_template(env, name, defaults={}, overrides={}, default_source=None):
    """
    If ``name` is say ``nginx.conf``, check fabric environment for
    ``nginx_conf_path`` and then ``nginx_conf_template_path``. If
    ``nginx_conf_path`` is set, return the contents of that file. If
    nginx_conf_template_path is set, return the contents of that file
    but with variable interpolation performed. Variable interpolation
    is performed using a derivative of the fabric environment defined
    using the supplied ``defaults`` and ``overrides`` using the
    ``_extend_env`` function below.

    Finally, if neither ``nginx_conf_path`` or
    ``nginx_conf_template_path`` are set, check the
    ``installed_files`` directory for ``nginx.conf`` and finally
    ``nginx.conf.template``.
    """
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
    """
    Create a new ``dict`` from fabric's ``env``, first adding defaults
    specified via ``defaults`` (if available). Finally, override
    anything in env, with values specified by ``overrides``.
    """
    new_env = {}
    for key, value in defaults.iteritems():
        new_env[key] = value
    for key, value in env.iteritems():
        new_env[key] = value
    for key, value in overrides.iteritems():
        new_env[key] = value
    return new_env


def _setup_conf_file(env, dest, name, defaults={}, overrides={}, default_source=None):
    conf_file_contents = _render_config_file_template(env, name, defaults, overrides, default_source)
    _write_to_file(conf_file_contents, dest, mode=0755)


def _add_to_profiles(line, profiles=[], use_sudo=True):
    """
    If it's not already there, append ``line`` to shell profiles files.
    By default, these are ``/etc/profile`` and ``/etc/bash.bashrc`` but can be
    overridden by providing a list of file paths to the ``profiles`` argument.
    """
    if not profiles:
        profiles = ['/etc/bash.bashrc', '/etc/profile']
    for profile in profiles:
        if not env.safe_contains(profile, line):
            env.safe_append(profile, line, use_sudo=use_sudo)


def install_venvburrito():
    """
    If not already installed, install virtualenv-burrito
    (https://github.com/brainsik/virtualenv-burrito) as a convenient
    method for installing and managing Python virtualenvs.
    """
    url = "https://raw.github.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh"
    if not env.safe_exists("$HOME/.venvburrito/startup.sh"):
        env.safe_run("curl -s {0} | $SHELL".format(url))
        # Add the startup script into the ubuntu user's bashrc
        _add_to_profiles(". $HOME/.venvburrito/startup.sh", ['/home/ubuntu/.bashrc'], use_sudo=False)


def _create_python_virtualenv(env, venv_name, reqs_file=None, reqs_url=None):
    """
    Using virtual-burrito, create a new Python virtualenv named ``venv_name``.
    Do so only if the virtualenv of the given name does not already exist.
    virtual-burrito installs virtualenvs in ``$HOME/.virtualenvs``.

    By default, an empty virtualenv is created. Python libraries can be
    installed into the virutalenv at the time of creation by providing a path
    to the requirements.txt file (``reqs_file``). Instead of providing the file,
    a url to the file can be provided via ``reqs_url``, in which case the
    requirements file will first be downloaded. Note that if the ``reqs_url``
    is provided, the downloaded file will take precedence over ``reqs_file``.
    """
    # First make sure virtualenv-burrito is installed
    install_venvburrito()
    activate_vburrito = ". $HOME/.venvburrito/startup.sh"

    def create():
        if "venv_directory" not in env:
            _create_global_python_virtualenv(env, venv_name, reqs_file, reqs_url)
        else:
            _create_local_python_virtualenv(env, venv_name, reqs_file, reqs_url)

    # TODO: Terrible hack here, figure it out and fix it.
    #   prefix or vburrito does not work with is_local or at least deployer+is_local
    if env.is_local:
        create()
    else:
        with prefix(activate_vburrito):
            create()


def _create_local_python_virtualenv(env, venv_name, reqs_file, reqs_url):
    """
    Use virtualenv directly to setup virtualenv in specified directory.
    """
    venv_directory = env.get("venv_directory")
    if not env.safe_exists(venv_directory):
        if reqs_url:
                env.safe_sudo("wget --output-document=%s %s" % (reqs_file, reqs_url))
        env.logger.debug("Creating virtualenv in directory %s" % venv_directory)
        env.safe_sudo("virtualenv --no-site-packages '%s'" % venv_directory)
        env.logger.debug("Activating")
        env.safe_sudo(". %s/bin/activate; pip install -r '%s'" % (venv_directory, reqs_file))


def _create_global_python_virtualenv(env, venv_name, reqs_file, reqs_url):
    """
    Use mkvirtualenv to setup this virtualenv globally for user.
    """
    if venv_name in env.safe_run("lsvirtualenv | grep {0} || true"
        .format(venv_name)):
        env.logger.info("Virtualenv {0} already exists".format(venv_name))
    else:
        with _make_tmp_dir():
            if reqs_file or reqs_url:
                if not reqs_file:
                    # This mean the url only is provided so 'standardize ' the file name
                    reqs_file = 'requirements.txt'
                cmd = "mkvirtualenv -r {0} {1}".format(reqs_file, venv_name)
            else:
                cmd = "mkvirtualenv {0}".format(venv_name)
            if reqs_url:
                env.safe_run("wget --output-document=%s %s" % (reqs_file, reqs_url))
            env.safe_run(cmd)
            env.logger.info("Finished installing virtualenv {0}".format(venv_name))


def _read_boolean(env, name, default):
    property_str = env.get(name, str(default))
    return property_str.upper() in ["TRUE", "YES"]
