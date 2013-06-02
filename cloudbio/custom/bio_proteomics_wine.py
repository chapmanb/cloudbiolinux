from fabric.api import put, cd
from fabric.contrib.files import exists

from shared import (_make_tmp_dir, _fetch_and_unpack, _write_to_file, _get_bin_dir)

import os


def install_proteomics_wine_env(env):
    script_src = env.get("setup_proteomics_wine_env_script")
    script_dest = "%s/bin/setup_proteomics_wine_env.sh" % env.get("system_install")
    if not exists(script_dest):
        put(script_src, script_dest, mode=0755, use_sudo=True)


def install_multiplierz(env):
    """
    Assumes your wine environment contains an install Python 2.6
    in C:\Python26.
    """
    wine_user = _get_wine_user(env)

    install_proteomics_wine_env(env)
    env.safe_sudo("setup_proteomics_wine_env.sh", user=wine_user)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _fetch_and_unpack("hg clone http://multiplierz.hg.sourceforge.net:8000/hgroot/multiplierz/multiplierz")
            with cd("multiplierz"):
                wine_prefix = _get_wine_prefix(env)
                env.safe_sudo("%s; wine %s/drive_c/Python26/python.exe setup.py install" % (_conf_wine(env), wine_prefix), user=wine_user)


def install_proteowizard(env):
    build_id = "83083"
    version = "3_0_4472"
    url = "http://teamcity.labkey.org:8080/repository/download/bt36/%s:id/pwiz-bin-windows-x86-vc100-release-%s.tar.bz2?guest=1" % (build_id, version)
    install_dir = env.get("install_dir")
    share_dir = "%s/share/proteowizard" % install_dir
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _fetch_and_unpack(url, need_dir=False)
            env.safe_sudo("cp -r . '%s'" % share_dir)
    proteowizard_apps = ["msconvert", "msaccess", "chainsaw", "msdiff", "mspicture", "mscat", "txt2mzml", "MSConvertGUI", "Skyline", "Topograph", "SeeMS"]
    for app in proteowizard_apps:
        setup_wine_wrapper(env, "%s/%s" % (share_dir, app))


def install_morpheus(env):
    url = "http://www.chem.wisc.edu/~coon/Downloads/Morpheus/latest/Morpheus.zip"  # TODO:
    install_dir = env.get("install_dir")
    share_dir = "%s/share/morpheus" % install_dir
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _fetch_and_unpack(url, need_dir=False)
            env.safe_sudo("cp -r Morpheus/* '%s'" % share_dir)
    proteowizard_apps = ["morpheus_cl.exe", "Morpheus.exe"]
    for app in proteowizard_apps:
        setup_wine_wrapper(env, "%s/%s" % (share_dir, app))


def setup_wine_wrapper(env, to):
    basename = os.path.basename(to)
    contents = """#!/bin/bash
setup_proteomics_wine_env.sh
export WINEPREFIX=$HOME/.wine-proteomics
wine %s "$@"
""" % to
    bin_dir = _get_bin_dir(env)
    dest = "%s/%s" % (bin_dir, basename)
    _write_to_file(contents, dest, 0755)


def _conf_wine(env):
    return "export WINEPREFIX=%s" % _get_wine_prefix(env)


def _get_wine_prefix(env):
    wine_user = _get_wine_user(env)
    return "~%s/.wine-proteomics" % wine_user


def _get_wine_user(env):
    return env.get("wine_user", env.get("user"))
