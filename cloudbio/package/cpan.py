"""Install perl packages using CPAN and cpanminus (cpanm).
"""
import os

from fabric.api import cd, settings

from cloudbio.flavor.config import get_config_file
from cloudbio.fabutils import find_cmd
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.custom import shared as cshared

def install_packages(env):
    config_file = get_config_file(env, "perl-libs.yaml")
    (packages, _) = _yaml_to_packages(config_file.base, subs_yaml_file=config_file.dist, namesort=False)
    cpanm_cmd = find_cmd(env, "cpanm", "--version")
    for package in packages:
        if package.count("==") > 1:
            _install_from_url(env, cpanm_cmd, package)
        else:
            _install_from_cpan(env, cpanm_cmd, package)

def _install_from_cpan(env, cpanm_cmd, package):
    """Install from CPAN using cpanm, handling special arguments.

    The simplest input string is just a package to install (like XML::Simple) but
    users can also specify build arguments and exports as additional items separated
    by ';'
    """
    parts = package.split(";")
    if len(parts) == 1:
        perl_lib = parts[0]
        args = ""
        exports = []
    elif len(parts) == 2:
        perl_lib, args = parts
        exports = []
    else:
        perl_lib, args = parts[:2]
        exports = parts[2:]
    export_strs = []
    for export in exports:
        export_strs.append("export " + export.format(system_install=env.system_install))
    export = " && ".join(export_strs) + " && " if export_strs else ""
    build_args = ("--build-args='%s'" % args) if args else ""
    env.safe_run("%s %s -i --notest --local-lib=%s %s '%s'" % (export, cpanm_cmd, env.system_install,
                                                                build_args, perl_lib))

def _install_from_url(env, cpanm_cmd, package):
    """Check version of a dependency and download and install with cpanm if not up to date.

    Packages installed via URL have the package name, target version and URL separated
    with '=='. They can also optionally have a build directory or dependency to remove.
    """
    parts = package.split("==")
    package, target_version, url = parts[:3]
    args = {}
    if len(parts) > 3:
        for key, value in (x.split("=") for x in parts[3:]):
            args[key] = value
    with settings(warn_only=True):
        cur_version = env.safe_run_output("export PERL5LIB=%s/lib/perl5:${PERL5LIB} && " % env.system_install +
                                          """perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' %s"""
                                          % package)
    if cur_version != target_version:
        with cshared._make_tmp_dir() as work_dir:
            with cd(work_dir):
                dl_dir = cshared._fetch_and_unpack(url)
                if args.get("build"):
                    dl_dir = os.path.join(dl_dir, args["build"])
                with cd(dl_dir):
                    if args.get("depremove"):
                        for fname in ["Makefile.PL", "MYMETA.json", "MYMETA.yml"]:
                            env.safe_run(r"""sed -i.bak -e '/^.*%s.*/s/^/#/' %s""" % (args["depremove"], fname))
                    env.safe_run("%s -i --notest --local-lib=%s ." % (cpanm_cmd, env.system_install))
