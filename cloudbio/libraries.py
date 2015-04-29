"""Installers for programming language specific libraries.
"""
import os

from fabric.api import env, cd, settings
from cloudbio import fabutils
from cloudbio.custom import shared

def r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.
    """
    with shared._make_tmp_dir() as tmp_dir:
        with cd(tmp_dir):
            # Create an Rscript file with install details.
            out_file = os.path.join(tmp_dir, "install_packages.R")
            _make_install_script(out_file, config)
            # run the script and then get rid of it
            rscript = fabutils.find_cmd(env, "Rscript", "--version")
            if rscript:
                env.safe_run("%s %s" % (rscript, out_file))
            else:
                env.logger.warn("Rscript not found; skipping install of R libraries.")
            env.safe_run("rm -f %s" % out_file)

def _make_install_script(out_file, config):
    if env.safe_exists(out_file):
        env.safe_run("rm -f %s" % out_file)
    env.safe_run("touch %s" % out_file)
    lib_loc = os.path.join(env.system_install, "lib", "R", "site-library")
    env.safe_sudo("mkdir -p %s" % lib_loc)
    with settings(warn_only=True):
        env.safe_sudo("chown -R %s %s" % (env.user, lib_loc))
    repo_info = """
    .libPaths(c("%s"))
    library(methods)
    cran.repos <- getOption("repos")
    cran.repos["CRAN" ] <- "%s"
    options(repos=cran.repos)
    source("%s")
    """ % (lib_loc, config["cranrepo"], config["biocrepo"])
    env.safe_append(out_file, repo_info)
    install_fn = """
    repo.installer <- function(repos, install.fn, pkg_name_fn) {
      %s
      maybe.install <- function(pname) {
        check_name <- ifelse(is.null(pkg_name_fn), pname, pkg_name_fn(pname))
        if (!(is.element(check_name, installed.packages()[,1])))
          install.fn(pname)
      }
    }
    """
    if config.get("update_packages", True):
        update_str = """
        update.packages(lib.loc="%s", repos=repos, ask=FALSE)
        """ % lib_loc
    else:
        update_str = "\n"
    env.safe_append(out_file, install_fn % update_str)
    std_install = """
    std.pkgs <- c(%s)
    std.installer = repo.installer(cran.repos, install.packages, NULL)
    lapply(std.pkgs, std.installer)
    """ % (", ".join('"%s"' % p for p in config['cran']))
    env.safe_append(out_file, std_install)
    if len(config.get("bioc", [])) > 0:
        bioc_install = """
        bioc.pkgs <- c(%s)
        bioc.installer = repo.installer(biocinstallRepos(), biocLite, NULL)
        lapply(bioc.pkgs, bioc.installer)
        """ % (", ".join('"%s"' % p for p in config['bioc']))
        env.safe_append(out_file, bioc_install)
    if config.get("cran-after-bioc"):
        std2_install = """
        std2.pkgs <- c(%s)
        lapply(std2.pkgs, std.installer)
        """ % (", ".join('"%s"' % p for p in config['cran-after-bioc']))
        env.safe_append(out_file, std2_install)
    if config.get("github"):
        dev_install = """
        library(devtools)
        github.pkgs <- c(%s)
        get_pkg_name <- function(orig) {
          unlist(strsplit(unlist(strsplit(orig, "/"))[2], "@"))[1]
        }
        github_installer = repo.installer(NULL, install_github, get_pkg_name)
        lapply(github.pkgs, github_installer)
        """ % (", ".join('"%s"' % p for p in config['github']))
        env.safe_append(out_file, dev_install)