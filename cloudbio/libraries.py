"""Installers for programming language specific libraries.
"""
import os

from fabric.api import env


def r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.
    """
    # Create an Rscript file with install details.
    out_file = "install_packages.R"
    if env.safe_exists(out_file):
        env.safe_run("rm -f %s" % out_file)
    env.safe_run("touch %s" % out_file)
    lib_loc = os.path.join(env.system_install, "lib", "R", "site-library")
    env.safe_sudo("mkdir -p %s" % lib_loc)
    repo_info = """
    .libPaths(c("%s"))
    cran.repos <- getOption("repos")
    cran.repos["CRAN" ] <- "%s"
    options(repos=cran.repos)
    source("%s")
    """ % (lib_loc, config["cranrepo"], config["biocrepo"])
    env.safe_append(out_file, repo_info)
    install_fn = """
    repo.installer <- function(repos, install.fn) {
      %s
      maybe.install <- function(pname) {
        if (!(pname %%in%% installed.packages()))
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
    std.installer = repo.installer(cran.repos, install.packages)
    lapply(std.pkgs, std.installer)
    """ % (", ".join('"%s"' % p for p in config['cran']))
    env.safe_append(out_file, std_install)
    if len(config.get("bioc", [])) > 0:
        bioc_install = """
        bioc.pkgs <- c(%s)
        bioc.installer = repo.installer(biocinstallRepos(), biocLite)
        lapply(bioc.pkgs, bioc.installer)
        """ % (", ".join('"%s"' % p for p in config['bioc']))
        env.safe_append(out_file, bioc_install)
    # run the script and then get rid of it
    env.safe_sudo("Rscript %s" % out_file)
    env.safe_run("rm -f %s" % out_file)
