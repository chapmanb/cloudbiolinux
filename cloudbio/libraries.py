"""Installers for programming language specific libraries.
"""

from fabric.api import *
from fabric.contrib.files import *

def r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.
    """
    # Create an Rscript file with install details.
    out_file = "install_packages.R"
    if exists(out_file):
        run("rm -f %s" % out_file)
    run("touch %s" % out_file)
    repo_info = """
    cran.repos <- getOption("repos")
    cran.repos["CRAN" ] <- "%s"
    options(repos=cran.repos)
    source("%s")
    """ % (config["cranrepo"], config["biocrepo"])
    append(out_file, repo_info)
    install_fn = """
    repo.installer <- function(repos, install.fn) {
      update.or.install <- function(pname) {
        if (pname %in% installed.packages())
          update.packages(lib.loc=c(pname), repos=repos, ask=FALSE)
        else
          install.fn(pname)
      }
    }
    """
    append(out_file, install_fn)
    std_install = """
    std.pkgs <- c(%s)
    std.installer = repo.installer(cran.repos, install.packages)
    lapply(std.pkgs, std.installer)
    """ % (", ".join('"%s"' % p for p in config['cran']))
    append(out_file, std_install)
    if len(config.get("bioc", [])) > 0:
        bioc_install = """
        bioc.pkgs <- c(%s)
        bioc.installer = repo.installer(biocinstallRepos(), biocLite)
        lapply(bioc.pkgs, bioc.installer)
        """ % (", ".join('"%s"' % p for p in config['bioc']))
        append(out_file, bioc_install)
    if config.get("update_packages", True):
        final_update = """
        update.packages(repos=biocinstallRepos(), ask=FALSE)
        update.packages(ask=FALSE)
        """
        append(out_file, final_update)
    # run the script and then get rid of it
    env.safe_sudo("Rscript %s" % out_file)
    run("rm -f %s" % out_file)
