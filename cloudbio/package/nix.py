"""Install software with the Nix package manager.
"""
from fabric.api import *
from fabric.contrib.files import *

from cloudbio.package.shared import _yaml_to_packages
from cloudbio.flavor.config import get_config_file

def _setup_nix_sources():
    if env.nixpkgs:
        target_info = run("uname -a")
        env.logger.info("Target: "+target_info)
        # find the target architecture, if not preset
        if not env.has_key("arch"):
          env.arch = run("uname -m")

     # first override the path
        append("/root/.bashrc", "export PATH=$HOME/.nix-profile/bin:$PATH", use_sudo=True)
        env.logger.info("Checking NixPkgs")
        if not exists("/nix/store"):
            # first time installation
            if not exists("/usr/bin/nix-env"):
               # install Nix (standard Debian release)
               nix_deb = "nix_0.16-1_"+env.arch+".deb"
               if not exists(nix_deb):
                   # run("wget http://hydra.nixos.org/build/565031/download/1/nix_0.16-1_i386.deb")
                   run("wget http://hydra.nixos.org/build/565048/download/1/"+nix_deb)
                   sudo("dpkg -i "+nix_deb)
        run("nix-channel --list")
        if run("nix-channel --list") == "":
            # Setup channel
            sudo("nix-channel --add http://nixos.org/releases/nixpkgs/channels/nixpkgs-unstable")
        sudo("nix-channel --update")
        # upgrade Nix to latest (and remove the older version, as it is much slower)
        sudo("nix-env -b -i nix")
        if exists("/usr/bin/nix-env"):
            env.logger.info("uninstall older Nix (Debian release)")
            sudo("dpkg -r nix")

def _nix_packages(to_install):
    """Install packages available via nixpkgs (optional)
    """
    if env.nixpkgs:
        env.logger.info("Update and install NixPkgs packages")
        pkg_config_file = get_config_file(env, "packages-nix.yaml").base
        sudo("nix-channel --update")
        # Retrieve final package names
        (packages, _) = _yaml_to_packages(pkg_config_file, to_install)
        packages = env.flavor.rewrite_config_items("packages", packages)
        for p in packages:
            sudo("nix-env -b -i %s" % p)
