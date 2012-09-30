"""Install software and configure package managers.
"""

from cloudbio.package.deb import (_apt_packages, _add_apt_gpg_keys,
                                  _setup_apt_automation, _setup_apt_sources)
from cloudbio.package.rpm import (_yum_packages, _setup_yum_bashrc,
                                  _setup_yum_sources)


def _configure_and_install_native_packages(env, pkg_install):
    """
    Setups up native package repositories, determines list
    of native packages to install, and installs them.
    """
    if env.distribution in ["debian", "ubuntu"]:
        _setup_apt_sources()
        _setup_apt_automation()
        _add_apt_gpg_keys()
        _apt_packages(pkg_install)
    elif env.distribution in ["centos", "scientificlinux"]:
        _setup_yum_sources()
        _yum_packages(pkg_install)
        _setup_yum_bashrc()
    else:
        raise NotImplementedError("Unknown target distribution")
