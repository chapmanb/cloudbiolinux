"""Automated installation on debian package systems with apt.
"""
from fabric.api import *
from fabric.contrib.files import *

from cloudbio.package.shared import _yaml_to_packages

def _apt_packages(to_install=None, pkg_list=None, pkg_config_file_path=None):
    """ 
    Install packages available via apt-get.
    Note that ``to_install`` and ``pkg_list`` arguments cannot be used simultaneously.

    :type to_install:  list
    :param to_install: A list of strings (ie, groups) present in the ``main.yaml``
                       config file that will be used to filter out the specific
                       packages to be installed.

    :type pkg_list:  list
    :param pkg_list: An explicit list of packages to install. No other files,
                     flavors, or editions are considered.

    :type pkg_config_file_path:  string
    :param pkg_config_file_path: Allows a custom path to be specified where
                                 ``packages.yaml`` and ``packages-[dist].yaml``
                                 are stored. Note that the file names cannot
                                 be customized, only the path.
    """
    env.logger.info("Update the system")
    sudo("apt-get update") # Always update
    if to_install is not None:
        env.logger.info("Will install all packages (as listed in packages*.yaml "
                "files in {0})".format(pkg_config_file_path))
        if pkg_config_file_path is None:
            pkg_config_file_path = env.config_dir
        pkg_config_file = os.path.join(pkg_config_file_path, "packages.yaml")
        subs_pkg_config_file = os.path.join(pkg_config_file_path, "packages-%s.yaml" %
                                       env.distribution)
        if not os.path.exists(subs_pkg_config_file): subs_pkg_config_file = None
        env.edition.apt_upgrade_system()
        # Retrieve final package names
        (packages, _) = _yaml_to_packages(pkg_config_file, to_install,
                                          subs_pkg_config_file)
        # At this point allow the Edition to rewrite the package list - 
        # this is shared within and between editions.
        # Ref:  https://github.com/chapmanb/cloudbiolinux/pull/10#issuecomment-1616423
        packages = env.edition.rewrite_config_items("packages", packages)

        # At this point allow the Flavor to rewrite the package list
        packages = env.flavor.rewrite_config_items("packages", packages)
    elif pkg_list is not None:
        env.logger.info("Will install specific packages: {0}".format(pkg_list))
        packages = pkg_list
    else:
        raise ValueError("Need a file with packages or a list of packages")
    # A single line install is much faster - note that there is a max
    # for the command line size, so we do 30 at a time
    group_size = 30
    i = 0
    env.logger.info("Installing %i packages" % len(packages))
    while i < len(packages):
        env.logger.info("Package install progress: {0}/{1}".format(i, len(packages)))
        sudo("apt-get -y --force-yes install %s" % " ".join(packages[i:i+group_size]))
        i += group_size
    sudo("apt-get clean")

def _add_apt_gpg_keys():
    """Adds GPG keys from all repositories
    """
    env.logger.info("Update GPG keys for repositories")
    standalone = [
        "http://archive.cloudera.com/debian/archive.key",
        'http://download.virtualbox.org/virtualbox/debian/oracle_vbox.asc'
    ]
    keyserver = [
            ("keyserver.ubuntu.com", "7F0CEB10"),
            ("keyserver.ubuntu.com", "E084DAB9"),
            ("keyserver.ubuntu.com", "D67FC6EAE2A11821"),
        ]
    standalone, keyserver = env.edition.rewrite_apt_keys(standalone, keyserver)
    for key in standalone:
        sudo("wget -q -O- %s | apt-key add -" % key)
    for url, key in keyserver:
        sudo("apt-key adv --keyserver %s --recv %s" % (url, key))

def _setup_apt_automation():
    """Setup the environment to be fully automated for tricky installs.

    Sun Java license acceptance:
    http://www.davidpashley.com/blog/debian/java-license

    MySQL root password questions; install with empty root password:
    http://snowulf.com/archives/540-Truly-non-interactive-unattended-apt-get-install.html

    Postfix, setup for no configuration. See more on issues here:
    http://www.uluga.ubuntuforums.org/showthread.php?p=9120196
    """
    interactive_cmd = "export DEBIAN_FRONTEND=noninteractive"
    if not contains(env.shell_config, interactive_cmd):
        append(env.shell_config, interactive_cmd)
    package_info = [
            "postfix postfix/main_mailer_type select No configuration",
            "postfix postfix/mailname string notusedexample.org",
            "mysql-server-5.1 mysql-server/root_password string '(password omitted)'",
            "mysql-server-5.1 mysql-server/root_password_again string '(password omitted)'",
            "sun-java6-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-bin shared/accepted-sun-dlj-v1-1 select true",
            "grub-pc grub2/linux_cmdline string ''",
            "grub-pc grub-pc/install_devices_empty boolean true",
            "acroread acroread/default-viewer boolean false",
            "rabbitmq-server rabbitmq-server/upgrade_previous note",
            ]
    package_info = env.edition.rewrite_apt_automation(package_info)
    cmd = ""
    for l in package_info:
        cmd += "echo %s | /usr/bin/debconf-set-selections ; " % l
    sudo(cmd)

def _setup_apt_sources():
    """Add sources for retrieving library packages.
       Using add-apt-repository allows processing PPAs (on Ubuntu)

       This method modifies the apt sources file.

       Uses python-software-properties, which provides an abstraction of apt repositories
    """

    # It may be sudo is not installed - which has fab fail - therefor
    # we'll try to install it by default, assuming we have root access
    # already (e.g. on EC2). Fab will fail anyway, otherwise.
    if not exists('/usr/bin/sudo') or not exists('/usr/bin/curl'):
        sudo('apt-get update')
        sudo('apt-get -y --force-yes install sudo curl')

    env.logger.debug("_setup_apt_sources " + env.sources_file + " " + env.edition.name)
    env.edition.check_packages_source()
    comment = "# This file was modified for "+ env.edition.name
    # Setup apt download policy (default is None)
    # (see also https://help.ubuntu.com/community/PinningHowto)
    preferences = env.edition.rewrite_apt_preferences([])
    if len(preferences):
        # make sure it exists, and is empty
        sudo("rm -f %s" % env.apt_preferences_file)
        sudo("touch %s" % env.apt_preferences_file)
        append(env.apt_preferences_file, comment, use_sudo=True)
        lines = "\n".join(preferences)
        env.logger.debug("Policy %s" % lines)
        # append won't duplicate, so we use echo
        sudo("/bin/echo -e \"%s\" >> %s" % (lines, env.apt_preferences_file))
        # check there is no error parsing the file
        env.logger.debug(sudo("apt-cache policy"))

    # Make sure a source file exists
    if not exists(env.sources_file):
        sudo("touch %s" % env.sources_file)
    # Add a comment
    if not contains(env.sources_file, comment):
        append(env.sources_file, comment, use_sudo=True)
    for source in env.edition.rewrite_apt_sources_list(env.std_sources):
        env.logger.debug("Source %s" % source)
        if source.startswith("ppa:"):
            sudo("apt-get install -y --force-yes python-software-properties")
            sudo("add-apt-repository '%s'" % source)
        elif not contains(env.sources_file, source): # FIXME: append never adds dups!
            append(env.sources_file, source, use_sudo=True)
