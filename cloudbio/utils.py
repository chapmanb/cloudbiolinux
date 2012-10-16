"""Utilities for logging and progress tracking.
"""
import logging
import os
import sys

from fabric.main import load_settings
from fabric.colors import yellow, red, green, magenta
from fabric.api import settings, hide, cd, run
from fabric.contrib.files import exists

from cloudbio.edition import _setup_edition
from cloudbio.distribution import _setup_distribution_environment
from cloudbio.flavor import Flavor
from cloudbio.flavor.config import get_config_file


class ColorFormatter(logging.Formatter):
    """ Format log message based on the message level
        http://stackoverflow.com/questions/1343227/can-pythons-logging-format-be-modified-depending-on-the-message-log-level
    """
    # Setup formatters for each of the levels
    err_fmt  = red("ERR [%(filename)s(%(lineno)d)] %(msg)s")
    warn_fmt  = magenta("WARN [%(filename)s(%(lineno)d)]: %(msg)s")
    dbg_fmt  = yellow("DBG [%(filename)s]: %(msg)s")
    info_fmt = green("INFO: %(msg)s")

    def __init__(self, fmt="%(name)s %(levelname)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt
        # Replace the original format with one customized by logging level
        if record.levelno == 10:   # DEBUG
            self._fmt = ColorFormatter.dbg_fmt
        elif record.levelno == 20: # INFO
            self._fmt = ColorFormatter.info_fmt
        elif record.levelno == 30: # WARN
            self._fmt = ColorFormatter.warn_fmt
        elif record.levelno == 40: # ERROR
            self._fmt = ColorFormatter.err_fmt
        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)
        # Restore the original format configured by the user
        self._fmt = format_orig
        return result 

def _setup_logging(env):
    env.logger = logging.getLogger("cloudbiolinux")
    env.logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # Use custom formatter
    ch.setFormatter(ColorFormatter())
    env.logger.addHandler(ch)

def _update_biolinux_log(env, target, flavor):
    """Updates the VM so it contains information on the latest BioLinux
       update in /var/log/biolinux.log.

       The latest information is appended to the file and can be used to see if
       an installation/update has completed (see also ./test/test_vagrant).
    """
    if not target:
        target = env.get("target", None)
        if not target:
            target = "unknown"
        else:
            target = target.name
    if not flavor:
        flavor = env.get("flavor", None)
        if not flavor:
            flavor = "unknown"
        else:
            flavor = flavor.name
    logfn = "/var/log/biolinux.log"
    info = "Target="+target+"; Edition="+env.edition.name+"; Flavor="+flavor
    env.logger.info(info)
    env.safe_sudo("date +\"%D %T - Updated "+info+"\" >> "+logfn)


def _configure_fabric_environment(env, flavor=None, fabricrc_loader=None):
    if not fabricrc_loader:
        fabricrc_loader = _parse_fabricrc

    _setup_flavor(env, flavor)
    fabricrc_loader(env)
    _setup_edition(env)
    _setup_distribution_environment()  # get parameters for distro, packages etc.
    _create_local_paths(env)


def _setup_flavor(env, flavor):
    """Setup a flavor, providing customization hooks to modify CloudBioLinux installs.

    Specify flavor as a name, in which case we look it up in the standard flavor
    directory (contrib/flavor/your_flavor), or as a path to a flavor directory outside
    of cloudbiolinux.
    """
    env.flavor = Flavor(env)
    env.flavor_dir = None
    if flavor:
        # setup the directory for flavor customizations
        if os.path.isabs(flavor):
            flavor_dir = flavor
        else:
            flavor_dir = os.path.join(os.path.dirname(__file__), "..", "contrib", "flavor", flavor)
        assert os.path.exists(flavor_dir), \
            "Did not find directory {0} for flavor {1}".format(flavor_dir, flavor)
        env.flavor_dir = flavor_dir
        # Load python customizations to base configuration if present
        py_flavor = "{0}flavor".format(os.path.split(os.path.realpath(flavor_dir)))
        flavor_custom_py = os.path.join(flavor_dir, "{0}.py".format(py_flavor))
        if os.path.exists(flavor_custom_py):
            sys.path.append(flavor_dir)
            mod = __import__(flavor_dir, fromlist=[py_flavor])
    env.logger.info("This is a %s" % env.flavor.name)


def _parse_fabricrc(env):
    """Defaults from fabricrc.txt file; loaded if not specified at commandline.
    """
    env.config_dir = os.path.join(os.path.dirname(__file__), "..", "config")
    if not env.has_key("distribution") and not env.has_key("system_install"):
        env.logger.info("Reading default fabricrc.txt")
        env.update(load_settings(get_config_file(env, "fabricrc.txt").base))


def _create_local_paths(env):
    """Expand any paths defined in terms of shell shortcuts (like ~).
    """
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):
        # This is the first point we call into a remote host - make sure
        # it does not fail silently by calling a dummy run
        env.logger.info("Now, testing connection to host...")
        test = run("pwd")
        # If there is a connection failure, the rest of the code is (sometimes) not
        # reached - for example with Vagrant the program just stops after above run
        # command.
        if test != None:
            env.logger.info("Connection to host appears to work!")
        else:
            raise NotImplementedError("Connection to host failed")
        env.logger.debug("Expand paths")
        if env.has_key("local_install"):
            if not exists(env.local_install):
                run("mkdir -p %s" % env.local_install)
            with cd(env.local_install):
                result = run("pwd")
                env.local_install = result
