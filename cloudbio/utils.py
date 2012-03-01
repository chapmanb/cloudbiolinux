"""Utilities for logging and progress tracking.
"""
import logging

from fabric.api import sudo
from fabric.colors import yellow, red, green, magenta

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
    sudo("date +\"%D %T - Updated "+info+"\" >> "+logfn)
