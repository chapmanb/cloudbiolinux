import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
logger.addHandler(ch)


class Edition:
    """Base class. Every edition derives from this
    """
    def __init__(self, env):
        self.name = "BioLinux base Edition"
        self.env = env
        self.include_oracle_virtualbox = True
        self.include_freenx = True

    def check_packages_source(self):
        """Override for check package definition file before updating
        """
        logger.debug("check_packages_source not implemented")
