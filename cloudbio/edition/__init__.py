import logging

class Edition:
    """Base class. Every edition derives from this
    """
    def __init__(self, env):
        self.name = "BioLinux base Edition"
        self.env = env
        self.include_oracle_virtualbox = True
        self.include_freenx = True
        self.logger = logging.getLogger(__name__)

    def check_packages_source(self):
        """Override for check package definition file before updating
        """
        self.logger.debug("check_packages_source not implemented")
