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
        self.env.logger.debug("check_packages_source not implemented")
