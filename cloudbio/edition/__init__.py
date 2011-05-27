"""An Edition reflects a base install, the default being BioLinux.

Editions are shared between multiple projects. To specialize an edition, create
a Flavor instead.

Other editions can be found in this directory
"""

class Edition:
    """Base class. Every edition derives from this
    """
    def __init__(self, env):
        self.name = "BioLinux base Edition"
        self.short_name = "biolinux"
        self.version = env.version
        self.env = env
        self.include_oracle_virtualbox = True
        self.include_freenx = True
        self.include_apt_automation = True
        self.include_force_upgrade = True
        self.include_hadoop = True

    def check_packages_source(self):
        """Override for check package definition file before updating
        """
        self.env.logger.debug("check_packages_source not implemented")
