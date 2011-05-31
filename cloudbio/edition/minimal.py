from fabric.api import *

from cloudbio.edition.base import Edition

class Minimal(Edition):
    """Minimal specialization of BioLinux
    """
    def __init__(self, env):
        # Edition.__init__(self,env) Probably not a good idea to call super from Minimal
        self.name = "Minimal Edition"
        self.short_name = "minimal"
        self.version = env.version
        self.env = env

    def check_packages_source(self):
        # Removes sources, just to be sure
        self.env.logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)

    def rewrite_apt_sources_list(self, sources):
        """Allows editions to modify the sources list. Minimal only
           uses the default packages
        """
        return []

    def rewrite_apt_automation(self, list):
        """Allows editions to modify the apt automation list
        """
        return []

    def rewrite_apt_keys(self, list):
        """Allows editions to modify key list"""
        return []

    def rewrite_apt_keyserver(self, list):
        """Allows editions to modify key list"""
        return []

    def apt_upgrade_system(self):
        """Do nothing"""
        env.logger.debug("Skipping forced system upgrade")

