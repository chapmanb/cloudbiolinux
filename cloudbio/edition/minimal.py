from fabric.api import *

from cloudbio.edition import Edition

class Minimal(Edition):
    """Minimal specialization of BioLinux
    """
    def __init__(self, env):
        # Edition.__init__(self,env) Probably not a good idea to call super from Minimal
        self.name = "Minimal Edition"
        self.short_name = "minimal"
        self.version = env.version
        self.env = env
        self.include_oracle_virtualbox = False
        self.include_freenx = False
        self.include_apt_automation = False
        self.force_upgrade = False
        self.include_hadoop = False
        self.is_ubuntu = False
        self.is_debian = False
        self.is_centos = False
        self.is_debian_derived = False
        if env.distribution == "ubuntu":
            self.is_ubuntu = True
            self.is_debian_derived = True
        elif env.distribution == "centos":
            self.centos = True
        elif env.distribution == "debian":
            self.is_debian = True
            self.is_debian_derived = True

    def check_packages_source(self):
        # Removes sources, just to be sure
        self.env.logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)

