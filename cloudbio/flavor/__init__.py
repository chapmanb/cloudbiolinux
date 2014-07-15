"""A Flavor reflects a specialization of a base install, the default being BioLinux.

If you want to create a new specialization (say for your own server), the
recommended procedure is to choose an existing base install (Edition) and write
a Flavor. When your Flavor is of interest to other users, it may be a good idea
to commit it to the main project (in ./contrib/flavor).

Other (main) flavors can be found in this directory and in ./contrib/flavors
"""

class Flavor:
    """Base class. Every flavor derives from this
    """
    def __init__(self, env):
        self.name = "Base Flavor - no overrides"
        self.short_name = None  # override this
        self.env = env
        self.check_distribution()

    def rewrite_config_items(self, name, items):
        """Generic hook to rewrite a list of configured items.

        Can define custom dispatches based on name: packages, custom,
        python, ruby, perl
        """
        return items

    def check_distribution(self):
        """Ensure the distribution matches an expected type for this edition.

        Base supports multiple distributions.
        """
        pass

    def check_packages_source(self):
        """Override for check package definition file before updating
        """
        pass

    def rewrite_apt_sources_list(self, sources):
        """Allows editions to modify the sources list
        """
        return sources

    def rewrite_apt_preferences(self, preferences):
        """Allows editions to modify the apt preferences policy file
        """
        return preferences

    def rewrite_apt_automation(self, package_info):
        """Allows editions to modify the apt automation list
        """
        return package_info

    def rewrite_apt_keys(self, standalone, keyserver):
        """Allows editions to modify key list"""
        return standalone, keyserver

    def apt_upgrade_system(self, env=None):
        """Upgrade system through apt - so this behaviour can be overridden
        """
        env.safe_sudo("apt-get -y --force-yes upgrade")

    def post_install(self):
        """Add scripts for starting FreeNX and CloudMan.
        """
        pass
