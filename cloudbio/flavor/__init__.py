"""A Flavor reflects a specialization of a base install, the default being BioLinux.

If you want to create a new specialization (say for your own server), the
recommended procedure is to choose an existing base install (Edition) and write
a Flavor. When your Flavor is of interest to other users, it may be a good idea
to commit it to the main project (in ./contrib/flavor).

Other (main) flavors can be found in this directory and in ./contrib/flavors
"""

from fabric.api import *

class Flavor:
    """Base class. Every flavor derives from this
    """
    def __init__(self, env):
        self.name = "Base Flavor - no overrides"
        self.env = env

    def rewrite_packages_list(self, list):
        return list

    def rewrite_python_egg_list(self, list):
        return list

    def rewrite_ruby_gem_list(self, list):
        return list

    def rewrite_perl_cpan_list(self, list):
        return list

    def rewrite_r_cran_list(self, list):
        return list

    def rewrite_custom_list(self, list):
        return list

    def post_install(self):
        """Post installation hook"""

env.flavor = Flavor(env)
