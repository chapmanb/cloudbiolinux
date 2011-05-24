"""A Flavor reflects a specialization of a base install, the default being BioLinux.

If you want to create a new specialization (say for your own server), the
recommended procedure is to choose an existing base install (Edition) and write
a Flavor. When your Flavor is of interest to other users, it may be a good idea
to commit it to the main project.

Other flavors can be found in this directory.
"""

class Flavor:
    """Base class. Every flavor derives from this
    """
    def __init__(self, env):
        self.name = "BioLinux base Flavor"
        self.env = env
