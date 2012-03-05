from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

from cloudbio.custom.shared import (_fetch_and_unpack)

import sys

class SealFlavor(Flavor):
	"""A flavour for installing Seal
	"""
	def __init__(self, env):
		Flavor.__init__(self,env)
		self.name = "Seal Flavor"

	def rewrite_config_items(self, name, packages):
		if name == 'packages':
			env.logger.info("Hello from SealFlavor.  Package list:\n")
			for package in packages:
				env.logger.info("Selected: "+name+" "+package)
			if sys.version_info < (2,7):
				# for versions of Python prior to 2.7 we need to add importlib
				# and argparse
				packages += [ 
					"python-importlib",
					"python-argparse"
				]
		return packages

	def post_install(self):
		env.logger.info("Starting post-install")
		pass

env.flavor = SealFlavor(env)
