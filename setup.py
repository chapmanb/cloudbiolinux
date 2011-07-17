#!/usr/bin/env python
"""Install cloudbiolinux install libraries for external use.

This is not needed for running fabric files, but is useful for
reusing code.
"""
from setuptools import setup, find_packages

setup(name = "cloudbiolinux",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "configure virtual (or real) machines with tools for biological analyses",
      license = "MIT",
      url = "http://cloudbiolinux.org",
      packages = find_packages(),
      scripts = [],
      install_requires = [
          "PyYAML >= 3.09",
          "fabric >= 1.1.1"]
      )
