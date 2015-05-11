"""Fabric deployment file to install genomic data on remote instances.

Designed to automatically download and manage biologically associated
data on cloud instances like Amazon EC2.

Fabric (http://docs.fabfile.org) manages automation of remote servers.

Usage:
    fab -i key_file -H servername -f data_fabfile.py install_data
"""
import os
import sys

from fabric.main import load_settings
from fabric.api import *
from fabric.contrib.files import *
from fabric.context_managers import path
try:
    import boto
except ImportError:
    boto = None

# preferentially use local cloudbio directory
for to_remove in [p for p in sys.path if p.find("cloudbiolinux-") > 0]:
    sys.path.remove(to_remove)
sys.path.append(os.path.dirname(__file__))

from cloudbio.utils import _setup_logging, _configure_fabric_environment
from cloudbio.biodata import genomes

# -- Host specific setup

env.remove_old_genomes = False

def setup_environment():
    """Setup environment with required data file locations.
    """
    _setup_logging(env)
    _add_defaults()
    _configure_fabric_environment(env, ignore_distcheck=True)

def _add_defaults():
    """Defaults from fabricrc.txt file; loaded if not specified at commandline.
    """
    env.config_dir = os.path.join(os.path.dirname(__file__), "config")
    conf_file = "tool_data_table_conf.xml"
    env.tool_data_table_conf_file = os.path.join(os.path.dirname(__file__),
                                                 "installed_files", conf_file)
    if not env.has_key("distribution"):
        config_file = os.path.join(env.config_dir, "fabricrc.txt")
        if os.path.exists(config_file):
            env.update(load_settings(config_file))

CONFIG_FILE = os.path.join(os.path.dirname(__file__), "config", "biodata.yaml")

def install_data(config_source=CONFIG_FILE):
    """Main entry point for installing useful biological data.
    """
    setup_environment()
    genomes.install_data(config_source)

def install_data_s3(config_source=CONFIG_FILE, do_setup_environment=True):
    """Install data using pre-existing genomes present on Amazon s3.
    """
    setup_environment()
    genomes.install_data_s3(config_source)

def install_data_rsync(config_source=CONFIG_FILE):
    """Install data using Galaxy rsync data servers.
    """
    setup_environment()
    genomes.install_data_rsync(config_source)

def install_data_ggd(recipe, organism):
    """Install data using Get Genomics Data (GGD) recipes.
    """
    setup_environment()
    from cloudbio.biodata import ggd, genomes
    genome_dir = os.path.join(genomes._make_genome_dir(), organism)
    recipe_file = os.path.join(os.path.dirname(__file__), "ggd-recipes", organism, "%s.yaml" % recipe)
    ggd.install_recipe(genome_dir, recipe_file)

def upload_s3(config_source=CONFIG_FILE):
    """Upload prepared genome files by identifier to Amazon s3 buckets.
    """
    setup_environment()
    genomes.upload_s3(config_source)
