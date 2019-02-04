import os
import json

from fabric.api import cd
from fabric.contrib import files
from fabric.state import _AttributeDict

from cloudbio.flavor.config import get_config_file
from utils import build_properties, upload_config, config_dir


# Code based heavily on fabric-provision. https://github.com/caffeinehit/fabric-provision

DEFAULTS = dict(
    path='/var/chef',
    data_bags=config_dir(os.path.join('chef', 'data_bags')),
    roles=config_dir(os.path.join('chef', 'roles')),
    cookbooks=config_dir(os.path.join('chef', 'cookbooks')),
    log_level='info',
    recipes=[],
    run_list=[],
    json={},
)

SOLO_RB = """
log_level            :%(log_level)s
log_location         STDOUT
file_cache_path      "%(path)s"
data_bag_path        "%(path)s/data_bags"
role_path            [ "%(path)s/roles" ]
cookbook_path        [ "%(path)s/cookbooks" ]
Chef::Log::Formatter.show_time = true
"""


class ChefDict(_AttributeDict):
    def add_recipe(self, recipe):
        self.run_list.append('recipe[{0}]'.format(recipe))

    def add_role(self, role):
        self.run_list.append('role[{0}]'.format(role))

    def _get_json(self):
        the_json = self['json'].copy()
        the_json['run_list'] = self['run_list']
        return the_json

    json = property(fget=_get_json)

chef = ChefDict(DEFAULTS)


def omnibus(env):
    """
    Install Chef from Opscode's Omnibus installer
    """
    ctx = {
        'filename': '%(path)s/install-chef.sh' % chef,
        'url': 'http://opscode.com/chef/install.sh',
    }
    if not files.exists(ctx['filename']):
        env.safe_sudo('wget -O %(filename)s %(url)s' % ctx)
        with cd(chef.path):
            env.safe_sudo('bash install-chef.sh')


def _chef_provision(env, _omnibus=True):
    env.safe_sudo('mkdir -p %(path)s' % chef)

    omnibus(env)

    config_files = {'node.json': json.dumps(chef.json),
                    'solo.rb': SOLO_RB % chef}
    upload_config(chef, config_folder_names=['cookbooks', 'data_bags', 'roles'], config_files=config_files)

    with cd(chef.path):
        env.safe_sudo('chef-solo -c solo.rb -j node.json')


def _configure_chef(env, chef):

    # Set node json properties
    node_json_path = get_config_file(env, "node_extra.json").base
    chef.json = _build_chef_properties(env, node_json_path)

    # Set whether to use the Opscode Omnibus Installer to load Chef.
    use_omnibus_installer_str = env.get("use_chef_omnibus_installer", "false")
    chef.use_omnibus_installer = use_omnibus_installer_str.upper() in ["TRUE", "YES"]


def _build_chef_properties(env, config_file):
    """
    Build python object representation of the Chef-solo node.json file from
    node_extra.json in config dir and the fabric environment.
    """

    json_properties = _parse_json(config_file)
    return build_properties(env, "chef", json_properties)


def _parse_json(filename):
    """ Parse a JSON file
        First remove comments and then use the json module package
        Comments look like :
            // ...
    """
    with open(filename) as f:
        lines = f.readlines()
        content = ''.join([line for line in lines if not line.startswith('//')])
        return json.loads(content)
