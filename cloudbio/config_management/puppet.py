from fabric.state import _AttributeDict
from fabric.api import cd

from utils import upload_config, config_dir, build_properties
from cloudbio.package.deb import _apt_packages
import os

DEFAULTS = dict(
    path='/var/puppet',
    log_level='info',
    modules=config_dir(os.path.join('puppet', 'modules'))
)

puppet = _AttributeDict(DEFAULTS)


def _puppet_provision(env, classes):
    env.safe_sudo('mkdir -p %(path)s' % puppet)
    manifest_body = "node default {\n%s\n}\n" % _build_node_def_body(env, classes)
    config_files = {"manifest.pp": manifest_body}
    upload_config(puppet, config_folder_names=["modules"], config_files=config_files)
    # TODO: Allow yum based install
    _apt_packages(pkg_list=["puppet"])
    with cd(puppet.path):
        env.safe_sudo("sudo puppet apply --modulepath=modules manifest.pp")


def _build_node_def_body(env, classes):
    contents = ""
    properties = build_properties(env, "puppet")
    contents += "\n".join(["$%s = '%s'" % (key, value.replace("'", r"\'")) for key, value in properties.iteritems()])
    contents += "\n"
    contents += "\n".join([_build_class_include(env, class_name) for class_name in classes])
    return contents


def _build_class_include(env, class_name):
    """
    If parentns::classname is included and fabric
    properties such as puppet_parentns__classname_prop = val1
    are set, the class included in puppet will be something like

    class { 'parentns::classname':
        prop => 'val1',
    }
    """
    include_def = "class { '%s': \n" % class_name
    property_prefix = _property_prefix(class_name)
    for name, value in env.iteritems():
        if name.startswith(property_prefix):
            property_name = name[len(property_prefix):]
            include_def += "  %s => '%s',\n" % (property_name, value)
    include_def += "\n}"
    return include_def


def _property_prefix(class_name):
    return "puppet_%s_" % class_name.replace("::", "__")
