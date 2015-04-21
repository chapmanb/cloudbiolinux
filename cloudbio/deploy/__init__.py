#!/usr/bin/env python

import os

from tempfile import tempdir
from subprocess import call
from inspect import getargspec

from cloudbio.utils import _setup_logging, _configure_fabric_environment, _parse_fabricrc
from cloudbio.biodata.genomes import install_data, install_data_s3, install_data_rsync
from cloudbio.galaxy import _setup_galaxy_env_defaults
from cloudbio.galaxy.utils import _chown_galaxy
from cloudbio.galaxy.tools import _install_tools

from fabfile import _perform_install, _install_custom

from .util import eval_template
from .volume import attach_volumes, make_snapshots, detach_volumes

import cloudbio.deploy.plugins

from fabric.main import load_settings
from fabric.api import put, run, env, settings, sudo


try:
    from .vmlauncher.transfer import FileTransferManager
    from .vmlauncher import build_vm_launcher
except ImportError:
    build_vm_launcher = None
    FileTransferManager = None

DEFAULT_CLOUDBIOLINUX_TARGET = None
DEFAULT_CLOUDBIOLINUX_FLAVOR = None


def deploy(options):
    _setup_logging(env)

    actions = _expand_actions(options.get("actions"))
    if options["vm_provider"] == "novm":
        vm_launcher = LocalVmLauncher(options)
    else:
        if not build_vm_launcher:
            raise ImportError("Require vmlauncher: https://github.com/jmchilton/vm-launcher")
        vm_launcher = build_vm_launcher(options)

    if _do_perform_action("list", actions):
        for node in vm_launcher.list():
            print "Active node with uuid %s <%s>" % (node.uuid, node)

    if _do_perform_action("destroy", actions):
        target_name = options["hostname"]
        for node in vm_launcher.list():
            node_name = node.name
            if node_name == target_name:
                vm_launcher.destroy(node)

    __invoke_plugin_actions(env, actions, "local_actions", [vm_launcher, options])

    # Do we have remaining actions requiring an vm?
    if len(actions) > 0:
        print 'Setting up virtual machine'
        vm_launcher.boot_and_connect()
        _setup_vm(options, vm_launcher, actions)


class LocalVmLauncher:
    """Provide a lightweight real machine, non-vm class for launching.
    """
    def __init__(self, options):
        self.options = options

    def get_ip(self):
        specified_hostname = self.options.get("hostname", None)
        hostname = specified_hostname or "localhost"
        return hostname

    def get_key_file(self):
        return None

    def boot_and_connect(self):
        pass

    def destroy(self):
        pass

    def get_user(self):
        return env.user

    def list(self):
        return []


def _setup_vm(options, vm_launcher, actions):
    destroy_on_complete = get_boolean_option(options, 'destroy_on_complete', False)
    try:
        ip = vm_launcher.get_ip()
        _setup_fabric(vm_launcher, ip, options)
        with settings(host_string=ip):
            _setup_cloudbiolinux(options)
            if 'attach_volumes' in actions:
                attach_volumes(vm_launcher, options)
            if 'max_lifetime' in options:
                seconds = options['max_lifetime']
                # Unclear why the sleep is needed, but seems to be otherwise
                # this doesn't work.
                run("bash -c 'nohup sudo shutdown -h %d &'; sleep 2" % seconds)

            configure_instance(options, actions)

            if 'transfer' in actions:
                transfer_files(options)

            __invoke_plugin_actions(env, actions, "ready_actions", [vm_launcher, options])

            if 'ssh' in actions:
                _interactive_ssh(vm_launcher)
            if 'attach_ip' in actions:
                vm_launcher.attach_public_ip()
            if 'snapshot_volumes' in actions:
                make_snapshots(vm_launcher, options)
            if 'detach_volumes' in actions:
                detach_volumes(vm_launcher, options)
            if 'package' in actions:
                name_template = vm_launcher.package_image_name()
                name = eval_template(env, name_template)
                vm_launcher.package(name=name)
            if not destroy_on_complete and hasattr(vm_launcher, "uuid"):
                print 'Your instance (%s) is waiting at http://%s' % (vm_launcher.uuid, ip)
    finally:
        if destroy_on_complete:
            vm_launcher.destroy()


def _expand_actions(actions):
    unique_actions = set()
    for simple_action in _possible_actions():
        if simple_action in actions:
            unique_actions.add(simple_action)
    compound_actions = __get_plugin_actions(env, "compound_actions")
    for compound_action in compound_actions.keys():
        if compound_action in actions:
            for compound_action_part in compound_actions[compound_action]:
                unique_actions.add(compound_action_part)
    return unique_actions


def _possible_actions():
    possible_actions = [ "list",
                         "destroy",
                         "transfer",
                         "purge_tools",
                         "setup_tools",
                         "setup_biodata",
                         "setup_ssh_key",
                         "package",
                         "setup_image",
                         "launch",  # Dummy action justs launches image
                         "install_biolinux",
                         "install_custom",
                         "ssh",
                         "attach_ip",
                         "snapshot_volumes",
                         "attach_volumes",
                         "detach_volumes",
                        ]
    for action_type in ["local_actions", "configure_actions", "ready_action"]:
        for action in  __get_plugin_actions(env, action_type):
            possible_actions.append(action)
    return possible_actions


def _do_perform_action(action, action_list):
    do_perform = action in action_list
    if do_perform:
        action_list.remove(action)
    return do_perform


def _setup_fabric(vm_launcher, ip, options):
    env.user = vm_launcher.get_user()
    env.hosts = [ip]
    env.key_filename = vm_launcher.get_key_file()
    env.disable_known_hosts = True


def _setup_cloudbiolinux(options):
    def fabricrc_loader(env):
        _setup_cloudbiolinux_fabric_properties(env, options)

    flavor = get_main_options_string(options, "flavor", DEFAULT_CLOUDBIOLINUX_FLAVOR)
    _configure_fabric_environment(env, flavor, fabricrc_loader=fabricrc_loader)
    _setup_image_user_data(env, options)


def _setup_cloudbiolinux_fabric_properties(env, options):
    fabricrc_file = get_main_options_string(options, "fabricrc_file", None)
    env.config_dir = os.path.join(os.path.dirname(__file__), "..", "..", "config")
    env.tool_data_table_conf_file = os.path.join(env.config_dir, "..",
                                                 "installed_files",
                                                 "tool_data_table_conf.xml")
    if fabricrc_file:
        env.update(load_settings(fabricrc_file))
    else:
        # Let cloudbiolinux find out default file based on flavor, dist, etc...
        _parse_fabricrc(env)
    overrides = options.get("fabricrc_overrides", {})
    for key, value in overrides.iteritems():
        # yaml parses bools, wouldn't be expected coming out of a fabricrc
        # file so replace everything with a string.
        if isinstance(value, bool):
            overrides[key] = str(value)
    env.update(overrides)
    _setup_galaxy_env_defaults(env)


def _setup_image_user_data(env, options):
    if "image_user_data" in options:
        env["image_user_data_dict"] = options["image_user_data"]


def purge_genomes():
    sudo("rm -rf %s" % env.data_files)


def configure_ssh_key(options):
    if "galaxy_ssh_key" in options:
        key_file = options["galaxy_ssh_key"]
        sudo("mkdir -p /home/%s/.ssh" % (env.galaxy_user))
        sudo("chmod 700 /home/%s/.ssh" % (env.galaxy_user))
        put(local_path=key_file,
            remote_path="/home/%s/.ssh/%s" % (env.galaxy_user, os.path.basename(key_file)),
            use_sudo=True,
            mode=0600)
        _chown_galaxy(env, "/home/%s/.ssh" % env.galaxy_user)


def setup_biodata(options):
    install_proc = install_data
    genome_source = options.get("genome_source", "default")
    install_proc = {
        "default": install_data,
        "S3": install_data_s3,
        "rsync": install_data_rsync,
    }[genome_source]
    if genome_source == "default":
        install_proc(options["genomes"], ["ggd", "s3", "raw"])
    else:
        install_proc(options["genomes"])


def configure_instance(options, actions):
    if "install_biolinux" in actions:
        install_biolinux(options)
    if "install_custom" in actions:
        install_custom(options)
    if "purge_tools" in actions:
        purge_tools()

    __invoke_plugin_actions(env, actions, "configure_actions", [options])

    if "setup_tools" in actions:
        install_tools(options["tools"])
    if "setup_biodata" in actions:
        setup_biodata(options)
    if "setup_ssh_key" in actions:
        configure_ssh_key(options)


def install_custom(options):
    package = options.get("package")
    _install_custom(package)


def install_biolinux(options):
    flavor = options.get("flavor", DEFAULT_CLOUDBIOLINUX_FLAVOR)
    target = options.get("target", DEFAULT_CLOUDBIOLINUX_TARGET)
    _perform_install(target=target, flavor=flavor, more_custom_add=options.get("custom_add", None))


def _interactive_ssh(vm_launcher):
    """ Launch an interactive SSH session to host described by vm_launcher object.
    """
    host = vm_launcher.get_ip()
    user = vm_launcher.get_user()
    key_file = vm_launcher.get_key_file()
    cmd = "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -i '%s' -l '%s' '%s'" % (key_file, user, host)
    call(cmd, shell=True)


def transfer_files(options):
    transfer_options = _build_transfer_options(options, "/mnt/uploaded_data", "galaxy")
    _do_transfer(transfer_options, options.get("files", []), options.get("compressed_files", []))


def _build_transfer_options(options, destination, user):
    transfer_options = {}
    transfer_options['compress'] = get_boolean_option(options, 'compress_transfers', True)
    transfer_options['num_compress_threads'] = int(get_main_options_string(options, 'num_compress_threads', '1'))
    transfer_options['num_transfer_threads'] = int(get_main_options_string(options, 'num_transfer_threads', '1'))
    transfer_options['num_decompress_threads'] = int(get_main_options_string(options, 'num_decompress_threads', '1'))
    transfer_options['chunk_size'] = int(get_main_options_string(options, 'transfer_chunk_size', '0'))
    transfer_options['transfer_retries'] = int(get_main_options_string(options, 'transfer_retries', '3'))
    transfer_options['local_temp'] = get_main_options_string(options, 'local_temp_dir', tempdir)
    transfer_options['destination'] = destination
    transfer_options['transfer_as'] = user
    return transfer_options


def _do_transfer(transfer_options, files, compressed_files=[]):
    if not FileTransferManager:
        raise ImportError("Require vmlauncher: https://github.com/jmchilton/vm-launcher")
    FileTransferManager(**transfer_options).transfer_files(files, compressed_files)


def purge_tools():
    env.safe_sudo("rm -rf %s" % env.install_dir)


def install_tools(tools_conf):
    """
    """
    _install_tools(env, tools_conf)


def get_boolean_option(options, name, default=False):
    if name not in options:
        return default
    else:
        return options[name]


def get_main_options_string(options, key, default=''):
    value = default
    if key in options:
        value = options[key]
    return value


def __invoke_plugin_actions(env, actions, action_type, provided_args):
    possible_actions = __get_plugin_actions(env, action_type)
    for action in list(actions):
        if action in possible_actions:
            __invoke_plugin_action(env, possible_actions[action], provided_args)
            actions.remove(action)


def __invoke_plugin_action(env, action_function, provided_args):
    arg_spec = getargspec(action_function).args
    args = [] if not arg_spec else provided_args
    action_function(*args)


def __get_plugin_actions(env, action_type):
    actions = {}
    for plugin_module in __get_plugin_modules(env):
        if hasattr(plugin_module, action_type):
            for action_name, action_function in getattr(plugin_module, action_type).iteritems():
                actions[action_name] = action_function
    return actions


def __get_plugin_modules(env):
    if not "plugin_modules" in env:
        unsorted_module_names = __get_plugin_module_names( )
        ## Load modules in reverse order to allow hierarchical overrides
        module_names = sorted(unsorted_module_names, reverse=True)
        modules = []
        for plugin_module_name in module_names:
            try:
                module = __import__(plugin_module_name)
                for comp in plugin_module_name.split(".")[1:]:
                    module = getattr(module, comp)
                modules.append(module)
            except BaseException, exception:
                exception_str = str(exception)
                message = "%s rule module could not be loaded: %s" % (plugin_module_name, exception_str)
                env.logger.warn(message)
                continue
        env.plugin_modules = modules
    return env.plugin_modules


def __get_plugin_module_names():
    plugin_module_dir = cloudbio.deploy.plugins.__path__[0]
    names = []
    for fname in os.listdir(plugin_module_dir):
        if not(fname.startswith("_")) and fname.endswith(".py"):
            rule_module_name = "cloudbio.deploy.plugins.%s" % fname[:-len(".py")]
            names.append( rule_module_name )
    return names
