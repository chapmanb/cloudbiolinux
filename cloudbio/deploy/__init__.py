#!/usr/bin/env python

import os

from tempfile import tempdir
from subprocess import call

from cloudbio.utils import _setup_logging, _configure_fabric_environment, _parse_fabricrc
from cloudbio.biodata.genomes import install_data, install_data_s3, install_data_rsync
from cloudbio.galaxy import _setup_galaxy_env_defaults
from cloudbio.galaxy.utils import _chown_galaxy
from cloudbio.package.deb import _apt_packages
from fabfile import _perform_install, _install_custom


from .cloudman import cloudman_launch, bundle_cloudman
from image import configure_MI
from tools import install_tools, purge_tools
from galaxy import setup_galaxy, refresh_galaxy, seed_database, seed_workflows, wait_for_galaxy, purge_galaxy
from .util import sudoers_append, wget, eval_template
from .volume import attach_volumes, make_snapshots, sync_cloudman_bucket, detach_volumes

from fabric.main import load_settings
from fabric.api import put, run, env, settings, sudo, cd, get
from fabric.context_managers import prefix
from fabric.colors import red

try:
    from vmlauncher.transfer import FileTransferManager
    from vmlauncher import build_vm_launcher
except ImportError:
    vmlauncher = None


DEFAULT_CLOUDBIOLINUX_TARGET = None
DEFAULT_CLOUDBIOLINUX_FLAVOR = None


def deploy(options):
    actions = _expand_actions(options.get("actions"))
    if options["vm_provider"] == "novm":
        vm_launcher = LocalVmLauncher(options)
    else:
        if not vmlauncher:
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

    if _do_perform_action("bundle_cloudman", actions):
        bundle_cloudman(options)

    if _do_perform_action("sync_cloudman_bucket", actions):
        sync_cloudman_bucket(vm_launcher, options)

    if _do_perform_action("cloudman_launch", actions):
        cloudman_launch(vm_launcher, options)

    # Do we have remaining actions requiring an vm.
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
        return self.options["hostname"]

    def get_key_file(self):
        return None

    def boot_and_connect(self):
        pass

    def destroy(self):
        pass

    def get_user(self):
        pass

    def list(self):
        return []

def _setup_vm(options, vm_launcher, actions):
    destroy_on_complete = get_boolean_option(options, 'destroy_on_complete', False)
    use_galaxy = get_boolean_option(options, 'use_galaxy', False)
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
            do_refresh_galaxy = get_boolean_option(options, 'refresh_galaxy', False)
            if do_refresh_galaxy:
                refresh_galaxy(env.galaxy_repository)
            if use_galaxy:
                copy_runtime_properties(ip, options)
            if 'transfer' in actions:
                transfer_files(options)
                # Upload local compressed genomes to the cloud image, obsecure option.
                do_upload_genomes = get_boolean_option(options, 'upload_genomes', False)
                if do_upload_genomes:
                    upload_genomes(options)
            if not _seed_at_configure_time(options) and use_galaxy:
                seed_database()
                seed_workflows(options)
            if 'transfer' in actions and use_galaxy:
                wait_for_galaxy()
                create_data_library_for_uploads(options)
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
                print 'Your Galaxy instance (%s) is waiting at http://%s' % (vm_launcher.uuid, ip)
    finally:
        if destroy_on_complete:
            vm_launcher.destroy()


def _expand_actions(actions):
    unique_actions = set()
    for simple_action in ["list",
                          "destroy",
                          "transfer",
                          "purge_galaxy",
                          "setup_galaxy",
                          "purge_tools",
                          "setup_tools",
                          "setup_biodata",
                          "purge_biodata",
                          "purge_genomes",  # *_genomes are deprecated galaxy-vm-launcher ways of handling biodata.
                          "setup_genomes",
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
                          # Special CloudMan tasks not tied to a particular remote instance
                          "cloudman_launch",
                          "sync_cloudman_bucket",
                          "bundle_cloudman",
                          ]:
        if simple_action in actions:
            unique_actions.add(simple_action)
    compound_actions = {"configure": ["setup_image", "setup_tools", "setup_genomes", "setup_galaxy", "setup_ssh_key"],
                        "reinstall_galaxy": ["purge_galaxy", "setup_galaxy"],
                        "reinstall_genomes": ["purge_genomes", "setup_genomes"],
                        "reinstall_biodata": ["purge_biodata", "setup_biodata"],
                        "reinstall_tools": ["purge_tools", "setup_tools"]}
    for compound_action in compound_actions.keys():
        if compound_action in actions:
            for compound_action_part in compound_actions[compound_action]:
                unique_actions.add(compound_action_part)
    return unique_actions


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
    _setup_logging(env)
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


def configure_smtp(options):
    if 'smtp_server' in options:
        smtp_server = options['smtp_server']
        username = options['smtp_user']
        password = options['smtp_password']
        conf_file_contents = """mailhub=%s
UseSTARTTLS=YES
AuthUser=%s
AuthPass=%s
FromLineOverride=YES
""" % (smtp_server, username, password)
        _apt_packages(pkg_list=["ssmtp"])
        sudo("""echo "%s" > /etc/ssmtp/ssmtp.conf""" % conf_file_contents)
        aliases = """root:%s:%s
galaxy:%s:%s
%s:%s:%s""" % (username, smtp_server, username, smtp_server, env.user, username, smtp_server)
        sudo("""echo "%s" > /etc/ssmtp/revaliases""" % aliases)


def configure_sudoers(options):
    if "sudoers_additions" in options:
        for addition in options["sudoers_additions"]:
            sudoers_append(addition)


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
        install_proc(options["genomes"], ["s3", "raw"])
    else:
        install_proc(options["genomes"])

## Deprecated galaxy-vm-launcher way of setting up biodata.
def setup_genomes(options):
    install_proc = install_data
    sudo("mkdir -p %s" % env.data_files)
    sudo("chown -R %s:%s %s" % (env.user, env.user, env.data_files))
    put("config/tool_data_table_conf.xml", "%s/tool_data_table_conf.xml" % env.galaxy_home)
    indexing_packages = ["bowtie", "bwa", "samtools"]
    path_extensions = ":".join(map(lambda package: "/opt/galaxyTools/tools/%s/default" % package, indexing_packages))
    with prefix("PATH=$PATH:%s" % path_extensions):
        if 'S3' == options['genome_source']:
            install_proc = install_data_s3
        install_proc(options["genomes"])
    if options.get("setup_taxonomy_data", False):
        setup_taxonomy_data()
    stash_genomes_where = get_main_options_string(options, "stash_genomes")
    if stash_genomes_where:
        stash_genomes(stash_genomes_where)


def setup_taxonomy_data():
    taxonomy_directory = os.path.join(env.data_files, "taxonomy")
    env.safe_sudo("mkdir -p '%s'" % taxonomy_directory, user=env.user)
    with cd(taxonomy_directory):
        taxonomy_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
        gi_taxid_nucl = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz"
        gi_taxid_prot = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz"
        wget(taxonomy_url)
        wget(gi_taxid_nucl)
        wget(gi_taxid_prot)
        run("gunzip -c taxdump.tar.gz | tar xvf -")
        run("gunzip gi_taxid_nucl.dmp.gz")
        run("gunzip gi_taxid_prot.dmp.gz")
        run("cat gi_taxid_nucl.dmp gi_taxid_prot.dmp > gi_taxid_all.dmp")
        run("sort -n -k 1 gi_taxid_all.dmp > gi_taxid_sorted.txt")
        run("rm gi_taxid_nucl.dmp gi_taxid_prot.dmp gi_taxid_all.dmp")
        run("cat names.dmp | sed s/[\\(\\)\\'\\\"]/_/g > names.temporary")
        run("mv names.dmp names.dmp.orig")
        run("mv names.temporary names.dmp")


def configure_instance(options, actions):
    if "setup_image" in actions:
        _configure_package_holds(options)
        configure_MI(env)
        configure_smtp(options)
        configure_sudoers(options)
    if "install_biolinux" in actions:
        install_biolinux(options)
    if "install_custom" in actions:
        install_custom(options)
    if "purge_tools" in actions:
        purge_tools()
    if "setup_tools" in actions:
        install_tools(options["tools"])
    if "purge_genomes" in actions or "purge_biodata" in actions:
        purge_genomes()
    if "setup_biodata" in actions:
        setup_biodata(options)
    if "setup_genomes" in actions:
        setup_genomes(options)
    if "purge_galaxy" in actions:
        purge_galaxy()
    if "setup_galaxy" in actions:
        seed = _seed_at_configure_time(options)
        setup_galaxy(options, seed=seed)
        if seed:
            seed_workflows(options)
    if "setup_ssh_key" in actions:
        configure_ssh_key(options)


def install_custom(options):
    package = options.get("package")
    _install_custom(package)


def install_biolinux(options):
    flavor = options.get("flavor", DEFAULT_CLOUDBIOLINUX_FLAVOR)
    target = options.get("target", DEFAULT_CLOUDBIOLINUX_TARGET)
    _perform_install(target=target, flavor=flavor)


def _indices_dir_name():
    indices_dir = env.data_files
    if indices_dir.endswith("/"):
        indices_dir = indices_dir[0:(len(indices_dir) - 1)]
    indices_dir_name = os.path.basename(indices_dir)
    return indices_dir_name


def _configure_package_holds(options):
    # No longer respected. TODO: Implement.
    if 'package_holds' in options:
        env.package_holds = options['package_holds']
    else:
        env.package_holds = None


def _cd_indices_parent():
    return cd(_indices_parent())


def _indices_parent():
    parent_dir = os.path.abspath(os.path.join(env.data_files, ".."))
    return parent_dir


def _interactive_ssh(vm_launcher):
    """ Launch an interactive SSH session to host described by vm_launcher object.
    """
    host = vm_launcher.get_ip()
    user = vm_launcher.get_user()
    key_file = vm_launcher.get_key_file()
    cmd = "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -i '%s' -l '%s' '%s'" % (key_file, user, host)
    call(cmd, shell=True)


def stash_genomes(where):
    with _cd_indices_parent():
        sudo("chown %s:%s ." % (env.user, env.user))
        indices_dir_name = _indices_dir_name()
        remote_compressed_indices = "%s.tar.gz" % indices_dir_name
        run("tar czvf %s %s" % (remote_compressed_indices, indices_dir_name))
        if where == 'download':
            get(remote_path=remote_compressed_indices,
                local_path="compressed_genomes.tar.gz")
        elif where == 'opt':
            sudo("cp %s /opt/compressed_genomes.tar.gz" % remote_compressed_indices)
        else:
            print(red("Invalid option specified for stash_genomes [%s] - valid values include download and opt." % where))


def upload_genomes(options):
    with _cd_indices_parent():
        sudo("chown %s:%s ." % (env.user, env.user))
        indices_dir_name = _indices_dir_name()
        _transfer_genomes(options)
        run("rm -rf %s" % indices_dir_name)
        run("tar xzvfm compressed_genomes.tar.gz")
        sudo("/etc/init.d/galaxy restart")


def transfer_files(options):
    transfer_options = _build_transfer_options(options, "/mnt/uploaded_data", "galaxy")
    _do_transfer(transfer_options, options.get("files", []), options.get("compressed_files", []))


def _transfer_genomes(options):
    # Use just transfer settings in YAML
    options = options['transfer']
    transfer_options = _build_transfer_options(options, _indices_parent(), env.user)
    transfer_options["compress"] = False
    _do_transfer(transfer_options, ["compressed_genomes.tar.gz"])


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
    if not vmlauncher:
        raise ImportError("Require vmlauncher: https://github.com/jmchilton/vm-launcher")
    FileTransferManager(**transfer_options).transfer_files(files, compressed_files)


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


def create_data_library_for_uploads(options):
    with cd(os.path.join(env.galaxy_home, "scripts", "api")):
        db_key_arg = get_main_options_string(options, 'db_key')
        transfer_history_name = get_main_options_string(options, 'transfer_history_name')
        transfer_history_api_key = get_main_options_string(options, 'transfer_history_api_key')
        cmd_template = 'python handle_uploads.py --api_key="%s" --db_key="%s" --history="%s" --history_api_key="%s" '
        galaxy_data = options["galaxy"]
        admin_user_api_key = galaxy_data["users"][0]["api_key"]
        cmd = cmd_template % (admin_user_api_key, db_key_arg, transfer_history_name, transfer_history_api_key)
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; %s'" % cmd, user="galaxy")


def copy_runtime_properties(fqdn, options):
    runtime_properties_raw = options.get("runtime_properties", {})
    runtime_properties = {"FQDN": fqdn}
    for runtime_property_raw in runtime_properties_raw:
        (name, value) = runtime_property_raw.split(":")
        runtime_properties[name] = value
    export_file = ""
    for (name, value) in runtime_properties.iteritems():
        export_file = "export %s=%s\n%s" % (name, value, export_file)
    sudo('mkdir -p %s' % env.galaxy_home)
    _chown_galaxy(env, env.galaxy_home)
    sudo("echo '%s' > %s/runtime_properties" % (export_file, env.galaxy_home), user=env.galaxy_user)


def _seed_at_configure_time(options):
    if 'seed_galaxy' in options:
        return options['seed_galaxy'] == 'configure'
    else:
        return True
