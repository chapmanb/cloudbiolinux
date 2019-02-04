"""
Deployer plugin containing actions related to older galaxy-vm-launcher functionality.
"""
from __future__ import print_function

import os
import time

from cloudbio.biodata.genomes import install_data, install_data_s3
from cloudbio.deploy import get_main_options_string, _build_transfer_options, _do_transfer, transfer_files, get_boolean_option
from cloudbio.deploy.util import wget, start_service, ensure_can_sudo_into, sudoers_append
from cloudbio.galaxy.utils import _chown_galaxy
from cloudbio.galaxy.tools import _setup_install_dir
from cloudbio.custom.galaxy import install_galaxy_webapp
from cloudbio.galaxy import _setup_users, _setup_xvfb, _install_nginx_standalone, _setup_postgresql
from cloudbio.package import _configure_and_install_native_packages
from cloudbio.package.deb import _apt_packages


from fabric.api import put, run, env, sudo, get, cd
from fabric.context_managers import prefix
from fabric.contrib.files import append, contains, exists


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
    """
    Setup up taxonomy data required by Galaxy. Need to find another place to put
    this, it is useful.
    """
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
            print("Invalid option specified for stash_genomes [%s] - valid values include download and opt." % where)


def upload_genomes(options):
    with _cd_indices_parent():
        sudo("chown %s:%s ." % (env.user, env.user))
        indices_dir_name = _indices_dir_name()
        _transfer_genomes(options)
        run("rm -rf %s" % indices_dir_name)
        run("tar xzvfm compressed_genomes.tar.gz")
        sudo("/etc/init.d/galaxy restart")


def purge_genomes():
    sudo("rm -rf %s" % env.data_files)


def _cd_indices_parent():
    return cd(_indices_parent())


def _indices_parent():
    parent_dir = os.path.abspath(os.path.join(env.data_files, ".."))
    return parent_dir


def _indices_dir_name():
    indices_dir = env.data_files
    if indices_dir.endswith("/"):
        indices_dir = indices_dir[0:(len(indices_dir) - 1)]
    indices_dir_name = os.path.basename(indices_dir)
    return indices_dir_name


def galaxy_transfer(vm_launcher, options):
    transfer_files(options)
    # Upload local compressed genomes to the cloud image, obsecure option.
    do_upload_genomes = get_boolean_option(options, 'upload_genomes', False)
    if do_upload_genomes:
        upload_genomes(options)
    if not _seed_at_configure_time(options):
        seed_database()
        seed_workflows(options)
    wait_for_galaxy()
    create_data_library_for_uploads(options)


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


def _seed_at_configure_time(options):
    if 'seed_galaxy' in options:
        return options['seed_galaxy'] == 'configure'
    else:
        return True


def copy_runtime_properties(vm_launcher, options):
    fqdn = vm_launcher.get_ip()
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


def _transfer_genomes(options):
    # Use just transfer settings in YAML
    options = options['transfer']
    transfer_options = _build_transfer_options(options, _indices_parent(), env.user)
    transfer_options["compress"] = False
    _do_transfer(transfer_options, ["compressed_genomes.tar.gz"])


def wait_for_galaxy():

    while not "8080" in run("netstat -lant"):
        # Check if galaxy has started
        print("Waiting for galaxy to start.")
        time.sleep(10)


def purge_galaxy():
    sudo("/etc/init.d/galaxy stop")
    sudo("rm -rf %s" % env.galaxy_home)
    init_script = "postgresql"
    # if env.postgres_version[0] < '9':
    #    # Postgres 8.4 had different name for script
    #    init_script = "postgresql-%s" % env.postgres_version
    sudo("/etc/init.d/%s restart" % init_script)
    sudo('psql  -c "drop database galaxy;"', user="postgres")
    sudo('psql  -c "create database galaxy;"', user="postgres")


def setup_galaxy(options):
    seed = _seed_at_configure_time(options)
    setup_galaxy(options, seed=seed)
    if seed:
        seed_workflows(options)


def _setup_galaxy(options, seed=True):
    """Deploy a Galaxy server along with some tools.
    """
    _setup_install_dir(env)  # Still needed? -John
    install_galaxy_webapp(env)
    #_fix_galaxy_permissions()
    _setup_shed_tools_dir()
    _setup_galaxy_log_dir()
    _migrate_galaxy_database()
    if seed:
        seed_database(options["galaxy"])
    _start_galaxy()


def _migrate_galaxy_database():
    with cd(env.galaxy_home):
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; python ./scripts/build_universe_config.py conf.d; python -ES ./scripts/fetch_eggs.py; ./create_db.sh'", user="galaxy")


def seed_database(galaxy_data):
    with cd(env.galaxy_home):
        sudo("rm -f seed.py")
        _setup_database_seed_file(galaxy_data)
        sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; python ./scripts/build_universe_config.py conf.d; python -ES ./scripts/fetch_eggs.py; python seed.py'", user="galaxy")


def seed_workflows(options):
    wait_for_galaxy()
    galaxy_data = options["galaxy"]
    with cd(os.path.join(env.galaxy_home, "workflows")):
        for user in galaxy_data["users"]:
            api_key = user["api_key"]
            workflows = None
            if "workflows" in user:
                workflows = user["workflows"]
            if not workflows:
                continue
            for workflow in workflows:
                sudo("bash -c 'export PYTHON_EGG_CACHE=eggs; bash import_all.sh %s %s'" % (api_key, workflow), user=env.galaxy_user)


def _setup_database_seed_file(galaxy_data):
    _seed_append("""from scripts.db_shell import *
from galaxy.util.bunch import Bunch
from galaxy.security import GalaxyRBACAgent
bunch = Bunch( **globals() )
bunch.engine = engine
# model.flush() has been removed.
bunch.session = db_session
# For backward compatibility with "model.context.current"
bunch.context = db_session
security_agent = GalaxyRBACAgent( bunch )
security_agent.sa_session = sa_session

def add_user(email, password, key=None):
    query = sa_session.query( User ).filter_by( email=email )
    if query.count() > 0:
        return query.first()
    else:
        user = User(email)
        user.set_password_cleartext(password)
        sa_session.add(user)
        sa_session.flush()

        security_agent.create_private_user_role( user )
        if not user.default_permissions:
            security_agent.user_set_default_permissions( user, history=True, dataset=True )

        if key is not None:
            api_key = APIKeys()
            api_key.user_id = user.id
            api_key.key = key
            sa_session.add(api_key)
            sa_session.flush()
        return user

def add_history(user, name):
    query = sa_session.query( History ).filter_by( user=user ).filter_by( name=name )
    if query.count() == 0:
        history = History(user=user, name=name)
        sa_session.add(history)
        sa_session.flush()
        return history
    else:
        return query.first()

""")
    i = 0
    for user in galaxy_data["users"]:
        username = user["username"]
        password = user["password"]
        api_key = user["api_key"]
        histories = None
        if "histories" in user:
            histories = user["histories"]
        user_object = "user_%d" % i
        _seed_append("""%s = add_user("%s", "%s", "%s")""" % (user_object, username, password, api_key))
        _import_histories(user_object, histories)
        i = i + 1


def _import_histories(user_object, histories):
    if not histories:
        return
    for history_name in histories:
        _import_history(user_object, history_name)


def _import_history(user_object, history_name):
    history_name_stripped = history_name.strip()
    if history_name_stripped:
        _seed_append("""add_history(%s, "%s")""" % (user_object, history_name_stripped))


def _seed_append(text):
    append("%s/seed.py" % env.galaxy_home, text, use_sudo=True)


def _start_galaxy():
    # Create directory to store galaxy service's pid file.
    _make_dir_for_galaxy("/var/lib/galaxy")
    start_service("galaxy")


def refresh_galaxy(target_galaxy_repo):
    _update_galaxy(target_galaxy_repo)
    sudo("/etc/init.d/galaxy restart", pty=False)


def _setup_galaxy_log_dir():
    _make_dir_for_galaxy("/var/log/galaxy")


def _setup_shed_tools_dir():
    _make_dir_for_galaxy("%s/../shed_tools" % env.galaxy_home)


def _make_dir_for_galaxy(path):
    sudo("mkdir -p '%s'" % path)
    _chown_galaxy(env, path)


def _update_galaxy(target_galaxy_repo):
    # Need to merge? -John
    hg_command = "hg pull %s; hg update" % target_galaxy_repo
    with cd(env.galaxy_home):
        sudo(hg_command, user=env.galaxy_user)


def refresh_galaxy_action(vm_launcher, options):
    refresh_galaxy(env.galaxy_repository)


def setup_image(options):
    _configure_package_holds(options)
    configure_MI(env)
    configure_smtp(options)
    configure_sudoers(options)


def _configure_package_holds(options):
    # No longer respected. TODO: Implement.
    if 'package_holds' in options:
        env.package_holds = options['package_holds']
    else:
        env.package_holds = None


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


def configure_MI(env):
    # Clean this next line up.
    _configure_and_install_native_packages(env, ["minimal", "cloudman", "galaxy"])
    # _update_system()
    _setup_users(env)
    _setup_xvfb(env)
    _required_programs(env)


# == required programs
def _required_programs(env):
    """ Install required programs """
    if not exists(env.install_dir):
        sudo("mkdir -p %s" % env.install_dir)
        sudo("chown %s %s" % (env.user, env.install_dir))

    # Setup global environment for all users
    install_dir = os.path.split(env.install_dir)[0]
    exports = ["export PATH=%s/bin:%s/sbin:$PATH" % (install_dir, install_dir),
               "export LD_LIBRARY_PATH=%s/lib" % install_dir]
    for e in exports:
        _ensure_export(e)
    # Install required programs
    _install_nginx_standalone(env)
    _start_nginx(env)
    _deploy_setup_postgresql(env)

    # Verify this is not needed.
    # _install_samtools()


def _ensure_export(command):
    if not contains('/etc/bash.bashrc', command):
        append('/etc/bash.bashrc', command, use_sudo=True)


def _start_nginx(env):
    galaxy_data = env.galaxy_data_mount
    env.safe_sudo("mkdir -p '%s'" % env.galaxy_data)
    _chown_galaxy(env, galaxy_data)
    start_service("nginx")


def _deploy_setup_postgresql(env):
    ensure_can_sudo_into("postgres")
    _setup_postgresql(env)


configure_actions = {"setup_image": setup_image,
                     "setup_genomes": setup_genomes,
                     "purge_genomes": purge_genomes,
                     "setup_galaxy": setup_galaxy,
                     "purge_galaxy": purge_galaxy,
                    }

ready_actions = {"galaxy_transfer": galaxy_transfer,
                 "refresh_galaxy": refresh_galaxy_action,
                 "copy_runtime_properties": copy_runtime_properties,
}

compound_actions = {"configure": ["setup_image", "setup_tools", "setup_genomes", "setup_galaxy", "setup_ssh_key"],
                    "reinstall_galaxy": ["purge_galaxy", "setup_galaxy"],
                    "reinstall_genomes": ["purge_genomes", "setup_genomes"],
                    "reinstall_tools": ["purge_tools", "setup_tools"]
}
