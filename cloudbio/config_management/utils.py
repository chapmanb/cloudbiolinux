from tempfile import mkdtemp
import os
from fabric.api import settings, local, put, sudo, cd
from fabric.contrib import files


def config_dir(relative_path):
    cloudbiolinux_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
    return os.path.join(cloudbiolinux_dir, "config", relative_path)


def build_properties(env, prefix, overrides={}):
    # Prefix will be either chef or puppet
    prefix = "%s_" % prefix
    # Clone fresh dictonary to modify
    overrides = dict(overrides)

    # Load fabric environment properties into properties.
    for key, value in env.iteritems():
        # Skip invalid properties.
        if key in overrides or not isinstance(value, str):
            continue

        if key.startswith(prefix):
            # If a property starts with chef_ assume it is meant for chef and
            # add without this prefix. So chef_apache_dir would be available
            # as apache_dir.
            overrides[key[len(prefix):]] = value
        else:
            # Otherwise, allow chef to access property anyway but prefix with
            # cloudbiolinux_ so it doesn't clash with anything explicitly
            # configured for chef.
            overrides["cloudbiolinux_%s" % key] = value
    return overrides


def upload_config(config, config_folder_names=[], config_files={}):
    """ Common code to upload puppet and chef config files
    to remote server.

    Heavily based on upload procedure from fabric-provision:
    https://github.com/caffeinehit/fabric-provision/blob/master/provision/__init__.py
    """
    names = config_folder_names + config_files.keys()
    ctx = dict(map(lambda name: (name, '%s/%s' % (config.path, name)), names))

    tmpfolder = mkdtemp()

    listify = lambda what: what if isinstance(what, list) else [what]

    for folder_name in config_folder_names:
        setattr(config, folder_name, listify(getattr(config, folder_name)))

    for folder_name in config_folder_names:
        local('mkdir %s/%s' % (tmpfolder, folder_name))

    def copyfolder(folder, what):
        if not os.path.exists(folder):
            os.makedirs(folder)

        with settings(warn_only=True):
            local('cp -r %(folder)s/* %(tmpfolder)s/%(what)s' % dict(
                    folder=folder,
                    tmpfolder=tmpfolder,
                    what=what))

    for what in config_folder_names:
        map(lambda f: copyfolder(f, what), getattr(config, what))

    folder_paths = " ".join(map(lambda folder_name: "./%s" % folder_name, config_folder_names))
    local('cd %s && tar -f config_dir.tgz -cz %s' % (tmpfolder, folder_paths))

    # Get rid of old files
    with settings(warn_only=True):
        map(lambda what: sudo("rm -rf '%s'" % ctx[what]), ctx.keys())

    # Upload
    put('%s/config_dir.tgz' % tmpfolder, config.path, use_sudo=True)

    with cd(config.path):
        sudo('tar -xf config_dir.tgz')

    for file, contents in config_files.iteritems():
        files.append(ctx[file], contents, use_sudo=True)
