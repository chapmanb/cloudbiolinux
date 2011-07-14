"""Build instructions associated with CloudMan.

http://wiki.g2.bx.psu.edu/Admin/Cloud

Adapted from Enis Afgan's code: https://bitbucket.org/afgane/mi-deployment
"""

cm_upstart = """
description     "Start CloudMan contextualization script"

start on runlevel [2345]

task
exec python %s/ec2autorun.py
"""
import os

REPO_ROOT_URL = "https://bitbucket.org/afgane/mi-deployment/raw/tip"

def _configure_ec2_autorun(env, use_repo_autorun=False):
    script = "ec2autorun.py"
    remote = os.path.join(env.install_dir, script)
    if use_repo_autorun:
        url = os.path.join(REPO_ROOT_URL, script)
        sudo("wget --output-document=%s %s" % (remote, url))
    else:
        install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
        put(os.path.join(install_file_dir, script), remote, mode=0777, use_sudo=True)
    # Create upstart configuration file for boot-time script
    cloudman_boot_file = 'cloudman.conf'
    with open( cloudman_boot_file, 'w' ) as f:
        print >> f, cm_upstart % env.install_dir
    remote_file = '/etc/init/%s' % cloudman_boot_file
    put(cloudman_boot_file, remote_file, use_sudo=777)
    os.remove(cloudman_boot_file)

def _cleanup_ec2(env):
    """Clean up any extra files after building.
    """
    env.logger.info("Cleaning up for EC2 AMI creation")
    fnames = [".bash_history", "/var/log/firstboot.done", ".nx_setup_done",
              "/var/crash/*", "%s/ec2autorun.py.log" % env.install_dir]
    for fname in fnames:
        sudo("rm -f %s" % fname)
    # Stop Apache from starting automatically at boot (it conflicts with Galaxy's nginx)
    sudo('/usr/sbin/update-rc.d -f apache2 remove')

    # RabbitMQ fails to start if its database is embedded into the image
    # because it saves the current IP address or host name so delete it now.
    # When starting up, RabbitMQ will recreate that directory.
    with settings(warn_only=True):
        sudo('/etc/init.d/rabbitmq-server stop')
        sudo('stop rabbitmq-server')
        sudo('/etc/init.d/rabbitmq-server stop')
    sudo('initctl reload-configuration')
    for db_location in ['/var/lib/rabbitmq/mnesia', '/mnesia']:
        if exists(db_location):
            sudo('rm -rf %s' % db_location)
    # remove existing ssh host key pairs
    # http://docs.amazonwebservices.com/AWSEC2/latest/UserGuide/index.html?AESDG-chapter-sharingamis.htm
    sudo("rm -f /etc/ssh/ssh_host_*")
