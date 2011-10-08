#!/usr/bin/env python
"""Automate the final steps of configuring the CloudBioLinux exome example.

- Updates configuration file with server details.
- Adds RabbitMQ user and virtual host with correct permissions.

Run the script with sudo or the root user.
"""
import os
import time
import shutil
import socket
import subprocess
import ConfigParser

import yaml

def main():
    config_dir = "/export/data/galaxy"
    amqp_config = os.path.join(config_dir, "universe_wsgi.ini")
    pp_config = os.path.join(config_dir, "post_process.yaml")
    wait_until_mounted(amqp_config)
    update_amqp_config(amqp_config, socket.getfqdn())
    amqp_user, amqp_pass = read_ampq_config(amqp_config)
    amqp_vhost = read_pp_config(pp_config)
    setup_rabbitmq(amqp_vhost, amqp_user, amqp_pass)

def setup_rabbitmq(vhost, user, passwd):
    """Add virtual host, user and password to RabbitMQ.
    """
    base_cl = ["rabbitmqctl"]
    subprocess.check_call(base_cl + ["add_user", user, passwd])
    subprocess.check_call(base_cl + ["add_vhost", vhost])
    subprocess.check_call(base_cl + ["set_permissions", "-p", vhost,
                                     user, '".*"', '".*"', '".*"'])

def read_pp_config(fname):
    """Read AMQP vhost from YAML configuration file.
    """
    with open(fname) as in_handle:
        config = yaml.load(in_handle)
    return config["distributed"]["rabbitmq_vhost"]

def read_ampq_config(fname):
    """Get AMQP username and password from configuration file
    """
    config = ConfigParser.ConfigParser()
    config.read(fname)
    return (config.get("galaxy_amqp", "userid"),
            config.get("galaxy_amqp", "password"))

def update_amqp_config(fname, hostname):
    """Update AMQP configuration with internal hostname.
    """
    orig_stat = os.stat(fname)
    backup_file = "{0}.bak".format(fname)
    shutil.move(fname, backup_file)
    with open(backup_file) as in_handle:
        with open(fname, "w") as out_handle:
            for line in in_handle:
                if line.startswith("host ="):
                    line = "host = {0}\n".format(hostname)
                out_handle.write(line)
    # make updated file readable by initial user
    os.chown(fname, orig_stat.st_uid, orig_stat.st_gid)

def wait_until_mounted(fname):
    """Wait up to 3 minutes for mounted directory with file.
    """
    max_tries = 36
    wait_sec = 5
    num_tries = 0
    while 1:
        if os.path.exists(fname):
            break
        elif num_tries > max_tries:
            raise ValueError("Did not find {f} after {s} seconds.".format(
                f=fname, s=max_tries*wait_sec))
        time.sleep(wait_sec)
        num_tries += 1
 
if __name__ == "__main__":
    main()
