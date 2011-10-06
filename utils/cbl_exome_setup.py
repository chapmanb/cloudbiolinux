#!/usr/bin/env python
"""Automate the final steps of configuring the CloudBioLinux exome example.

- Updates configuration file with server details.
- Adds RabbitMQ user and virtual host with correct permissions.
"""
import os
import shutil
import socket
import subprocess
import ConfigParser

import yaml

def main():
    config_dir = "/export/data/galaxy"
    amqp_config = os.path.join(config_dir, "universe_wsgi.ini")
    pp_config = os.path.join(config_dir, "post_process.yaml")
    update_amqp_config(amqp_config, socket.getfqdn())
    amqp_user, amqp_pass = read_ampq_config(amqp_config)
    amqp_vhost = read_pp_config(pp_config)
    setup_rabbitmq(amqp_vhost, amqp_user, amqp_pass)

def setup_rabbitmq(vhost, user, passwd):
    """Add virtual host, user and password to RabbitMQ.
    """
    base_cl = ["sudo", "rabbitmqctl"]
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
    backup_file = "{0}.bak".format(fname)
    shutil.move(fname, backup_file)
    with open(backup_file) as in_handle:
        with open(fname, "w") as out_handle:
            for line in in_handle:
                if line.startswith("host ="):
                    line = "host = {0}\n".format(hostname)
                out_handle.write(line)
 
if __name__ == "__main__":
    main()
