from __future__ import print_function
import os
import time

from libcloud.compute.ssh import SSHClient
from libcloud.compute.base import NodeImage, NodeSize
from libcloud.compute.types import Provider
from libcloud.compute.providers import get_driver
import six

# Ubuntu 10.04 LTS (Lucid Lynx) Daily Build [20120302]
DEFAULT_AWS_IMAGE_ID = "ami-0bf6af4e"
DEFAULT_AWS_SIZE_ID = "m1.large"
DEFAULT_AWS_AVAILABILITY_ZONE = "us-west-1"

from fabric.api import local, env, sudo, put, run


class VmLauncher:

    def __init__(self, driver_options_key, options):
        self.driver_options_key = driver_options_key
        self.options = options
        self.__set_and_verify_key()

    def __set_and_verify_key(self):
        key_file = self.options.get('key_file', None)
        if not key_file:
            key_file = self._driver_options()['key_file']
        # Expand tildes in path
        self.key_file = os.path.expanduser(key_file)
        if not os.path.exists(self.key_file):
            raise Exception("Invalid or unspecified key_file option: %s" % self.key_file)

    def _get_driver_options(self, driver_option_keys):
        driver_options = {}
        for key in driver_option_keys:
            if key in self._driver_options():
                driver_options[key] = self._driver_options()[key]
        return driver_options

    def _driver_options(self):
        return self.options[self.driver_options_key]

    def get_key_file(self):
        return self.key_file

    def boot_and_connect(self):
        conn = self._connect_driver()
        node = self._boot()  # Subclasses should implement this, and return libcloud node like object
        self.conn = conn
        self.node = node
        self.uuid = node.uuid
        self.connect(conn)

    def _connect_driver(self):
        if not getattr(self, 'conn', None):
            self.conn = self._get_connection()
        return self.conn

    def _wait_for_node_info(self, f):
        initial_value = f(self.node)
        if initial_value:
            return self._parse_node_info(initial_value)
        while True:
            time.sleep(10)
            refreshed_node = self._find_node()
            refreshed_value = f(refreshed_node)
            if refreshed_value and not refreshed_value == []:
                return self._parse_node_info(refreshed_value)

    def _parse_node_info(self, value):
        if isinstance(value, six.string_types):
            return value
        else:
            return value[0]

    def _find_node(self):
        nodes = self.conn.list_nodes()
        node_uuid = self.node.uuid
        for node in nodes:
            if node.uuid == node_uuid:
                return node

    def destroy(self, node=None):
        self._connect_driver()
        if node == None:
            node = self.node
        self.conn.destroy_node(node)

    def __get_ssh_client(self):
        ip = self.get_ip()  # Subclasses should implement this
        key_file = self.get_key_file()
        ssh_client = SSHClient(hostname=ip,
                               port=self.get_ssh_port(),
                               username=self.get_user(),
                               key=key_file)
        return ssh_client

    def get_user(self):
        return "ubuntu"

    def get_ssh_port(self):
        return 22

    def connect(self, conn, tries=5):
        i = 0
        while i < tries:
            try:
                ssh_client = self.__get_ssh_client()
                conn._ssh_client_connect(ssh_client=ssh_client, timeout=60)
                return
            except:
                i = i + 1

    def list(self):
        self._connect_driver()
        return self.conn.list_nodes()

    def _boot(self):
        conn = self.conn
        boot_new = True
        last_instance_path = None
        if 'use_existing_instance' in self._driver_options():
            boot_new = False
            instance_id = self._driver_options()['use_existing_instance']
            if instance_id == "__auto__":
                last_instance_path = ".vmlauncher_last_instance_%s" % self.driver_options_key
                if not os.path.exists(last_instance_path):
                    boot_new = True
                else:
                    instance_id = open(last_instance_path, "r").read()
            if not boot_new:
                nodes = conn.list_nodes()
                nodes_with_id = [node for node in nodes if node.uuid == instance_id]
                if not nodes_with_id:
                    err_msg_template = "Specified use_existing_instance with instance id %s, but no such instance found."
                    raise Exception(err_msg_template % instance_id)
                node = nodes_with_id[0]
        if boot_new:
            node = self._boot_new(conn)
            if last_instance_path:
                open(last_instance_path, "w").write(node.uuid)
        return node

    def _image_from_id(self, image_id=None):
        image = NodeImage(id=image_id, name="", driver="")
        return image

    def _get_image_id(self, image_id=None):
        if not image_id:
            if 'image_id' in self._driver_options():
                image_id = self._driver_options()['image_id']
            else:
                image_id = self._default_image_id()
        return image_id

    def _default_image_id(self):
        return None

    def _get_default_size_id(self):
        return None

    def _get_size_id_option(self):
        return "size_id"

    def _size_from_id(self, size_id):
        size = NodeSize(id=size_id, name="", ram=None, disk=None, bandwidth=None, price=None, driver="")
        return size

    def _get_size_id(self, size_id=None):
        if not size_id:
            size_id_option = self._get_size_id_option()
            if size_id_option in self._driver_options():
                size_id = self._driver_options()[size_id_option]
            else:
                size_id = self._get_default_size_id()
        return size_id

    def _boot_new(self, conn):
        hostname = self.options.get("hostname", "vm_launcher_instance")
        node = self.create_node(hostname)
        return node

    def access_id(self):
        return self._driver_options()["access_id"]

    def secret_key(self):
        return self._driver_options()["secret_key"]

    def package_image_name(self):
        name = self._driver_options()["package_image_name"] or "cloudbiolinux"
        return name

    def package_image_description(self, default=""):
        description = self._driver_options().get("package_image_description", default)
        return description


class VagrantConnection:
    """'Fake' connection type to mimic libcloud's but for Vagrant"""

    def _ssh_client_connect(self, ssh_client):
        pass

    def destroy_node(self, node=None):
        local("vagrant halt")

    def list_nodes(self):
        return [VagrantNode()]


class VagrantNode:

    def __init__(self):
        self.name = "vagrant"
        self.uuid = "vagrant"


class VagrantVmLauncher(VmLauncher):
    """Launches vagrant VMs."""

    def _get_connection():
        return VagrantConnection()

    def __init__(self, driver_options_key, options):
        if not 'key_file' in options:
            options['key_file'] = os.path.join(os.environ["HOME"], ".vagrant.d", "insecure_private_key")
        VmLauncher.__init__(self, driver_options_key, options)
        self.uuid = "test"

    def _boot(self):
        local("vagrant up")
        return VagrantNode()

    def get_ip(self):
        return "33.33.33.11"

    def get_user(self):
        return "vagrant"

    def package(self, **kwds):
        local("vagrant package")


class OpenstackVmLauncher(VmLauncher):
    """ Wrapper around libcloud's openstack API. """

    def get_ip(self):
        return self._wait_for_node_info(lambda node: node.public_ips + node.private_ips)

    def _get_size_id_option(self):
        return "flavor_id"

    def create_node(self, hostname, image_id=None, size_id=None, **kwds):
        image_id = self._get_image_id()
        image = self._image_from_id(image_id)
        size_id = self._get_size_id()
        size = self._size_from_id(size_id)
        if 'ex_keyname' not in kwds:
            kwds['ex_keyname'] = self._driver_options()['ex_keyname']

        node = self.conn.create_node(name=hostname,
                                     image=image,
                                     size=size,
                                     **kwds)
        return node

    def _get_connection(self):
        driver = get_driver(Provider.OPENSTACK)
        openstack_username = self._driver_options()['username']
        openstack_api_key = self._driver_options()['password']

        driver_option_keys = ['host',
                              'secure',
                              'port',
                              'ex_force_auth_url',
                              'ex_force_auth_version',
                              'ex_force_base_url',
                              'ex_tenant_name']

        driver_options = self._get_driver_options(driver_option_keys)
        conn = driver(openstack_username,
                      openstack_api_key,
                      **driver_options)
        return conn

    def package(self, **kwds):
        name = kwds.get("name", self.package_image_name())
        self.conn.ex_save_image(self.node, name)

    def attach_public_ip(self, public_ip=None):
        if not public_ip:
            public_ip = self._driver_options()["public_ip"]
        self.conn._node_action(self.node, "addFloatingIp", address=public_ip)


class EucalyptusVmLauncher(VmLauncher):

    def get_ip(self):
        return self._wait_for_node_info(lambda node: node.public_ips)

    def _get_connection(self):
        driver = get_driver(Provider.EUCALYPTUS)
        driver_option_keys = ['secret',
                              'secure',
                              'port',
                              'host',
                              'path']

        driver_options = self._get_driver_options(driver_option_keys)
        ec2_access_id = self.access_id()
        conn = driver(ec2_access_id, **driver_options)
        return conn

    def create_node(self, hostname, image_id=None, size_id=None, **kwds):
        image_id = self._get_image_id()
        image = self._image_from_id(image_id)
        size_id = self._get_size_id()
        size = self._size_from_id(size_id)
        if 'ex_keyname' not in kwds:
            kwds['ex_keyname'] = self._driver_options()["keypair_name"]
        node = self.conn.create_node(name=hostname,
                                     image=image,
                                     size=size,
                                     **kwds)
        return node


class Ec2VmLauncher(VmLauncher):

    def get_ip(self):
        return self._wait_for_node_info(lambda node: node.extra['dns_name'])

    def boto_connection(self):
        """
        Establish a boto library connection (for functionality not available in libcloud).
        """
        import boto.ec2
        region = boto.ec2.get_region(self._availability_zone())
        ec2_access_id = self.access_id()
        ec2_secret_key = self.secret_key()
        return region.connect(aws_access_key_id=ec2_access_id, aws_secret_access_key=ec2_secret_key)

    def boto_s3_connection(self):
        from boto.s3.connection import S3Connection
        ec2_access_id = self.access_id()
        ec2_secret_key = self.secret_key()
        return S3Connection(ec2_access_id, ec2_secret_key)

    def _default_image_id(self):
        return DEFAULT_AWS_IMAGE_ID

    def package(self, **kwds):
        package_type = self._driver_options().get('package_type', 'default')
        if package_type == "create_image":
            self._create_image(**kwds)
        else:
            self._default_package(**kwds)

    def _create_image(self, **kwds):
        ec2_conn = self.boto_connection()
        instance_id = run("curl --silent http://169.254.169.254/latest/meta-data/instance-id")

        if "name" not in kwds:
            name = self.package_image_name()
        else:
            name = kwds["name"]

        if "description" not in kwds:
            description = self.package_image_description(default="")
        else:
            description = kwds["description"]

        image_id = ec2_conn.create_image(instance_id, name=name, description=description)
        if self._driver_options().get("make_public", False):
            ec2_conn.modify_image_attribute(image_id, attribute='launchPermission', operation='add', groups=['all'])

    def _default_package(self, **kwds):
        env.packaging_dir = "/mnt/packaging"
        sudo("mkdir -p %s" % env.packaging_dir)
        self._copy_keys()
        self._install_ec2_tools()
        self._install_packaging_scripts()

    def _install_ec2_tools(self):
        sudo("apt-add-repository ppa:awstools-dev/awstools")
        sudo("apt-get update")
        sudo('export DEBIAN_FRONTEND=noninteractive; sudo -E apt-get install ec2-api-tools ec2-ami-tools -y --force-yes')

    def _install_packaging_scripts(self):
        user_id = self._driver_options()["user_id"]
        bundle_cmd = "sudo ec2-bundle-vol -k %s/ec2_key -c%s/ec2_cert -u %s" % \
            (env.packaging_dir, env.packaging_dir, user_id)
        self._write_script("%s/bundle_image.sh" % env.packaging_dir, bundle_cmd)

        bucket = self._driver_options()["package_bucket"]
        upload_cmd = "sudo ec2-upload-bundle -b %s -m /tmp/image.manifest.xml -a %s -s %s" % \
            (bucket,  self.access_id(), self.secret_key())
        self._write_script("%s/upload_bundle.sh" % env.packaging_dir, upload_cmd)

        name = self.package_image_name()

        manifest = "image.manifest.xml"
        register_cmd = "sudo ec2-register -K %s/ec2_key -C %s/ec2_cert %s/%s -n %s" % (env.packaging_dir, env.packaging_dir, bucket, manifest, name)
        self._write_script("%s/register_bundle.sh" % env.packaging_dir, register_cmd)

    def _write_script(self, path, contents):
        full_contents = "#!/bin/bash\n%s" % contents
        sudo("echo '%s' > %s" % (full_contents, path))
        sudo("chmod +x %s" % path)

    def _copy_keys(self):
        ec2_key_path = self._driver_options()["x509_key"]
        ec2_cert_path = self._driver_options()["x509_cert"]
        put(ec2_key_path, "%s/ec2_key" % env.packaging_dir, use_sudo=True)
        put(ec2_cert_path, "%s/ec2_cert" % env.packaging_dir, use_sudo=True)

    def _availability_zone(self):
        if "availability_zone" in self._driver_options():
            availability_zone = self._driver_options()["availability_zone"]
        else:
            availability_zone = DEFAULT_AWS_AVAILABILITY_ZONE
        return availability_zone

    def _get_default_size_id(self):
        return DEFAULT_AWS_SIZE_ID

    def _get_location(self):
        availability_zone = self._availability_zone()
        locations = self.conn.list_locations()
        for location in locations:
            if location.availability_zone.name == availability_zone:
                break
        return location

    def create_node(self, hostname, image_id=None, size_id=None, location=None, **kwds):
        self._connect_driver()
        image_id = self._get_image_id(image_id)
        image = self._image_from_id(image_id)

        size_id = self._get_size_id(size_id)
        size = self._size_from_id(size_id)

        if not location:
            location = self._get_location()

        if not "ex_keyname" in kwds:
            keyname = self._driver_options()["keypair_name"]
            kwds["ex_keyname"] = keyname

        node = self.conn.create_node(name=hostname,
                                     image=image,
                                     size=size,
                                     location=location,
                                     **kwds)
        return node

    def attach_public_ip(self, public_ip=None):
        if not public_ip:
            public_ip = self._driver_options()["public_ip"]
        self.conn.ex_associate_addresses(self.node, public_ip)

    def _get_connection(self):
        driver = get_driver(Provider.EC2)
        ec2_access_id = self.access_id()
        ec2_secret_key = self.secret_key()
        conn = driver(ec2_access_id, ec2_secret_key)
        return conn


def build_vm_launcher(options):
    provider_option_key = 'vm_provider'
    # HACK to maintain backward compatibity on vm_host option
    if not 'vm_provider' in options and 'vm_host' in options:
        print("Using deprecated 'vm_host' setting, please change this to 'vm_provider'")
        provider_option_key = 'vm_host'
    driver = options.get(provider_option_key, 'aws')   # Will just fall back on EC2
    driver_options_key = driver
    if driver in options:
        # Allow multiple sections or providers per driver (e.g. aws-project-1).
        # Assume the driver is just the provider name unless the provider
        # section sets an explict driver option. In above example,
        # the aws-project-1 would have to have a "driver: 'aws'" option
        # set.
        provider_options = options.get(driver)
        driver = provider_options.get('driver', driver)
    driver_classes = {'openstack': OpenstackVmLauncher,
                      'vagrant': VagrantVmLauncher,
                      'eucalyptus': EucalyptusVmLauncher}
    driver_class = driver_classes.get(driver, Ec2VmLauncher)
    vm_launcher = driver_class(driver_options_key, options)
    return vm_launcher
