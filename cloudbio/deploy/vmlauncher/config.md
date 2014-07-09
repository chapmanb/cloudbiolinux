# Configuring Cloud Parameters

Currently four different virtual machine providers are implemented: `aws`
(default), `openstack`, `eucalyptus` (partial support), and `vagrant`. Request
for supporting additional cloud infrastructures can be created here
https://github.com/jmchilton/vm-launcher/issues/new or pull requests are
always welcome. The `vm-launcher` project is built heavily on Apache
[libcloud], so support should be implemented at that level first, though
dozens of cloud providers are currently implemented.

## aws 

This cloud supports the following compute parameters: `access_id`,
`secret_key`, `size_id`, `image_id`, `availability_zone`.

The aws driver supports two packaging modes, this is the code that is called
when the `package` action is executed. By default, package will cause some
scripts to be created on the remote server to aid in packaging, if however
`package_type` parameter is set to `create_image`, the Amazon EC2 CreateImage
(http://support.smartbear.com/viewarticle/22739/) operation will be used to
automatically package the target instance. `create_image` mode can only be
used for EBS backed instances, which is why the other more complex mode is the
default.

When `create_image` is enabled, the additional packaging parameters include
`package_image_name`, `package_image_description`, and `make_public`.

For the default packaging mode, many additional parameters related to S3 must
be set including `x509_cert`, `x509_key`, `user_id`, and `package_bucket`.

## openstack

OpenStack may be targetted using either the native OpenStack APIs or using the
EC2 compatibility layer (e.g. what boto does). To target the EC2 compatibility
use the `eucalyptus` driver, this `openstack` driver targets the native API.

This driver allows the following parameters `username`, `password`, `host`,
`secure` (boolean), `port`, `ex_force_uth_url`, `ex_force_base_url`,
`ex_force_auth_version`, `ex_tenant_name`, `flavor_id`, `image_id`,
`keypair_name`, and `package_image_name`.

## eucalyptus

Support for the eucalyptus driver is somewhat experimental at this time and automated packaging is not available. This driver can be configured via the following options:: `secret`, `secure`, `port`, `host`, `path`, `size_id`.

## vagrant

The vagrant driver supports no additional parameters, a precise64 box
should be configured though this can be tweaks by adjusting the file
`Vagrantfile`.

[libcloud]: http://libcloud.apache.org/



