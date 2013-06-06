# CloudBioLinux Deployer CloudMan QuickStart

As far as I can determine there is no current documentation on how to build
CloudMan instances from scratch. Thus I am collecting my unofficial notes on
how to do this here - spefically using the CloudBioLinux deployer.

You will need to navigate the AWS management console and obtain the following
information.

* Your AWS Access ID and secret key (`access_id`, `secret_key`)
* Ubuntu EBS-backed AMI ID to target. This writeup was tested with ami-9b85eef2 (12.04.2 (64-bit) in us-east-1)
* Image size to use (e.g. m1-small)
* Availibity zone (e.g. us-east-1)
* You will need to setup a bucket to store your snaps file, here you will need the bucket name.
* You will need to setup two volumes in your target availibity zone, one for
  Galaxy tools and data (perhaps 20Gb for testing) and one for galaxyIndices. Here you will need the volume ids.
* Generate a private a key (e.g. galaxy1.pem) and copy it into keys directory (or anywhere really), 
  also note the keypair_name corresponding to the key.

Create a directory (e.g. `/home/mary/marys_cloudman_bucket_contents`). Copy
the files from an existing CloudMan bucket here (e.g. http://s3.amazonaws.com
/cloudman-dev).

It is not really important how you download these files, but one quick option
is to use `s3cmd` tool:

    % sudo apt-get install s3cmd  # Or your OS's package manager
    % mkdir /home/mary/marys_cloudman_bucket_contents
    % s3cmd -r get s3://cloudman-dev /home/mary/marys_cloudman_bucket_contents

Here you can replace the CloudMan source (i.e. `cm.tar.gz`) or any of these
files to match the customized setup you would like. In particular you are
going to want to create a custom snaps.yaml file. Here is a simple outline
that we will fill out as we good.

    version: 1
    clouds:
      - name: amazon
        regions:
        - deployments:
          - name: GalaxyCloud
            filesystems:
            - name: galaxy
              roles: galaxyTools,galaxyData
              snap_id: snap-XXXXXXXXXXX
              mount_point: /mnt/galaxy
            - name: galaxyIndices
              roles: galaxyIndices
              snap_id: snap-XXXXXXXXXXXX
              mount_point: /mnt/galaxyIndices
            default_mi: ami-XXXXXXXXXXXXX
            bucket: marys_cloudman_bucket
          name: us-east-1

Immediately this template can be updated to reflect the bucket created above
and the availibity zone you are targetting. We can update teh snap_id's and
the default_mi after creating them.

Copy and modify `settings.yaml`:

    % git clone git://github.com/chapmanb/cloudbiolinux.git
    % cd cloudbiolinux/deploy
    % cp settings-sample-cm.yaml settings.yaml
    % vim settings.yaml # or your favorite editor

Carefully scan through `settings.yaml` and change the properties marked as requiring
change. The word `UPDATE` in the comments indicates properties of special
interest that either don't have reasonable defaults or have reasonable
defaults but that I have deemed highly likely to be overridden.

Now you can use the CloudBioLinux deployer to launch an image, attach volumes,
install biolinux, take needed snapshots, and package the whole thing up:

    % ./deploy.sh --action=launch
    % ./deploy.sh --action=attach_volumes
    % ./deploy.sh --action=install_biolinux --flavor=cloudman/cloudman_and_galaxy
    % ./deploy.sh --action=snapshot_volumes
    % ./deploy.sh --action=detach_volumes
    % ./deploy.sh --action=package

If at any point in the above process you need to interactively inspect the
state of the instance being configured you can do this via the following command:

    % ./deploy.sh --action=ssh

Once a CloudMan AMI has been created, update `snaps.yaml` in your bucket
directory (e.g. `/home/mary/marys_cloudman_bucket_contents`) to reflect the
`snap_id`s and AMI created. These should all be available via the AWS
management console or by reviewing the output of the steps above.

Finally, you can upload your new bucket and launch a test CloudMan instance:

    % ./deploy.sh --action=sync_cloudman_bucket
    % ./deploy.sh --action=cloudman_launch

## Customizing

The above example uses the `cloudman/cloudman_and_galaxy` CloudBioLinux
flavor, but there are additional flavors of CloudBioLinux available. Please
consult [this page][1] 
and choose the most appropriate flavor:

### Customizing Galaxy

Installing a customized Galaxy is as simple as overriding the
`galaxy_repository` variable in the `fabricrc_overrides` section of the
`settings.yaml`.

### Customizing Tools

Out of the box, CloudBioLinux can be configured to install dozens of
bioinformatic packages out of the box and adding additional packages is fairly
straight forward. One simply need to create a CloudBioLinux flavor that
configures which such packages are installed and specify that flavor (either
in the command-line as shown above or in `settings.yaml`).

Your custom flavor should include the `cloudman` packages. If your flavor
additionally includes `galaxy` (as the flavor `cloudman_and_galaxy` shown
above) packages and `install_tool_dependencies` is set to `True` in
`settings.yaml` - CloudBioLinux will setup a tool dependencies directory for
Galaxy. This allows multiple versions of an application to be installed in
isolation.

When enabled, the list of tools and versions that is installed can be found in
``cloudbiolinux/contrib/flavor/cloudman/tools.yaml <https://github.com/chapmanb/cloudbiolinux/blob/master/contrib/flavor/cloudman/tools.yaml>``. One can
modify that file directly or specify an entirely new file by setting the
``galaxy_tools_conf`` property in the `fabric_overrides` section of `settings.yaml`.

### Customizing CloudMan

CloudMan is downloaded from the bucket you specify and installed at system
startup. Hence one can simply place a customized version of CloudMan (tarred
up and named `cm.tar.gz`) in the bucket.

If `cloudman_repository`, `bucket_source`, and `bucket_default` are set in the
`cloudman` section of `settings.yaml`, then one can execute the following
command to quickly tar up the local copy of CloudMan (in
`cloudman_repository`) and update your target bucket.

    % ./deploy.sh --action=bundle_cloudman --action=sync_cloudman_bucket

[1]: https://github.com/chapmanb/cloudbiolinux/tree/master/contrib/flavor/cloudman