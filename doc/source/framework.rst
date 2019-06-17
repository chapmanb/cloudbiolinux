
Using the CloudBioLinux Build Framework
---------------------------------------

-----------------------
Obtaining CloudBioLinux
-----------------------

CloudBioLinux can be obtained using `git <http://git-scm.com/>`_.

::

    % git clone https://github.com/chapmanb/cloudbiolinux.git
    % cd cloudbiolinux

-------------------------
Overview of the Framework
-------------------------


------------------------
Building Cloud Instances
------------------------

When building CloudBioLinux instances for the cloud, ``deploy/deploy.sh``
contains a script to automate cloud interactions and the installation of
CloudBioLinux. This script will require a system Python be available on your
system, but should otherwise install its own dependencies.

Before using the deployer you will need to create a settings file describing
your cloud credentials and connection information as well as any tweaks you
would like to make to the:

::

    % cd deploy
    % cp settings-sample-minimal.yaml settings.yaml

Before updating settings.yaml you will need to navigate the AWS management
console and obtain the following information.

* Your AWS Access ID and secret key (`access_id`, `secret_key`)
* Ubuntu EBS-backed AMI ID to target. This writeup was tested with 
  `ami-9b85eef2` (12.04.2 (64-bit) in us-east-1).
* Image size to use (e.g. m1-small)
* Availability zone (e.g. us-east-1)

Carefully scan through `settings.yaml` and change the properties marked as
requiring change. The word `UPDATE` in the comments indicates properties of
special interest that either don't have reasonable defaults or have reasonable
defaults but that I have deemed highly likely to be overridden. For this simple
example the only changes you will need to make are in the ``aws`` section.

Once you have updated ``settings.yaml``, launch and cloud instance and
configure it with the following command:

::

    % ./deploy.sh --action=install_biolinux --flavor=minimal

This command will configure CloudBioLinux with a minimal set of CloudBioLinux
packages. The set of packages that is installed is controlled by the
``--flavor`` command. More sophisticated setups that require using Amazon EBS
volumes and S3 buckets such as CloudMan clusters require additional
configuration as outlined below.

You can SSH into the newly created cloud instance with the command:

::

    % ./deploy.sh --action=ssh

-----------------------------------
Building CloudMan Enabled Instances
-----------------------------------

Before continuing, delete your previous instance and the file
``.vmlauncher_last_instance_aws``. TODO: Add action for this.

Building a more sophisticated CloudBioLinux image integrating tools such as
CloudMan requires additional settings. Please start by copying ``settings-
sample-cm.yaml`` to ``settings.yaml`` and repopulating the ``aws`` section
options. 

To fill out the remaining options found in this file, you will need to return
the AWS management console and do the following:

* You will need to setup a bucket to store your snaps file, here you will need the bucket name.
* You will need to setup two volumes in your target availability zone, one for
  Galaxy tools and data (perhaps 20Gb for testing) and one for galaxyIndices. Here you will need the volume ids.
* Generate a private a key (e.g. galaxy1.pem) and copy it into keys directory (or anywhere really), 
  also note the keypair_name corresponding to the key.

Next you will want to setup a directory to contain the S3 bucket contents that
will eventually be used by CloudMan to configure your cluster. Create a
directory (e.g. `/home/mary/marys_cloudman_bucket_contents`). Copy the files
from an existing CloudMan bucket here (e.g. http://s3.amazonaws.com /cloudman-
dev).

It is not really important how you download these files, but one quick option
is to use `s3cmd` tool:

::

    % sudo apt-get install s3cmd  # Or your OS's package manager
    % mkdir /home/mary/marys_cloudman_bucket_contents
    % s3cmd -r get s3://cloudman-dev /home/mary/marys_cloudman_bucket_contents

Here you can replace the CloudMan source (i.e. `cm.tar.gz`) or any of these
files to match the customized setup you would like. In particular you are
going to want to create a custom snaps.yaml file. Here is a simple outline
that we will fill out as we good.

::

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
and the availability zone you are targetting. We can update the snap_id's and
the default_mi after creating them.

Reopen ``settings.yaml`` and fill out the remaining properties, including the
volume ids you just created and the name of the bucket you used.

The following set of commands will now launch a new cloud instance, attach
and format tool and data volumes for CloudMan, build CloudBioLinux, snapshot
these volumes, and package the image.

::

    % ./deploy.sh --action=launch
    % ./deploy.sh --action=attach_volumes
    % ./deploy.sh --action=install_biolinux --flavor=cloudman/cloudman_and_galaxy
    % ./deploy.sh --action=snapshot_volumes
    % ./deploy.sh --action=detach_volumes
    % ./deploy.sh --action=package

Once a CloudMan AMI has been created, update `snaps.yaml` in your bucket
directory (e.g. `/home/mary/marys_cloudman_bucket_contents`) to reflect the
`snap_id`s and AMI created. These should all be available via the AWS
management console or by reviewing the output of the steps above.

Finally, you can upload your new bucket and launch a test CloudMan instance:

::

    % ./deploy.sh --action=sync_cloudman_bucket
    % ./deploy.sh --action=cloudman_launch

This last action (``cloudman_launch``) requires uncommenting the following
lines and updating the bucket name:

::

    #image_user_data:
    #  bucket_default: marys_cloudman_bucket

