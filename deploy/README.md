# CloudBioLinux Deployer

This CloudBioLinux deployer has grown out of the galaxy-vm-launcher and
can be used to launch cloud virtual machines, configure them with
Galaxy, and seed it with input data, genomes, workflows, etc.... More
recently, actions for installing CloudBioLinux and launching CloudMan
have been added.

## Prerequisites

The `deploy.sh` script should install the needed dependencies in a Python
virutal environment using venvburrito and doesn't require special permissions
as long as `python`, `easy_install`, and `git` are available.

## Specify settings

All deploy actions first require the existence of a setting file. 

    cp settings-sample-oldgalaxyvmlauncher.yaml settings.yaml

This file has numerous settings to customize how the deployer acts. Be
default, the deployer will target Amazon Web Services and `key_file`,
`access_id`, and `private_key` in the aws section of of this file must
be specified.

The argument `--settings=/path/to/custom_settings.yaml` may be passed
to `deploy.sh` to specify a custom path for this settings file.

## Configuring Galaxy

    ./deploy.sh --action=configure --action=transfer file1 file2 file3

When called this way, deploy.sh will launch a VM, configure Galaxy,
tools, and genomes. Once Galaxy is ready, it will transfer each of
provided input files to the newly launched VM, and use the Galaxy REST
API to add them to a Galaxy data library (and optionally a
history). Once all of that is complete, it will print a URL to the
screen telling the operator where to find the new Galaxy instance.

This does not install CloudMan, Galaxy is configured to run at startup
by an init script. A more traditional CloudMan workflow can be
achieved using the `install_biolinux` action described next.

## Installing CloudBioLinux

    ./deploy.sh --action=install_biolinux --action=package

This mode will launch an instance, install CloudBioLinux (a flavor can
be specified in setting.yaml), and `package` (see settings.yaml for
more details) the resulting virtual image.

## Additional Actions

The actions show above can be combined in different manners, for
instance `configure` and `package` can be used to configure a Galaxy
instance and package so that later `transfer` can be used without
requiring a full configure. Alternatively, `install_biolinux` can be
followed up with `transfer` to install CloudBioLinux and start
analyzing data without requiring Galaxy (be sure to set `use_galaxy:
False` in settings.yaml in this case).

You can see running instances on target cloud with this command: 

    ./deploy.sh --action=list

You can also destroy all running instances with this command:
    
    ./deploy.sh --action=destroy

If an existing CloudBioLinux image bundled with CloudMan has been
created and its image id set as `image_id` in the `cloudman` section
of `setting.yaml`, then this image can be launched for testing with:

    ./deploy.sh --action=launch_cloudman

The full list of actions can be found in `cloudbio/deploy/__init__.py`
and includes:

* `list`
* `destroy`
* `transfer`
* `destroy`
* `transfer`,
* `purge_galaxy`
* `setup_galaxy`
* `purge_tools`
* `setup_tools`
* `purge_genomes`
* `setup_genomes`
* `setup_ssh_key`
* `package`
* `setup_image`
* `launch` - Dummy action justs launches instance
* `install_biolinux`
* `cloudman_launch`

Additional composite actions are shortcuts for multiple actions - these include:

* `configure` - `setup_image`, `setup_tools`, `setup_genomes`, `setup_ssh_key`
* `reinstall_galaxy` - `purge_galaxy` and `setup_galaxy`
* `reinstall_genomes` - `purge_genomes` and `setup_genomes`
* `reinstall_tools` - `purge_tools` and `setup_tools`

## Configuring Cloud Provider

Cloud interactions are managed via the [vm-launcher] project, full
information on configuring different cloud providers can be found
[here][vm-launcher-config]

In brief, there are few different options for where to create the
VMs. Amazon EC2 is the default target, but it can also target
Eucalyptus or OpenStack based clouds. The ruby package `vagrant` can
be used to target virtual instances on your own machine.

[vm-launcher-config]: https://github.com/jmchilton/vm-launcher/blob/master/config.md
