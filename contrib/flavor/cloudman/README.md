# CloudMan Flavors

This document briefly describes the CloudMan/[Galaxy][10] flavors of
CloudBioLinux available and how to configure them. For quick, small
modifications to CloudMan your best bet is probably to modify an existing
instance, documentation on how to do this can be found [here][9].
Documentation for building a CloudMan image and corresponding volumes and
bucket from scratch using the CloudBioLinux deployer can be found [here][8].

## Core Flavors

This directory contains package lists used to setup [CloudMan][1] and/or CloudMan
and [Galaxy][2]. A minimal CloudMan machine can be configured with the following command:

         fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman

A slightly more flushed out instance can be installed via the command:

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:flavor=cloudman/cloudman

Finally, CloudMan and Galaxy can be installed together via the command:

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:flavor=cloudman/cloudman_and_galaxy

You can additionally configure your CloudMan and Galaxy instance by specifying
a configuration file: "-c <your_fabricrc.txt>" in the above command.

A subset of the parameters you may override via this mechanism includes (see
``config/fabricrc.txt`` for a full list):

* `galaxy_user` (default `galaxy`): User of Galaxy webapp
* `galaxy_home` (default `/mnt/galaxyTools/galaxy-central`): Galaxy installation directory
* `galaxy_tools_dir` (default `/mnt/galaxyTools/tools`): Galaxy tool's directory, formally called `install_dir` in mi-deployment's `tools_fabfile.py` file.
* `galaxy_loc_files` (default `/mnt/galaxyIndices`): Directory for installation of Galaxy loc files.
* `galaxy_repository` (default `https://bitbucket.org/galaxy/galaxy-central/`): This is the mercurial repository of Galaxy to install into.
* `galaxy_preconfigured_repository` (default `False`): If this is `True`, CloudBioLinux will not tweak the repository pulled in via mercurial for CloudMan. This is for applications which prepackage the needed changes to Galaxy directly into the Mercurial repository.
* `galaxy_cloud` (default `True`): If this is `True`, the Galaxy installation will be preconfigured for SGE and contain CloudMan branding.

## Galaxy Extensions

### Galaxy-P Flavors

There are two additional flavors assembled for the [Galaxy-P][4] project. 

Below are some recommended fabricrc.txt overrides for Galaxy-P.
  
* `dist_name = precise`
* `galaxy_repository = https://bitbucket.org/galaxyp/galaxyp-central/`
* `galaxy_tool_conf = /path/to/cloudbiolinux/contrib/flavor/proteomics/galaxyp/galaxy_tools.yaml`

Side Note: Galaxy-P is actively developed, tested, and deployed with Ubuntu 12.04
LTS, so this is what to target for best results. Feel free to let me
(jmchilton@gmail.com) know if there are issues or if you would like
help deploying in other environments and I will attempt to help in
whatever way I can.


#### The Galaxy-P Server Flavor

This flavor is used to build the internal and [public][5] Galaxy-P servers
hosted on the OpenStack cloud at the Minnesota Supercomputing Institute,
though can of course be used to build your own Galaxy-P environment. This
flavor is a sort of stripped down proteomics environment for Galaxy.

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:flavor=cloudman/cloudman_and_galaxyp

This flavor additionally requires a wine environment to be packaged using the
[proteomics-wine-env][7] project and made available to CloudBioLinux via the
`setup_proteomics_wine_env_script` fabricrc property.

#### The Galaxy-P Desktop Flavor

The Galaxy-P desktop flavor builds on the server flavor with additional
desktop tools and programming libraries designed to make it a broadly useful
environment for mass spec data analysis even outside the context of Galaxy
while also stripping out the wine related tools to avoid potential legal
concerns associated with redistributing such tools.

An Amazon AMI of this flavor will be available for launch on MSI's
[BioCloudCentral instance][6].

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:flavor=cloudman/cloudman_desktop_and_galaxyp


[1]: http://usecloudman.org/
[2]: http://usegalaxy.org/
[3]: http://cloudbiolinux.org/
[4]: http://getgalaxyp.org/
[5]: https://usegalaxyp.org/
[6]: https://biocloudcentral.msi.umn.edu/
[7]: https://github.com/jmchilton/proteomics-wine-env
[8]: https://github.com/chapmanb/cloudbiolinux/blob/master/deploy/cloudman.md
[9]: http://wiki.galaxyproject.org/CloudMan/CustomizeGalaxyCloud
[10]: http://galaxyproject.org
