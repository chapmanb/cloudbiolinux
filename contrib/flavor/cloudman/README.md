This directory contains package lists used to setup [CloudMan][1] and/or CloudMan
and [Galaxy][2]. A minimal CloudMan machine can be configured with the following command:

         fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman

A slightly more flushed out instance can be installed via the command:

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:flavor=cloudman/cloudman

Finally, CloudMan and Galaxy can be installed together via the command:

        fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:falvor=cloudman/cloudman_and_galaxy

You can additionally configure your CloudMan and Galaxy instance by specifing
a configuration file: "-c <your_fabricrc.txt>" in the above command.

A subset of the parameters you may override via this mechanism includes:

* `galaxy_user` (default `galaxy`): User of Galaxy webapp
* `galaxy_home` (default `/mnt/galaxyTools/galaxy-central`): Galaxy installation directory
* `galaxy_tools_dir` (default `/mnt/galaxyTools/tools`): Galaxy tool's directory, formally called `install_dir` in mi-deployment's `tools_fabfile.py` file.
* `galaxy_loc_files` (default `/mnt/galaxyIndices`): Directory for installation of Galaxy loc files.
* `galaxy_repository` (default `https://bitbucket.org/galaxy/galaxy-central/`): This is the mercurial repository of Galaxy to install into.
* `galaxy_preconfigured_repository` (default `False`): If this is `True`, CloudBioLinux will not tweak the repository pulled in via mercurial for CloudMan. This is for applications which prepackage the needed changes to Galaxy directly into the Mercurial repository.
* `galaxy_cloud` (default `True`): If this is `True`, the Galaxy installation will be preconfigured for SGE and contain CloudMan branding.

[1]: http://usecloudman.org/
[2]: http://usegalaxy.org/
[3]: http://cloudbiolinux.com/
