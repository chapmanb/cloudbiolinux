This document is meant to layout work that is to be done and/or that
has been done in merging mi-deployment and galaxy-vm-launcher into
cloudbiolinux.

## mi-deployment

* data_fabfile.py - Seems all the work was already been done.

* ec2autorun.py - Ported with the same default behavior as of 8/12/12,
  but optional user data extensions to handle more galaxy-vm-launcher 
  style use cases.

* instance-to-ebs-ami.sh - No work done?

* mi-fabile.py - 

  * apt package installations - Ported over and update to date as of
    8/12/12.

  * setting up users - Ported over and update to date as of 8/12/12.

  * install nginx - Ported over at some point. TODO: Determine if the
    cloudbiolinux stuff is up-to-date. 
    
    Extension: The cloudbiolinux version has been extended to allow
    parameterization of the Galaxy webapp installation directory.

  * Configure SGE - Ported over. TODO: Determine if the cloudbiolinux
    stuff is up-to-date.

  * Install setuptools - This doesn't seem to be needed with
    cloudbiolinux, at least with Ubuntu 12.04.

  * Install s3fs - TODO: Port this functionality.
 
  * _configure_postgresql - Is this really needed? Have the same code in 
    galaxy-vm-launcher, but Cloudman seems to working without this cleanup. 

  * _install_proftpd - Ported over, up-to-date as of 8/12/12. Should move 
    configuration files to installed_files instead of mi-deployment.

  * _configure_ec2_autorun - Ported over and up-to-date as of 8/12/12.

  * _configure_sge - Ported over and up-to-date as of 8/12/12.

  * _configure_galaxy_env - TODO: Port this functionality.

  * _configure_nfs - Ported over and seems functionally equivalent as of 8/12/12.

  * _configure_logrotate - Files should be moved into cloudbiolinux/installed_files,
    but this is functionally equivalent as of 8/12/12.

  * _save_image_conf_support - TODO: Port this functionality.
  
  * _configure_xvfb - TODO: Port this functionality. 

  * _configure_bash - Ported over as part of cloudbio.cloudman's _setup_env method
    and up-to-date 

* nginx.conf - Ported over (8/12/12), extended to allow 
  parameterization of galaxy_home path.
* nginx_errdoc.tar.gz - Right now this is fetched from mi-deployment, 
  should be stored in cloudbiolinux.
* tools_fabfile.py - Ported over and up-to-date as of 8/12/12.
  Functionality split into installing Galaxy and installing actual tool dependencies.

  To install Galaxy the way tools_fabfile does, install cloudman via:
  fab -f fabfile.py -i <key> -H ubuntu@<IP> install_biolinux:packagelist=./contrib/cloudman/cloudman_and_galaxy.yaml  
  instead of :
  fab -f fabfile.py -i <key> -H ubuntu@<IP> install_custom:cloudman
  
  Several extension points on top of mi-deployment have been added
  (configured via fabricrc options), but the default behavior should
  be the same. These extensions include, allowing deployer to set:

  * galaxy_repository: Which repo to point at.
  * galaxy_preconfigured_repository: Override default of true to skip the 
    tweaking of the installed galaxy repository that tools_fabfile does. The
    galaxy-vm-launcher approach is to include these changesets 
    https://bitbucket.org/jmchilton/cloud-galaxy-dist
    in the repository Galaxy being installed from.
  * galaxy_conf_directory: Use this work 
    https://bitbucket.org/galaxy/galaxy-central/pull-request/44/allow-usage-of-directory-of-configuration
    to allow overridding specific galaxy universe_wsgi.ini properties at image configuration time.
    Any fabric environment properties of the form: galaxy_universe_XXXX=YYYY will end up at runtime
    as Galaxy universe_wsgi.ini properties of the form XXXX=YYYY inside [app:main]. My extensions to CloudMan also allow setting properties via this mechanism from user data at startup time. By default, cloudbiolinux properties will have a priority of 200 and cloudman properties a priority of 400, so if there are conflicts, CloudMan's startup time properties will override those of CloudBioLinux's. 
  
  The remaining functionality of installing Galaxy application and
  R/Bioconductor dependencies is available, but turned off by
  default. To enable these, update the fabric enviornment to set
  galaxy_install_dependencies=true and galaxy_install_r_packages. This
  functionality has been extended to dynamically read which packages
  and versions to install from the file: contrib/cloudman/tools.yaml.
  Multiple versions of the same tool can be installed via this
  mechanism. I have updated various software versions and tweaked the
  install methods to work under Ubuntu 12.04 for all tools.

* volume_manipulations_fab.py - No work has been done on this. This
  file has not been updated in a while, is it still useful?

* xvfb_default and xvfb_init - TODO: Port these files and mi-deployments' 
  configure_xvfb functionality over.

* conf_files/proftp.conf, conf_files/proftp.initd - Right now this is fetched from mi-  deployment, should be stored in cloudbiolinux.

* conf_files/vimrc - Right now this is fetched from mi-deployment, 
  should be stored in cloudbiolinux.

* conf_files/apps.yaml, conf_files/config.yaml - The cloudman and
  galaxy dependencies have been ported over and will be installed
  using cloudbiolinux's similar package management configuration
  mechanisms. The extensions (wrf, graphlab) have not been.

* conf_files/cloudman.logrotate  - File still be fetched from mi-deployment, 
  this should be read from installed_files instead.

* copy_snap/local_to_ebs_fab.py - No work has been done on this. This
  file has not been updated in a while, is it still useful?

* tools/ - No work on the custom tools has been done. Are there more
  cloudbiolinux-y ways to handle these?

## galaxy-vm-launcher

* lib/genomes.py - Uses cloudbiolinux alread. TODO: Just delete out of 
  galaxy-vm-launcher.

* lib/image.py - 

  Setting up users and installing packages is good to go.

  install_nginx in cloudbio/custom/cloudman should be refactored to a common shared place, but still invoked in install_cloudman. A seperate nginx conf file should be optional for galaxy-vm-launcher or we should find a good way to make the cloudman 
  paths optional (and enabled by default).

  _configure_postgresql - I don't think this is needed except maybe allowing 
  overridding of postgresql.conf etc file.

  _init_postgresql_data, _configure_nginx_service, start_nginx - This should all
  be ported over but disabled by default (for cloudman compatiablity).

  configure_xvfb - TODO: Configure this.

* lib/tools.py - Mostly merged already. galaxy-vm-launcher should be refactored 
  to use cloudbio/galaxy/tools.

* lib/galaxy.py - 

  * Setting up galaxy and options (functionality mostly available in cloudbiolinux now gvl needs to be refactored to use it.)

  * Setting up galaxy service, log, and database migrations. Functionality could be move to cloudbio/galaxy and disabled by default (not compatiable with cloudman).

  * Seeding galaxy with data. Functionality should remain in the galaxy-vm-launcher, 
  not really compatiable with cloudman. galaxy-vm-launcher should be refactored to
  use blend though and possibily some fork of blend with additional methods that 
  directly interact with the database the way galaxy-vm-launcher does.

* lib/deploy.py - File should largely remain as is.
  
  setup_taxonomy_data - Could be moved into cloudbio/biodata. Does cloudman use 
  the metagenomics tools? How is this currently being configured?

## Concrete TODO List:

### mi-deployment migration:

* *TODETERMINE*: Is install_nginx in cloudbiolinux up-to-date with mi-deployment?
* *TODETERMINE*: Is SGE configuration in cloudbiolinux up-to-date with mi-deployment?
* *TODETERMINE*: Is setuptools install in mi-deployment needed with cloudbiolinux? (Seems no. -John)
* *TODO*: Port mi-deployment s3fs install functionality to cloudbiolinux.
* *TODETERMINE*: Is _configure_postgresql functionality needed? Seems to work without it.
* *TODO*: Port mi-deployment _configure_galaxy_env functionality to cloudbiolinux.
* *TODO*: Port mi-deployment _save_image_conf_support functionality to cloudbiolinx
* *TODO*: Port mi-deployment _configure_xvfb functionality to cloudbiolinx
* *TODO*: Move required files for _configure_logrotate to cloudbiolinux installed_files and update setup procedure accordingly.
* *TODO*: Move required files for proftpd to cloudbiolinux installed_files and update setup procedure accordingly.
* *TODO*: Move nginx_errdoc.tar.gz to cloudbiolinux installed_files and update setup procedure accordingly.
* *TODO*: Move required files for vimrc to cloudbiolinux installed_files and update setup procedure accordingly.
* *TODO*: Port mi-deployment volume_manipulations_fab.py functionality to cloudbiolinux (if makes sense ).
* *TODO*: Port mi-deployment instance-to-ebs-ami.sh functionality to cloudbiolinux (if makes sense ).
* *TODO*: Port mi-deployment copy_snap/local_to_ebs_fab.py functionality to cloudbiolinux (if makes sense ).
* *TODO*: Port mi-deployment wrf, graphlab, and tools/* functionality to cloudbiolinux (if makes sense ).

### galaxy-vm-launcher migration:

* *TODO*: Refactor install_nginx in cloudbio.custom.cloudman so it can be used by 
galaxy-vm-launcher
* *TODO*: Determine and implement good way to make cloudman specific parts of nginx.conf optional.
* *TODO*: Move these procedures into cloudbiolinux - _init_postgresql_data, _configure_nginx_service, start_nginx (all disabled by default to ensure cloudman compat.).
* *TODO*: Refactor deploy.py to use cloudbio/galaxy/tools instead of gvl/lib/tools.py
* *TODO*: Refactor galaxy.py to install galaxy via cloudbiolinux methods.
* *TODO*: Move this functionality into cloudbiolinx - setting up galaxy init service, log, and database migrations (all disabled by default to ensure cloudman compat.)
* *TODETERMINE*: Could we move setup_taxonomy_data from gvl/lib/deploy.py into cloudbio/biodata somewhere? How is cloudman being configuring this data?
