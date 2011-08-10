= Hacking BioLinux tips and tricks

The BioLinux tools allow building a full environment for Bioinformatics. The
design allows for flexible targets (Editions) and specializations (Flavors).

Please read the README and ./doc/vagrant documentation that come with the
BioLinux source tree first (see http://github.com/chapmanb/cloudbiolinux).

== Start with the Minimal edition

The Minimal edition is the smallest common denominator of all Editions, as it
installs the minimum of packages to bootstrap a full install. Once the
vagrant box is up and running, Minimal is invoked by

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/minimal/fabricrc_debian.txt install_biolinux:packagelist=$source/contrib/minimal/main.yaml

where $source points to your biolinux source tree. In fact, the testing script
in ./test/test_vagrant does exactly this! Try:

          cd $my_vms
          $source/test/test_vagrant --help

and the actual run:

          $source/test/test_vagrant

The main.yaml file ascertains the major editors are included, as well remote
access, version control, and the basic build system (gcc and friends). Note the
Minimal edition overwrites the (apt) sources file to make sure there are no
conflicts with user settings.

== Adding install packages

To expand on the package list you can define your own main.yaml, and pass that
in. In your main.yaml file add the meta-packages listed in
config/packages.yaml. Invoke your new package list with

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/minimal/fabricrc_debian.txt install_biolinux:packagelist=myproject/main.yaml

It is that simple!

If packages.yaml is not complete, you may suggest changing its contents in the
main repository. The alternative is to create your own flavor, which we will do
in a minute. The same strategy holds for the other definitions in the ./config
directory, such as for Ruby gems, Python eggs, Perl CPAN, R-CRAN etc.

== Define a Flavor

For a cross language Bio* project performance test I needed to create a special
version of BioLinux that would pull in a list of scripts and some additional
packages.  Starting from an existing edition (in this case the Minimum edition,
but it also works on top of BioNode and BioLinux editions), I created a new
flavor in ./contrib/flavor/pjotrp/biotest/biotestflavor.py, named BioTestFlavor
(note you also need an empty __init__.py file).  A Flavor class overrides the
Flavor methods defined in ./cloudbio/flavor/__init__.py, in particular
rewrite_config_items, a generic hook to rewrite a list of configured items (the
package lists), and post_install, a post installation hook. To see what packages
your Flavor wants to install, simply override rewrite_config_items, and add a print
statement. For example:

    def rewrite_config_items(self, name, packages):
        for package in packages:
          env.logger.info("Selected: "+name+" "+package)
        return packages

The flavor is itself is found through a fabricrc.txt file. The main package list may be
in a new main.yaml file.  Kicking it into submission:

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/flavor/pjotrp/biotest/fabricrc_debian.txt install_biolinux:packagelist=$source/contrib/flavor/pjotrp/biotest/main.yaml

The flavor module itsefl sets env.flavor on loading the module (this can only
happen once). For more examples see the files in ./contrib/flavor.

== Distribute a VirtualBox

With vagrant a box can be exported with

          vagrant package

== Flavor: change default sources (apt, yum, rpm)

Note: NYI

BioLinux creates a (default, or edition based) list of package sources. These
sources can be overridden by the Flavor.rewrite_apt_sources_list method - which
should return a new list.

== Flavor: install additional packages

The primary way of adding new packages is by creating a new main.yaml file, as
discussed above in ''Define a flavor''. In addition a flavor can define a
method: BioLinux creates a (default, or edition based) list of packages. These
sources can be overridden by the Flavor.rewrite_config_list method - which
should return a new list.

== Flavor: filter packages

To filter/remove packages from the default list, use rewrite_config_list to filter existing
meta packages.

== Flavor: rewrite Ruby gem, Perl CPAN, Python egg, R CRAN lists

The function rewrite_config_list also allows rewriting package lists for Ruby, Python, R,
Perl etc. The general idea is that the main
Editions define the inclusion of the main languages, and pull in Bio* related
packages. To override this behaviour use the rewrite functions, e.g.

    def rewrite_config_list(self, name, list):
      if name == 'ruby':
        return [ 'bio' ]
      return list

only allows the BioRuby 'bio' gem to be installed. This happens at the time your
meta main.yaml reads

    libraries:
      - ruby-libs

and pulls in ruby-libs.yaml. One ruby-libs.yaml is shared to make sure all
editions are up-to-date. Likewise for all the other yaml files. The
configuration options are at the main.yaml (meta-package) level, and by using
rewrite methods at Edition and Flavor levels.

== Flavor: install special software

BioLinux comes with a bag of tricks to install special software outside the
main package system. There are methods for checking out source repositories,
and building software. There are methods for accessing public data resources
(such as Amazon S3). These are so called custom installs which are defined in
custom.yaml. Each of these can be pulled in and are configured by code in the
./cloudbio/custom/ directory.  The method fetches names from custom.yaml that
delegate to a method in the custom/name.py program.  These mechanisms are
shared between BioLinux editions.

But, importantly, it is easy to role your own custom methods using a Flavor!
This mechanism can also be used to automatically run post-install software,
such as puppet, cfruby and chef.

For example, you can tell your flavor to clone a git repository, and execute
a script by adding a post_install method to your flavor. E.g.

            def post_install(self):
                env.logger.info("Starting post-install")
                if exists('Scalability'):
                    with cd('Scalability'):
                       run('git pull')
                else:
                   _fetch_and_unpack("git clone git://github.com/pjotrp/Scalability.git")
                # Now run a post installation routine
                run('./Scalability/scripts/hello.sh')

You can run only the post_install (convinient for testing!) using the post_install
target, e.g.

         fab -H hostname -f $source/fabfile.py -c  $flavor/fabricrc_debian.txt install_biolinux:packagelist=$flavor/main.yaml,target=post_install

For a full Flavor example see

    https://github.com/pjotrp/cloudbiolinux/blob/master/contrib/flavor/pjotrp/biotest/biotestflavor.py

= Tips and tricks

== Tip for checking BioLinux installation effects

To see the what a BioLinux install does to your system, store the settings of
the original (untouched) state of a VM:

1. Make a dump of the current installed package list

                 dpkg -l > dpkg-original-list.txt

2. Store the /etc tree - one way is to use git in /etc

After running BioLinux you can see what has been done to your system by diffing
against the package list, and checking /etc.

== Use the testing framework to create new Flavors

BioLinux comes with a testing framework in ./test. The frame work creates a new VM on
a local machine. You can add tests, to check if a VM is complete. See the main README file for
more information.
