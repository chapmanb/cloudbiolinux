"""Install proteomics tools not currently packaged.
"""

import os
import re

from fabric.api import cd
from fabric.context_managers import prefix

from shared import (_if_not_installed, _make_tmp_dir,
                    _get_install, _make_copy,
                    _java_install, _symlinked_java_version_dir,
                    _get_bin_dir, _get_install_subdir,
                    _fetch_and_unpack,
                    _create_python_virtualenv,
                    _get_bitbucket_download_url,
                    _write_to_file)
from cloudbio.galaxy.utils import _chown_galaxy

# Tools from Tabb lab are only available via TeamCity builds that
# and the artifacts eventually are deleted (I think), storing versions
# for CloudBioLinux at getgalaxyp.msi.umn.edu for safe keeping.
PROTEOMICS_APP_ARCHIVE_URL = "http://getgalaxyp.msi.umn.edu/downloads"


# TODO: Define TPP install root
@_if_not_installed("xinteract")
def install_transproteomic_pipeline(env):
    """
    """
    ## version should be of form X.X.X-codename
    default_version = "4.6.1-occupy"
    version = env.get("tool_version", default_version)
    version_parts = re.match("(\d\.\d)\.(\d)-(.*)", version)
    major_version = version_parts.group(1)
    revision = version_parts.group(2)
    codename = version_parts.group(3)
    if revision == "0":
        download_rev = ""
    else:
        download_rev = ".%s" % revision
    download_version = ("%s%s" % (major_version, download_rev))
    url_pieces = (major_version, codename, revision, download_version)
    url = 'http://sourceforge.net/projects/sashimi/files/Trans-Proteomic Pipeline (TPP)/TPP v%s (%s) rev %s/TPP-%s.tgz' % url_pieces

    def _chdir_src(work_cmd):
        def do_work(env):
            src_dir = "trans_proteomic_pipeline/src" if version == "4.6.1-occupy" else "src"
            with cd(src_dir):
                env.safe_append("Makefile.config.incl", "TPP_ROOT=%s/" % env["system_install"])
                env.safe_append("Makefile.config.incl", "TPP_WEB=/tpp/")
                env.safe_append("Makefile.config.incl", "XSLT_PROC=/usr/bin/xsltproc")
                env.safe_append("Makefile.config.incl", "CGI_USERS_DIR=${TPP_ROOT}cgi-bin")
                work_cmd(env)
        return do_work

    def _make(env):
        env.safe_run("make")
        env.safe_sudo("make install")
    _get_install(url, env, _chdir_src(_make))


@_if_not_installed("omssacl")
def install_omssa(env):
    default_version = "2.1.9"
    version = env.get("tool_version", default_version)
    url = 'ftp://ftp.ncbi.nih.gov/pub/lewisg/omssa/%s/omssa-%s.linux.tar.gz' % (version, version)
    env.safe_sudo("mkdir -p '%s'" % env["system_install"])
    ## OMSSA really wants mods.xml, usermods.xml, etc... in the same directory
    ## so just copying everything there.
    _get_install(url, env, _make_copy(find_cmd="ls --color=never -1", do_make=False))


@_if_not_installed("OpenMSInfo")
def install_openms(env):
    """
    See comments above, working on getting this to compile from source. In
    the meantime installing from deb will have to do.
    """
    default_version = "1.10.0"
    version = env.get("tool_version", default_version)
    dot_version = version[0:version.rindex('.')]
    url = 'http://downloads.sourceforge.net/project/open-ms/OpenMS/OpenMS-%s/OpenMS-%s.tar.gz' % (dot_version, version)

    def _make(env):
        with cd("contrib"):
            env.safe_run("cmake -DINSTALL_PREFIX=%s ." % env.get('system_install'))
            env.safe_run("make")
        env.safe_run("cmake -DINSTALL_PREFIX=%s ." % env.get('system_install'))
        env.safe_run("make")
        env.safe_sudo("make install")
    _get_install(url, env, _make)


@_if_not_installed("LTQ-iQuant")
def install_tint_proteomics_scripts(env):
    default_version = "1.19.19"
    version = env.get("tool_version", default_version)
    url = "http://artifactory.msi.umn.edu/simple/ext-release-local/msi/umn/edu/tint-proteomics-scripts/%s/tint-proteomics-scripts-%s.zip" % (version, version)

    def install_fn(env, install_dir):
        env.safe_sudo("mv * '%s'" % install_dir)
        bin_dir = _get_bin_dir(env)
        for script in ["ITraqScanSummarizer", "LTQ-iQuant", "LTQ-iQuant-cli", "MgfFormatter"]:
            env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, script), bin_dir))
        env.safe_sudo("chmod +x '%s'/*" % bin_dir)

    _java_install("tint-proteomics-scripts", version, url, env, install_fn)


@_if_not_installed("ms2preproc")
def install_ms2preproc(env):
    default_version = "2009"
    version = env.get("tool_version", default_version)
    get_cmd = 'wget "http://software.steenlab.org/ms2preproc/ms2preproc.zip" -O ms2preproc.zip'

    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            env.safe_run(get_cmd)
            env.safe_run("unzip ms2preproc.zip")
            with cd("ms2preproc"):
                env.safe_run("mv ms2preproc-r2821-x86_64 ms2preproc-x86_64")
                env.safe_run("chmod +x ms2preproc-x86_64")
                install_dir = _get_bin_dir(env)
                env.safe_sudo("mv ms2preproc-x86_64 '%s'/ms2preproc" % install_dir)


@_if_not_installed("MZmine")
def install_mzmine(env):
    default_version = "2.10"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/mzmine/mzmine2/%s/MZmine-%s.zip" % (version, version)

    def install_fn(env, install_dir):
        ## Enhanced MZmine startup script that works when used a symbolic link and tailored for CloudBioLinux.
        _get_gist_script(env, "https://gist.github.com/jmchilton/5474421/raw/15f3b817fa82d5f5e2143ee08bd248efee951d6a/MZmine")
        # Hack for multi-user environment.
        env.safe_sudo("chmod -R o+w conf")
        env.safe_sudo("mv * '%s'" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "MZmine"), os.path.join(bin_dir, "MZmine")))

    _java_install("mzmine2", version, url, env, install_fn)


@_if_not_installed("SearchGUI")
def install_searchgui(env):
    default_version = "1.13.1"
    version = env.get("tool_version", default_version)
    url = "http://searchgui.googlecode.com/files/SearchGUI-%s_mac_and_linux.zip" % version

    def install_fn(env, install_dir):
        dir_name = "SearchGUI-%s_mac_and_linux" % version
        env.safe_sudo("tar -xf %s.tar" % dir_name)
        with cd(dir_name):
            _get_gist_script(env, "https://gist.github.com/jmchilton/5002161/raw/dc9fa36dd0e6eddcdf43cd2b659e4ecee5ad29df/SearchGUI")
            _get_gist_script(env, "https://gist.github.com/jmchilton/5002161/raw/b97fb4d9fe9927de1cfc5433dd1702252e9c0348/SearchCLI")
            # Fix known bug with SearchGUI version 1.12.2
            env.safe_sudo("find -iname \"*.exe\" -exec rename s/.exe// {} \;")
            # Hack for multi-user environment.
            env.safe_sudo("chmod -R o+w resources")
            env.safe_sudo("mv * '%s'" % install_dir)
            bin_dir = os.path.join(env.get("system_install"), "bin")
            env.safe_sudo("mkdir -p '%s'" % bin_dir)
            env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "SearchGUI"), os.path.join(bin_dir, "SearchGUI")))
            env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "SearchCLI"), os.path.join(bin_dir, "SearchCLI")))

    _unzip_install("SearchGUI", version, url, env, install_fn)


@_if_not_installed("psm_eval")
def install_psm_eval(env):
    default_version = "0.1.0"
    version = env.get("tool_version", default_version)
    url = "git clone https://github.com/jmchilton/psm-eval.git"

    def install_fn(env, install_dir):
        env.safe_sudo("cp -r psm-eval/* '%s'" % install_dir)
        _create_python_virtualenv(env, "psme", "%s/requirements.txt" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "psm_eval"), os.path.join(bin_dir, "psm_eval")))

    _unzip_install("psm_eval", version, url, env, install_fn)


@_if_not_installed("PeptideShaker")
def install_peptide_shaker(env):
    default_version = "0.20.1"
    version = env.get("tool_version", default_version)
    url = "http://peptide-shaker.googlecode.com/files/PeptideShaker-%s.zip" % version

    def install_fn(env, install_dir):
        _get_gist_script(env, "https://gist.github.com/jmchilton/5002161/raw/f1fe76d6e6eed99a768ed0b9f41c2d0a6a4b24b7/PeptideShaker")
        _get_gist_script(env, "https://gist.github.com/jmchilton/5002161/raw/8a17d5fb589984365284e55a98a455c2b47da54f/PeptideShakerCLI")
        # Hack for multi-user environment.
        env.safe_sudo("chmod -R o+w resources")
        env.safe_sudo("mv * '%s'" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "PeptideShaker"), os.path.join(bin_dir, "PeptideShaker")))
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "PeptideShakerCLI"), os.path.join(bin_dir, "PeptideShakerCLI")))

    _java_install("PeptideShaker", version, url, env, install_fn)


def _get_gist_script(env, url):
    name = url.split("/")[-1]
    env.safe_sudo("wget '%s'" % url)
    env.safe_sudo("chmod +x '%s'" % name)


@_if_not_installed("Mayu")
def install_mayu(env):
    default_version = "1.06"
    version = env.get("tool_version", default_version)
    url = "http://proteomics.ethz.ch/muellelu/web/LukasReiter/Mayu/package/Mayu.zip"

    def install_fn(env, install_dir):
        share_dir = _get_install_subdir(env, "share")
        env.safe_sudo("mv Mayu '%s'" % share_dir)
        bin_dir = _get_bin_dir(env)
        executable = "%s/Mayu" % bin_dir
        env.safe_sudo("""echo '#!/bin/bash\ncd %s/Mayu; perl Mayu.pl \"$@\"' > %s """ % (share_dir, executable))
        env.safe_sudo("chmod +x '%s'" % executable)

    _unzip_install("mayu", version, url, env, install_fn)


def install_pride_inspector(env):
    default_version = "1.3.0"
    version = env.get("tool_version", default_version)
    url = "http://pride-toolsuite.googlecode.com/files/pride-inspector-%s.zip" % version

    def install_fn(env, install_dir):
        _get_gist_script(env, "https://gist.github.com/jmchilton/5474788/raw/6bcffd8680ec0e0301af44961184529a1f76dd3b/pride-inspector")
        # Hack for multi-user environment.
        env.safe_sudo("chmod -R o+w log config")
        env.safe_sudo("mv * '%s'" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "pride-inspector"), os.path.join(bin_dir, "pride-inspector")))

    _unzip_install("pride_inspector", version, url, env, install_fn, "PRIDE_Inspector")


def install_pride_converter2(env):
    default_version = "2.0.17"
    version = env.get("tool_version", default_version)
    url = "http://pride-converter-2.googlecode.com/files/pride-converter-%s-bin.zip" % version

    def install_fn(env, install_dir):
        _get_gist_script(env, "https://gist.github.com/jmchilton/5475119/raw/4e9135ada5114ba149f3ebc8965aee242bfc776f/pride-converter")
        # Hack for multi-user environment.
        env.safe_sudo("mkdir log; chmod o+w log")
        env.safe_sudo("mv * '%s'" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' %s" % (os.path.join(install_dir, "pride-converter"), os.path.join(bin_dir, "pride-converter")))

    _unzip_install("pride_converter2", version, url, env, install_fn, ".")


def _unzip_install(pname, version, url, env, install_fn, dir_name="."):
    install_dir = _symlinked_java_version_dir(pname, version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                _fetch_and_unpack(url, need_dir=False)
                with cd(dir_name):
                    install_fn(env, install_dir)


@_if_not_installed("SuperHirnv03")
def install_superhirn(env):
    default_version = "0.03"
    version = env.get("tool_version", default_version)
    url = "https://github.com/jmchilton/SuperHirn/zipball/%s/SuperHirn.zip" % version

    def _chdir(work_cmd):
        def do_work(env):
            with cd("SuperHirnv03/make"):
                work_cmd(env)
        return do_work

    _get_install(url, env, _chdir(_make_copy(find_cmd="find -perm -100 -name 'SuperHirn*'")))


@_if_not_installed("percolator")
def install_percolator(env):
    default_version = "2_04"
    version = env.get("tool_version", default_version)
    url = "https://github.com/downloads/percolator/percolator/percolator_%s_full_src.tar.gz" % version

    def make(env):
        with cd(".."):
            env.safe_run("env")
            env.safe_run("cmake -DCMAKE_INSTALL_PREFIX='%s' . " % env.system_install)
            env.safe_run("make -j8")
            env.safe_sudo("make install")

    _get_install(url, env, make)


@_if_not_installed("PepNovo")
def install_pepnovo(env):
    default_version = "20120423"
    version = env.get("tool_version", default_version)
    url = "http://proteomics.ucsd.edu/Downloads/PepNovo.%s.zip" % version

    def install_fn(env, install_dir):
        with cd("src"):
            env.safe_run("make")
            env.safe_sudo("mkdir -p '%s/bin'" % env.system_install)
            env.safe_sudo("mkdir -p '%s/share/pepnovo'" % env.system_install)
            env.safe_sudo("mv PepNovo_bin '%s/bin/PepNovo'" % env.system_install)
            env.safe_sudo("cp -r '../Models' '%s/share/pepnovo'" % env.system_install)

    _unzip_install("pepnovo", version, url, env, install_fn)


@_if_not_installed("crux")
def install_crux(env):
    default_version = "1.39"
    version = env.get("tool_version", default_version)
    url = "http://noble.gs.washington.edu/proj/crux/download/crux_%s-x86_64-Linux.zip" % version

    def _move(env):
        bin_dir = _get_bin_dir(env)
        env.safe_sudo("mv bin/* '%s'" % (bin_dir))

    _get_install(url, env, _move)


@_if_not_installed("Fido")
def install_fido(env):
    version = "2011"
    url = 'http://noble.gs.washington.edu/proj/fido/fido.tar.gz'

    # Adapted from Jorrit Boekel's mi-deployment fork
    # https://bitbucket.org/glormph/mi-deployment-protoeimcs
    def _chdir_src(work_cmd):
        def do_work(env):
            with cd("src/cpp"):
                env.safe_append('tmpmake', 'SHELL=/bin/bash')
                env.safe_append('tmpmake', 'prefix=%s' % env.get("system_install"))
                env.safe_append('tmpmake', 'CPPFLAGS=-Wall -ffast-math -march=x86-64 -pipe -O4 -g')
                env.safe_run('cat makefile |grep BINPATH -A 9999 >> tmpmake')
                env.safe_run('cp tmpmake makefile')
                work_cmd(env)
        return do_work

    _get_install(url, env, _chdir_src(_make_copy(find_cmd="find ../../bin -perm -100 -name 'Fido*'")))


def install_ipig(env):
    """ This tool is installed in Galaxy's jars dir """
    # This galaxy specific download probable doesn't belong in this file.
    default_version = "r5"
    version = env.get("tool_version", default_version)
    url = 'http://downloads.sourceforge.net/project/ipig/ipig_%s.zip' % version
    pkg_name = 'ipig'
    install_dir = os.path.join(env.galaxy_jars_dir, pkg_name)
    install_cmd = env.safe_sudo if env.use_sudo else env.safe_run
    install_cmd("mkdir -p %s" % install_dir)
    with cd(install_dir):
        install_cmd("wget %s -O %s" % (url, os.path.split(url)[-1]))
        install_cmd("unzip -u %s" % (os.path.split(url)[-1]))
        install_cmd("rm %s" % (os.path.split(url)[-1]))
        install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, install_dir))


def install_peptide_to_gff(env):
    default_version = "master"
    version = env.get("tool_version", default_version)
    repository = "hg clone https://jmchilton@bitbucket.org/galaxyp/peptide_to_gff"

    def install_fn(env, install_dir):
        env.safe_sudo("cp -r peptide_to_gff/* '%s'" % install_dir)
        _create_python_virtualenv(env, "peptide_to_gff", "%s/requirements.txt" % install_dir)
        bin_dir = os.path.join(env.get("system_install"), "bin")
        env.safe_sudo("mkdir -p '%s'" % bin_dir)
        env.safe_sudo("ln -s '%s' '%s'" % (os.path.join(install_dir, "peptide_to_gff"), os.path.join(bin_dir, "peptide_to_gff")))

    _unzip_install("peptide_to_gff", version, repository, env, install_fn)


def install_galaxy_protk(env):
    """This method installs Ira Cooke's ProtK framework. Very galaxy specific,
    can only be installed in context of custom Galaxy tool.


    By default this will install ProtK from rubygems server, but if
    env.protk_version is set to <version>@<url> (e.g.
    1.1.5@https://bitbucket.org/iracooke/protk-working) the
    gem will be cloned with hg and installed from source.
    """
    if not env.get('galaxy_tool_install', False):
        from cloudbio.custom.galaxy import _prep_galaxy
        _prep_galaxy(env)
    default_version = "1.2.2"
    version = env.get("tool_version", default_version)
    version_and_revision = version
    install_from_source = version_and_revision.find("@") > 0
    # e.g. protk_version = 1.1.5@https://bitbucket.org/iracooke/protk-working
    if install_from_source:
        (version, revision) = version_and_revision.split("@")
        url = _get_bitbucket_download_url(revision, "https://bitbucket.org/iracooke/protk")
    else:
        version = version_and_revision

    ruby_version = "1.9.3"
    force_rvm_install = False
    with prefix("HOME=~%s" % env.galaxy_user):
        def rvm_exec(env, cmd="", rvm_cmd="use", with_gemset=False):
            target = ruby_version if not with_gemset else "%s@%s" % (ruby_version, "protk-%s" % version)
            prefix = ". $HOME/.rvm/scripts/rvm; rvm %s %s; " % (rvm_cmd, target)
            env.safe_sudo("%s %s" % (prefix, cmd), user=env.galaxy_user)
        if not env.safe_exists("$HOME/.rvm") or force_rvm_install:
            env.safe_sudo("curl -L get.rvm.io | bash -s stable; source ~%s/.rvm/scripts/rvm" % (env.galaxy_user), user=env.galaxy_user)
            rvm_exec(env, rvm_cmd="install")
            rvm_exec(env, cmd="rvm gemset create protk-%s" % version)
        if not install_from_source:
            # Typical rubygem install
            rvm_exec(env, "gem install  --no-ri --no-rdoc protk -v %s" % version, with_gemset=True)
        else:
            with cd("~%s" % env.galaxy_user):
                env.safe_sudo("rm -rf protk_source; hg clone '%s' protk_source" % url, user=env.galaxy_user)
                rvm_exec(env, "cd protk_source; gem build protk.gemspec; gem install protk", with_gemset=True)

        protk_properties = {}
        ## ProtK can set these up itself, should make that an option.
        protk_properties["tpp_root"] = os.path.join(env.galaxy_tools_dir, "transproteomic_pipeline", "default")
        protk_properties['openms_root'] = "/usr"  # os.path.join(env.galaxy_tools_dir, "openms", "default", "bin")
        ### Assumes omssa, blast, and transproteomic_pipeline CBL galaxy installs.
        protk_properties['omssa_root'] = os.path.join(env.galaxy_tools_dir, "omssa", "default", "bin")
        protk_properties['blast_root'] = os.path.join(env.galaxy_tools_dir, "blast+", "default")
        protk_properties['pwiz_root'] = os.path.join(env.galaxy_tools_dir, "transproteomic_pipeline", "default", "bin")
        # Other properties: log_file, blast_root
        env.safe_sudo("mkdir -p \"$HOME/.protk\"", user=env.galaxy_user)
        env.safe_sudo("mkdir -p \"$HOME/.protk/Databases\"", user=env.galaxy_user)
        import  yaml
        _write_to_file(yaml.dump(protk_properties), "/home/%s/.protk/config.yml" % env.galaxy_user, "0755")

        rvm_exec(env, "protk_setup.rb galaxyenv", with_gemset=True)

        install_dir = os.path.join(env.galaxy_tools_dir, "galaxy_protk", version)
        env.safe_sudo("mkdir -p '%s'" % install_dir)
        _chown_galaxy(env, install_dir)
        env.safe_sudo('ln -s -f "$HOME/.protk/galaxy/env.sh" "%s/env.sh"' % install_dir, user=env.galaxy_user)
        with cd(install_dir):
            with cd(".."):
                env.safe_sudo("ln -s -f '%s' default" % version)


@_if_not_installed("myrimatch")
def install_myrimatch(env):
    default_version = "2.1.131"
    _install_tabb_tool(env, default_version, "myrimatch-bin-linux-x86_64-gcc41-release", ["myrimatch"])


@_if_not_installed("pepitome")
def install_pepitome(env):
    default_version = "1.0.45"
    _install_tabb_tool(env, default_version, "pepitome-bin-linux-x86_64-gcc41-release", ["pepitome"])


@_if_not_installed("directag")
def install_directag(env):
    default_version = "1.3.62"
    _install_tabb_tool(env, default_version, "directag-bin-linux-x86_64-gcc41-release", ["adjustScanRankerScoreByGroup", "directag"])


@_if_not_installed("tagrecon")
def install_tagrecon(env):
    default_version = "1.4.63"
    # TODO: Should consider a better way to handle the unimod xml and blosum matrix.
    _install_tabb_tool(env, default_version, "tagrecon-bin-linux-x86_64-gcc41-release", ["tagrecon", "unimod.xml", "blosum62.fas"])


@_if_not_installed("idpQonvert")
def install_idpqonvert(env):
    default_version = "3.0.475"
    version = env.get("tool_version", default_version)
    url = "%s/idpQonvert_%s" % (PROTEOMICS_APP_ARCHIVE_URL, version)
    env.safe_run("wget --no-check-certificate -O %s '%s'" % ("idpQonvert", url))
    env.safe_run("chmod 755 idpQonvert")
    env.safe_sudo("mkdir -p '%s/bin'" % env["system_install"])
    env.safe_sudo("mv %s '%s/bin'" % ("idpQonvert", env["system_install"]))
    env.safe_sudo("chmod +x '%s/bin/idpQonvert'" % env["system_install"])


def _install_tabb_tool(env, default_version, download_name, exec_names):
    version = env.get("tool_version", default_version)
    url = "%s/%s-%s.tar.bz2" \
        % (PROTEOMICS_APP_ARCHIVE_URL, download_name, version.replace(".", "_"))
    _fetch_and_unpack(url, False)
    env.safe_sudo("mkdir -p '%s/bin'" % env["system_install"])
    for exec_name in exec_names:
        env.safe_sudo("mv %s '%s/bin'" % (exec_name, env["system_install"]))
