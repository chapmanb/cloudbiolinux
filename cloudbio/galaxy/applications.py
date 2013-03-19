"""
This file is largely derived from a similar file in mi-deployment written Dr.
Enis Afgan.

https://bitbucket.org/afgane/mi-deployment/src/8cba95baf98f/tools_fabfile.py

Long term it will be best to install these packages for Galaxy via the Tool
Shed, however many of these tools are not yet in the tool shed and the tool
shed installation is not currently available via the Galaxy API. Until such a
time as that is available, Galaxy dependencies may be installed via these
functions.

I have taken a first crack at harmonizing this with the rest of CloudBioLinux.
Wasn't able to reuse fastx_toolkit, tophat, cufflinks.

"""
import os

from fabric.api import sudo, run, cd
from fabric.contrib.files import exists

from cloudbio.custom.shared import _make_tmp_dir, _if_not_installed, _set_default_config
from cloudbio.custom.shared import _get_install, _configure_make, _make_copy, _fetch_and_unpack, _get_bin_dir


@_if_not_installed(None)
def install_fastx_toolkit(env):
    version = env.tool_version
    gtext_version = "0.6.1"
    url_base = "http://hannonlab.cshl.edu/fastx_toolkit/"
    fastx_url = "%sfastx_toolkit-%s.tar.bz2" % (url_base, version)
    gtext_url = "%slibgtextutils-%s.tar.bz2" % (url_base, gtext_version)
    pkg_name = 'fastx_toolkit'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % gtext_url)
            run("tar -xjvpf %s" % (os.path.split(gtext_url)[-1]))
            install_cmd = sudo if env.use_sudo else run
            with cd("libgtextutils-%s" % gtext_version):
                run("./configure --prefix=%s" % (install_dir))
                run("make")
                install_cmd("make install")
            run("wget %s" % fastx_url)
            run("tar -xjvpf %s" % os.path.split(fastx_url)[-1])
            with cd("fastx_toolkit-%s" % version):
                run("export PKG_CONFIG_PATH=%s/lib/pkgconfig; ./configure --prefix=%s" % (install_dir, install_dir))
                run("make")
                install_cmd("make install")


## TODO: Rework to use more of custom enhancements
@_if_not_installed("maq")
def install_maq(env):
    version = env["tool_version"]
    url = "http://downloads.sourceforge.net/project/maq/maq/%s/maq-%s.tar.bz2" \
            % (version, version)
    _get_install(url, env, _configure_make)


@_if_not_installed("macs14")
def install_macs(env):
    from cloudbio.custom.bio_nextgen  import install_macs as cbl_install_macs
    install_dir = env.system_install
    cbl_install_macs(env)
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("echo 'PYTHONPATH=%s/lib/python%s/site-packages:$PYTHONPATH' >> %s/env.sh" % (env.python_version, install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("tophat")
def install_tophat(env):
    version = env.tool_version
    url = 'http://tophat.cbcb.umd.edu/downloads/tophat-%s.Linux_x86_64.tar.gz' % version
    pkg_name = "tophat"
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
                install_cmd("mv * %s" % install_dir)


@_if_not_installed("cufflinks")
def install_cufflinks(env):
    version = env.tool_version
    url = 'http://cufflinks.cbcb.umd.edu/downloads/cufflinks-%s.Linux_x86_64.tar.gz' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
                install_cmd("mv * %s" % install_dir)


@_if_not_installed("megablast")
def install_megablast(env):
    version = env.tool_version
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/%s/blast-%s-x64-linux.tar.gz' % (version, version)
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd('blast-%s/bin' % version):
                    install_cmd("mv * %s" % install_dir)


@_if_not_installed("blastn")
def install_blast(env):
    version = env.tool_version
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/%s/ncbi-blast-%s-x64-linux.tar.gz' % (version[:-1], version)
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd('ncbi-blast-%s/bin' % version):
                bin_dir = _get_bin_dir(env)
                install_cmd("mv * '%s'" % bin_dir)


@_if_not_installed("sputnik")
def install_sputnik(env):
    version = env.tool_version
    url = 'http://bitbucket.org/natefoo/sputnik-mononucleotide/downloads/sputnik_%s_linux2.6_x86_64' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O sputnik %s" % url)
            install_cmd("mv sputnik %s" % install_dir)


@_if_not_installed("taxonomy2tree")
def install_taxonomy(env):
    version = env.tool_version
    url = 'http://bitbucket.org/natefoo/taxonomy/downloads/taxonomy_%s_linux2.6_x86_64.tar.gz' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
                install_cmd("mv * %s" % install_dir)


@_if_not_installed("add_scores")
def install_add_scores(env):
    version = env.tool_version
    url = 'http://bitbucket.org/natefoo/add_scores/downloads/add_scores_%s_linux2.6_x86_64' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O add_scores %s" % url)
            install_cmd("mv add_scores %s" % install_dir)


@_if_not_installed("HYPHY")
def install_hyphy(env):
    version = env.tool_version
    url = 'http://www.datam0nk3y.org/svn/hyphy'
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("svn co -r %s %s src" % (version, url))
            run("mkdir -p build/Source/Link")
            run("mkdir build/Source/SQLite")
            run("cp src/trunk/Core/*.{h,cp,cpp} build/Source")
            run("cp src/trunk/HeadlessLink/*.{h,cpp} build/Source/SQLite")
            run("cp src/trunk/NewerFunctionality/*.{h,cpp} build/Source/")
            run("cp src/SQLite/trunk/*.{c,h} build/Source/SQLite/")
            run("cp src/trunk/Scripts/*.sh build/")
            run("cp src/trunk/Mains/main-unix.cpp build/Source/main-unix.cxx")
            run("cp src/trunk/Mains/hyphyunixutils.cpp build/Source/hyphyunixutils.cpp")
            run("cp -R src/trunk/{ChartAddIns,DatapanelAddIns,GeneticCodes,Help,SubstitutionClasses,SubstitutionModels,TemplateBatchFiles,TopologyInference,TreeAddIns,UserAddins} build")
            run("rm build/Source/preferences.cpp")
            with cd("build"):
                run("bash build.sh SP")
            install_cmd("mv build/* %s" % install_dir)
    _update_default(env, install_dir)


@_if_not_installed(None)
def install_gatk(env):
    version = env.tool_version
    url = 'ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/GenomeAnalysisTK-%s.tar.bz2' % version
    pkg_name = 'gatk'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
        install_cmd("mkdir -p %s/bin" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O gatk.tar.bz2 %s" % url)
            run("tar -xjf gatk.tar.bz2")
            install_cmd("cp GenomeAnalysisTK-%s/*.jar %s/bin" % (version, install_dir))
    # Create shell script to wrap jar
    sudo("echo '#!/bin/sh' > %s/bin/gatk" % (install_dir))
    sudo("echo 'java -jar %s/bin/GenomeAnalysisTK.jar $@' >> %s/bin/gatk" % (install_dir, install_dir))
    sudo("chmod +x %s/bin/gatk" % install_dir)
    # env file
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)
    # Link jar to Galaxy's jar dir
    jar_dir = os.path.join(env.galaxy_jars_dir, pkg_name)
    if not exists(jar_dir):
        install_cmd("mkdir -p %s" % jar_dir)
    tool_dir = os.path.join(env.galaxy_tools_dir, pkg_name, 'default', 'bin')
    install_cmd('ln --force --symbolic %s/*.jar %s/.' % (tool_dir, jar_dir))
    install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, jar_dir))


@_if_not_installed("srma.jar")
def install_srma(env):
    version = env.tool_version
    mirror_info = "?use_mirror=voxel"
    url = 'http://downloads.sourceforge.net/project/srma/srma/%s/srma-%s.jar' \
            % (version[:3], version)
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            install_cmd("mv srma-%s.jar %s" % (version, install_dir))
            install_cmd("ln -f -s srma-%s.jar %s/srma.jar" % (version, install_dir))
    sudo("touch %s/env.sh" % install_dir)
    _update_default(env, install_dir)


@_if_not_installed("BEAM2")
def install_beam(env):
    url = 'http://www.stat.psu.edu/~yuzhang/software/beam2.tar'
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar xf %s" % (os.path.split(url)[-1]))
            install_cmd("mv BEAM2 %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("pass2")
def install_pass(env):
    url = 'http://www.stat.psu.edu/~yuzhang/software/pass2.tar'
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar xf %s" % (os.path.split(url)[-1]))
            install_cmd("mv pass2 %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("lps_tool")
def install_lps_tool(env):
    version = env.tool_version
    url = 'http://www.bx.psu.edu/miller_lab/dist/lps_tool.%s.tar.gz' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar zxf %s" % (os.path.split(url)[-1]))
            install_cmd("./lps_tool.%s/MCRInstaller.bin -P bean421.installLocation=\"%s/MCR\" -silent" % (version, install_dir))
            install_cmd("mv lps_tool.%s/lps_tool %s" % (version, install_dir))
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("echo 'MCRROOT=%s/MCR/v711; export MCRROOT' >> %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("plink")
def install_plink(env):
    version = env.tool_version
    url = 'http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-%s-x86_64.zip' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("unzip %s" % (os.path.split(url)[-1]))
            install_cmd("mv plink-%s-x86_64/plink %s" % (version, install_dir))
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed(None)
def install_fbat(env):
    version = env.tool_version
    url = 'http://www.biostat.harvard.edu/~fbat/software/fbat%s_linux64.tar.gz' % version.replace('.', '')
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar zxf %s" % (os.path.split(url)[-1]))
            install_cmd("mv fbat %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("Haploview_beta.jar")
def install_haploview(env):
    url = 'http://www.broadinstitute.org/ftp/pub/mpg/haploview/Haploview_beta.jar'
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            install_cmd("mv %s %s" % (os.path.split(url)[-1], install_dir))
            install_cmd("ln -s %s %s/haploview.jar" % (os.path.split(url)[-1], install_dir))
    _update_default(env, install_dir)


@_if_not_installed("eigenstrat")
def install_eigenstrat(env):
    version = env.tool_version
    url = 'http://www.hsph.harvard.edu/faculty/alkes-price/files/EIG%s.tar.gz' % version
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar zxf %s" % (os.path.split(url)[-1]))
            install_cmd("mv bin %s" % install_dir)
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)



@_if_not_installed("augustus")
def install_augustus(env):
    version = env.tool_version
    url = "http://bioinf.uni-greifswald.de/augustus/binaries/augustus.%s.tar.gz" % version
    install_dir = env.system_install
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _fetch_and_unpack(url, need_dir=False)
            env.safe_sudo("mkdir -p '%s'" % install_dir)
            env.safe_sudo("mv * '%s'" % install_dir)
    env.safe_sudo("echo 'PATH=%s/bin:%s/scripts:$PATH' > %s/env.sh" % (install_dir, install_dir, install_dir))
    env.safe_sudo("echo 'export AUGUSTUS_CONFIG_PATH=%s/config' >> %s/env.sh" % (install_dir, install_dir))
    _update_default(env, install_dir)


@_if_not_installed("SortSam.jar")
def install_picard(env):
    version = env.tool_version
    mirror_info = "?use_mirror=voxel"
    url = 'http://downloads.sourceforge.net/project/picard/picard-tools/%s/picard-tools-%s.zip' % (version, version)
    pkg_name = 'picard'
    install_dir = env.system_install
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("unzip %s" % (os.path.split(url)[-1]))
            install_cmd("mv picard-tools-%s/*.jar %s" % (version, install_dir))
    _update_default(env, install_dir)
    # set up the jars directory
    jar_dir = os.path.join(env.galaxy_jars_dir, 'picard')
    if not exists(jar_dir):
        install_cmd("mkdir -p %s" % jar_dir)
    tool_dir = os.path.join(env.galaxy_tools_dir, pkg_name, 'default')
    install_cmd('ln --force --symbolic %s/*.jar %s/.' % (tool_dir, jar_dir))
    install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, jar_dir))


@_if_not_installed("fastqc")
def install_fastqc(env):
    """ This tool is installed in Galaxy's jars dir """
    version = env.tool_version
    url = 'http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v%s.zip' % version
    pkg_name = 'FastQC'
    install_dir = os.path.join(env.galaxy_jars_dir)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with cd(install_dir):
        install_cmd("wget %s -O %s" % (url, os.path.split(url)[-1]))
        install_cmd("unzip -u %s" % (os.path.split(url)[-1]))
        install_cmd("rm %s" % (os.path.split(url)[-1]))
        with cd(pkg_name):
            install_cmd('chmod 755 fastqc')
        install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, pkg_name))


def _update_default(env, install_dir):
    sudo("touch %s/env.sh" % install_dir)
    sudo("chmod +x %s/env.sh" % install_dir)
    _set_default_config(env, install_dir)

#@if_tool_not_found()
#def install_emboss(env):
#    version = env.tool_version
#    url = 'ftp://emboss.open-bio.org/pub/EMBOSS/old/%s/EMBOSS-%s.tar.gz' % (version, version)
#    pkg_name = 'emboss'
#    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
#    install_cmd = sudo if env.use_sudo else run
#    if not exists(install_dir):
#        install_cmd("mkdir -p %s" % install_dir)
#    with _make_tmp_dir() as work_dir:
#        with cd(work_dir):
#            run("wget %s" % url)
#            run("tar -xvzf %s" % os.path.split(url)[-1])
#            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
#                run("./configure --prefix=%s" % install_dir)
#                run("make")
#                install_cmd("make install")
#    phylip_version = '3.6b'
#    url = 'ftp://emboss.open-bio.org/pub/EMBOSS/old/%s/PHYLIP-%s.tar.gz' % (version, phylip_version)
#    with _make_tmp_dir() as work_dir:
#        with cd(work_dir):
#            run("wget %s" % url)
#            run("tar -xvzf %s" % os.path.split(url)[-1])
#            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
#                run("./configure --prefix=%s" % install_dir)
#                run("make")
#                install_cmd("make install")

