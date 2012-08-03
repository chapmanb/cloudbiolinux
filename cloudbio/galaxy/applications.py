"""
This file is largely derived from a similar file in mi-deployment written Dr.
Enis Afgan.

https://bitbucket.org/afgane/mi-deployment/src/8cba95baf98f/tools_fabfile.py

Long term it will be best to install these packages for Galaxy via the Tool
Shed, however many of these tools are not yet in the tool shed and the tool
shed installation is not currently available via the Galaxy API. Until such a
time as that is available, Galaxy dependencies may be installed via the
install_tools function defined below.

There is obviously huge overlap with stuff found in this file and other files
in CloudBioLinux. Once the Galaxy Tool Shed has matured, if this file continues
to prove useful, this file should be harmonized with the rest of CloudBioLinux.
"""
import os
import time
import re

from fabric.api import sudo, run, cd
from fabric.contrib.files import exists, settings

from cloudbio.custom.shared import _make_tmp_dir


def if_tool_not_found():
    def argcatcher(func):
        def decorator(*args, **kwargs):
            env = args[0]
            version = args[1]
            package_name = func.__name__[len("install_"):]
            env.logger.debug("Checking the existence of %s" % package_name)
            env_file = os.path.join(env.galaxy_tools_dir, package_name, version)
            if not exists(env_file):
                return func(*args, **kwargs)
            else:
                env.logger.debug("Tool %s already appears to be installed (%s found), skipping." % (package_name, env_file))
        return decorator
    return argcatcher


@if_tool_not_found()
def install_ucsc_tools(env, version="default"):
    """Install useful executables from UCSC.
    """
    from datetime import date
    version = date.today().strftime('%Y%m%d')
    url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
    pkg_name = 'ucsc_tools'
    tools = ["liftOver", "twoBitToFa", "wigToBigWig"]
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    for tool in tools:
        with cd(install_dir):
            if not exists(tool):
                install_cmd = sudo if env.use_sudo else run
                install_cmd("wget %s%s" % (url, tool))
                install_cmd("chmod 755 %s" % tool)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_bowtie(env, version):
    """Install the bowtie short read aligner."""
    mirror_info = "?use_mirror=cdnetworks-us-2"
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie/%s/bowtie-%s-src.zip" % (version, version)
    pkg_name = 'bowtie'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("unzip %s" % os.path.split(url)[-1])
            with cd("bowtie-%s" % version):
                run("make")
                for fname in run("find -perm -100 -name 'bowtie*'").split("\n"):
                    install_cmd("mv -f %s %s" % (fname.strip(), install_dir))
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_bwa(env, version):
    mirror_info = "?use_mirror=cdnetworks-us-1"
    url = "http://downloads.sourceforge.net/project/bio-bwa/bwa-%s.tar.bz2" % (
            version)
    pkg_name = 'bwa'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("tar -xjvpf %s" % (os.path.split(url)[-1]))
            with cd("bwa-%s" % version):
                run("make")
                install_cmd("mv bwa %s" % install_dir)
                install_cmd("mv solid2fastq.pl %s" % install_dir)
                install_cmd("mv qualfa2fq.pl %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_samtools(env, version):
    vext = ""
    mirror_info = "?use_mirror=cdnetworks-us-1"
    url = "http://downloads.sourceforge.net/project/samtools/samtools/%s/" \
            "samtools-%s%s.tar.bz2" % (version, version, vext)
    pkg_name = 'samtools'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("tar -xjvpf %s" % (os.path.split(url)[-1]))
            with cd("samtools-%s%s" % (version, vext)):
                run("sed -i.bak -r -e 's/-lcurses/-lncurses/g' Makefile")
                run("make")
                for install in ["samtools", "misc/maq2sam-long"]:
                    install_cmd("mv -f %s %s" % (install, install_dir))
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_fastx_toolkit(env, version):
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
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_maq(env, version):
    #version = "0.7.1"
    mirror_info = "?use_mirror=cdnetworks-us-1"
    url = "http://downloads.sourceforge.net/project/maq/maq/%s/maq-%s.tar.bz2" \
            % (version, version)
    pkg_name = 'maq'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("tar -xjvpf %s" % (os.path.split(url)[-1]))
            install_cmd = sudo if env.use_sudo else run
            with cd("maq-%s" % version):
                run("./configure --prefix=%s" % (install_dir))
                run("make")
                install_cmd("make install")
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_bfast(env, version):
    major_version_regex = "\d+\.\d+\.\d+"
    major_version = re.search(major_version_regex, version).group(0)
    url = "http://downloads.sourceforge.net/project/bfast/bfast/%s/bfast-%s%s.tar.gz"\
            % (major_version, version)
    pkg_name = 'bfast'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, "%s%s" % (version))
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % (url))
            run("tar -xzvpf %s" % (os.path.split(url)[-1]))
            install_cmd = sudo if env.use_sudo else run
            with cd("bfast-%s%s" % (version)):
                run("./configure --prefix=%s" % (install_dir))
                run("make")
                install_cmd("make install")
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_abyss(env, version):
    url = "http://www.bcgsc.ca/downloads/abyss/abyss-%s.tar.gz" % version
    pkg_name = 'abyss'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % (os.path.split(url)[-1]))
            install_cmd = sudo if env.use_sudo else run
            with cd("abyss-%s" % version):
                # Get boost first
                run("wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2")
                run("tar jxf boost_1_47_0.tar.bz2")
                run("ln -s boost_1_47_0/boost boost")
                run("rm boost_1_47_0.tar.bz2")
                # Get back to abyss
                run("./configure --prefix=%s --with-mpi=/opt/galaxy/pkg/openmpi" % install_dir)
                run("make")
                install_cmd("make install")
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_velvet(env, version):
    url = "http://www.ebi.ac.uk/~zerbino/velvet/velvet_%s.tgz" % version
    pkg_name = "velvet"
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd("velvet_%s" % version):
                run("make")
                for fname in run("find -perm -100 -name 'velvet*'").split("\n"):
                    with settings(warn_only=True):
                        tmp_cmd = "mv -f %s %s" % (fname, install_dir)
                        install_cmd(tmp_cmd)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_macs(env, version):
    url = "https://github.com/downloads/taoliu/MACS/MACS-%s.tar.gz" % version
    major_version = version.split("-")[0]
    pkg_name = "macs"
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            install_cmd = sudo if env.use_sudo else run
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd("MACS-%s" % major_version):
                install_cmd("python setup.py install --prefix %s" % install_dir)
                # SUDO leaves thing in build as root, which causes _make_tmp_dir
                # to not unwind properly if next rm is not present.
                install_cmd("rm -rf build")
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("echo 'PYTHONPATH=%s/lib/python%s/site-packages:$PYTHONPATH' >> %s/env.sh" % (env.python_version, install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_tophat(env, version):
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
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_cufflinks(env, version):
    version = '1.1.0'
    url = 'http://cufflinks.cbcb.umd.edu/downloads/cufflinks-%s.Linux_x86_64.tar.gz' % version
    pkg_name = "cufflinks"
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
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_megablast(env, version):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/%s/blast-%s-x64-linux.tar.gz' % (version, version)
    pkg_name = 'megablast'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd('blast-%s/bin' % version):
                    install_cmd("mv * %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_blast(env, version):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/%s/ncbi-blast-%s-x64-linux.tar.gz' % (version[:-1], version)
    pkg_name = 'blast'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd('ncbi-blast-%s/bin' % version):
                    install_cmd("mv * %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_sputnik(env, version):
    url = 'http://bitbucket.org/natefoo/sputnik-mononucleotide/downloads/sputnik_%s_linux2.6_x86_64' % version
    pkg_name = 'sputnik'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O sputnik %s" % url)
            install_cmd("mv sputnik %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh %s/sputnik" % (install_dir, install_dir))
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_taxonomy(env, version):
    url = 'http://bitbucket.org/natefoo/taxonomy/downloads/taxonomy_%s_linux2.6_x86_64.tar.gz' % version
    pkg_name = 'taxonomy'
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
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_add_scores(env, version):
    url = 'http://bitbucket.org/natefoo/add_scores/downloads/add_scores_%s_linux2.6_x86_64' % version
    pkg_name = 'add_scores'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O add_scores %s" % url)
            install_cmd("mv add_scores %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh %s/add_scores" % (install_dir, install_dir))
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_emboss(env, version):
    url = 'ftp://emboss.open-bio.org/pub/EMBOSS/old/%s/EMBOSS-%s.tar.gz' % (version, version)
    pkg_name = 'emboss'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
                run("./configure --prefix=%s" % install_dir)
                run("make")
                install_cmd("make install")
    phylip_version = '3.6b'
    url = 'ftp://emboss.open-bio.org/pub/EMBOSS/old/%s/PHYLIP-%s.tar.gz' % (version, phylip_version)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd(os.path.split(url)[-1].split('.tar.gz')[0]):
                run("./configure --prefix=%s" % install_dir)
                run("make")
                install_cmd("make install")
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_hyphy(env, version):
    url = 'http://www.datam0nk3y.org/svn/hyphy'
    pkg_name = 'hyphy'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
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
    sudo("touch %s/env.sh" % install_dir)
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_lastz(env, version):
    url = 'http://www.bx.psu.edu/~rsharris/lastz/older/lastz-%s.tar.gz' % version
    pkg_name = 'lastz'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s" % url)
            run("tar -xvzf %s" % os.path.split(url)[-1])
            with cd('lastz-distrib-%s' % version):
                run("sed -i -e 's/GCC_VERSION == 40302/GCC_VERSION >= 40302/' src/quantum.c")
                run("sed -i -e 's/-Werror//' src/Makefile")
                run("make")
                install_cmd("make LASTZ_INSTALL=%s install" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_perm(env, version):
    url = 'http://perm.googlecode.com/files/PerM_Linux64%28noOpenMp%29.gz'
    pkg_name = 'perm'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget -O PerM.gz %s" % url)
            run("gunzip PerM.gz")
            install_cmd("mv PerM %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh %s/PerM" % (install_dir, install_dir))
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_gatk(env, version):
    url = 'ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/GenomeAnalysisTK-%s.tar.bz2' % version
    pkg_name = 'gatk'
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
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
    sudo("chmod +x %s/env.sh" % install_dir)
    # default link
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))
    # Link jar to Galaxy's jar dir
    jar_dir = os.path.join(env.galaxy_jars_dir, pkg_name)
    if not exists(jar_dir):
        install_cmd("mkdir -p %s" % jar_dir)
    tool_dir = os.path.join(env.galaxy_tools_dir, pkg_name, 'default', 'bin')
    install_cmd('ln --force --symbolic %s/*.jar %s/.' % (tool_dir, jar_dir))
    install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, jar_dir))


@if_tool_not_found()
def install_srma(env, version):
    mirror_info = "?use_mirror=voxel"
    url = 'http://downloads.sourceforge.net/project/srma/srma/%s/srma-%s.jar' \
            % (version[:3], version)
    pkg_name = 'srma'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            install_cmd("mv srma-%s.jar %s" % (version, install_dir))
            install_cmd("ln -f -s srma-%s.jar %s/srma.jar" % (version, install_dir))
    sudo("touch %s/env.sh" % install_dir)
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_beam(env, version):
    url = 'http://www.stat.psu.edu/~yuzhang/software/beam2.tar'
    pkg_name = 'beam'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar xf %s" % (os.path.split(url)[-1]))
            install_cmd("mv BEAM2 %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_pass(env, version):
    url = 'http://www.stat.psu.edu/~yuzhang/software/pass2.tar'
    pkg_name = 'pass'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar xf %s" % (os.path.split(url)[-1]))
            install_cmd("mv pass2 %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_lps_tool(env, version):
    url = 'http://www.bx.psu.edu/miller_lab/dist/lps_tool.%s.tar.gz' % version
    pkg_name = 'lps_tool'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
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
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_plink(env, version):
    url = 'http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-%s-x86_64.zip' % version
    pkg_name = 'plink'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("unzip %s" % (os.path.split(url)[-1]))
            install_cmd("mv plink-%s-x86_64/plink %s" % (version, install_dir))
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_fbat(env, version):
    url = 'http://www.biostat.harvard.edu/~fbat/software/fbat%s_linux64.tar.gz' % version.replace('.', '')
    pkg_name = 'fbat'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar zxf %s" % (os.path.split(url)[-1]))
            install_cmd("mv fbat %s" % install_dir)
    sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_haploview(env, version):
    url = 'http://www.broadinstitute.org/ftp/pub/mpg/haploview/Haploview_beta.jar'
    pkg_name = 'haploview'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            install_cmd("mv %s %s" % (os.path.split(url)[-1], install_dir))
            install_cmd("ln -s %s %s/haploview.jar" % (os.path.split(url)[-1], install_dir))
    sudo("touch %s/env.sh" % install_dir)
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_eigenstrat(env, version):
    url = 'http://www.hsph.harvard.edu/faculty/alkes-price/files/EIG%s.tar.gz' % version
    pkg_name = 'eigenstrat'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            run("tar zxf %s" % (os.path.split(url)[-1]))
            install_cmd("mv bin %s" % install_dir)
    sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_mosaik(env, version):
    url = "http://mosaik-aligner.googlecode.com/files/Mosaik-%s-Linux-x64.tar.bz2" % version
    pkg_name = 'mosaik'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s -O %s" % (url, os.path.split(url)[-1]))
            install_cmd("tar -xjvpf %s -C %s" % (os.path.split(url)[-1], install_dir))
    with cd(install_dir):
        with cd("mosaik-aligner"):
            install_cmd("rm -rf data/ MosaikTools/ src/")
        install_cmd("mv mosaik-aligner/* .")
        install_cmd("rm -rf mosaik-aligner")
    install_cmd("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
    install_cmd("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_freebayes(env, version="default"):
    version = time.strftime("%Y-%m-%d")  # set version to today's date considering it's a repo
    url = "git://github.com/ekg/freebayes.git"
    pkg_name = 'freebayes'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            install_cmd("git clone --recursive %s" % url)
            with cd("freebayes"):
                install_cmd("make")
                install_cmd("mv bin/* %s" % install_dir)
    install_cmd("echo 'PATH=%s:$PATH' > %s/env.sh" % (install_dir, install_dir))
    install_cmd("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))


@if_tool_not_found()
def install_picard(env, version):
    mirror_info = "?use_mirror=voxel"
    url = 'http://downloads.sourceforge.net/project/picard/picard-tools/%s/picard-tools-%s.zip' % (version, version)
    pkg_name = 'picard'
    install_dir = os.path.join(env.galaxy_tools_dir, pkg_name, version)
    install_cmd = sudo if env.use_sudo else run
    if not exists(install_dir):
        install_cmd("mkdir -p %s" % install_dir)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget %s%s -O %s" % (url, mirror_info, os.path.split(url)[-1]))
            run("unzip %s" % (os.path.split(url)[-1]))
            install_cmd("mv picard-tools-%s/*.jar %s" % (version, install_dir))
    sudo("touch %s/env.sh" % install_dir)
    sudo("chmod +x %s/env.sh" % install_dir)
    install_dir_root = os.path.join(env.galaxy_tools_dir, pkg_name)
    if env.galaxy_update_default:
        sudo('ln --symbolic --no-dereference --force %s %s/default' % (install_dir, install_dir_root))
    else:
        sudo('if [ ! -d %s/default ]; then ln -s %s %s/default; fi' % (install_dir_root, install_dir, install_dir_root))
    # set up the jars directory
    jar_dir = os.path.join(env.galaxy_jars_dir, 'picard')
    if not exists(jar_dir):
        install_cmd("mkdir -p %s" % jar_dir)
    tool_dir = os.path.join(env.galaxy_tools_dir, pkg_name, 'default')
    install_cmd('ln --force --symbolic %s/*.jar %s/.' % (tool_dir, jar_dir))
    install_cmd('chown --recursive %s:%s %s' % (env.galaxy_user, env.galaxy_user, jar_dir))


@if_tool_not_found()
def install_fastqc(env, version):
    """ This tool is installed in Galaxy's jars dir """
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
