"""Install next gen sequencing analysis tools not currently packaged.
"""
import os
import re

from fabric.api import *
from fabric.contrib.files import *
import yaml

from shared import (_if_not_installed, _make_tmp_dir,
                    _get_install, _get_install_local, _make_copy, _configure_make,
                    _java_install, _python_cmd,
                    _symlinked_java_version_dir, _fetch_and_unpack, _python_make,
                    _get_lib_dir, _get_include_dir, _apply_patch)
from cloudbio.custom import shared, versioncheck

from cloudbio import libraries
from cloudbio.flavor.config import get_config_file


@_if_not_installed(["twoBitToFa", "gtfToGenePred"])
def install_ucsc_tools(env):
    """Useful executables from UCSC.

    todo: install from source to handle 32bit and get more programs
    http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
    """
    tools = ["liftOver", "faToTwoBit", "bedToBigBed",
             "bigBedInfo", "bigBedSummary", "bigBedToBed",
             "bedGraphToBigWig", "bigWigInfo", "bigWigSummary",
             "bigWigToBedGraph", "bigWigToWig",
             "fetchChromSizes", "wigToBigWig", "faSize", "twoBitInfo",
             "twoBitToFa", "faCount", "gtfToGenePred"]
    url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
    _download_executables(env, url, tools)


@_if_not_installed("blat")
def install_kent_tools(env):
    """

    Please note that the Blat source and executables are freely available for
    academic, nonprofit and personal use. Commercial licensing information is
    available on the Kent Informatics website (http://www.kentinformatics.com/).
    """
    tools = ["blat", "gfClient", "gfServer"]
    url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/"
    _download_executables(env, url, tools)


def _download_executables(env, base_url, tools):
    install_dir = shared._get_bin_dir(env)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            for tool in tools:
                final_tool = os.path.join(install_dir, tool)
                if not env.safe_exists(final_tool) and shared._executable_not_on_path(tool):
                    shared._remote_fetch(env, "%s%s" % (base_url, tool))
                    env.safe_sudo("cp -f %s %s" % (tool, install_dir))
                    final_path = os.path.join(install_dir, tool)
                    env.safe_sudo("chmod uga+rx %s" % final_path)

# --- Alignment tools
def install_featurecounts(env):
    """
    featureCounts from the subread package for counting reads mapping to
    genomic features
    """
    default_version = "1.4.4"
    version = env.get("tool_version", default_version)
    if versioncheck.up_to_date(env, "featureCounts", version, stdout_flag="Version"):
        return
    platform = "MacOS" if env.distribution == "macosx" else "Linux"
    url = ("http://downloads.sourceforge.net/project/subread/"
           "subread-%s/subread-%s-%s-x86_64.tar.gz"
           % (version, version, platform))
    _get_install(url, env, _make_copy("find . -type f -perm -100 -name 'featureCounts'",
                                      do_make=False))


@_if_not_installed("bowtie")
def install_bowtie(env):
    """The bowtie short read aligner.
    http://bowtie-bio.sourceforge.net/index.shtml
    """
    default_version = "1.0.0"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie/%s/" \
          "bowtie-%s-src.zip" % (version, version)
    _get_install(url, env, _make_copy("find . -perm -100 -name 'bowtie*'"))

@_if_not_installed("bowtie2")
def install_bowtie2(env):
    """bowtie2 short read aligner, with gap support.
    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    """
    default_version = "2.1.0"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/%s/" \
          "bowtie2-%s-source.zip" % (version, version)
    _get_install(url, env, _make_copy("find . -perm -100 -name 'bowtie2*'"))

@_if_not_installed("bfast")
def install_bfast(env):
    """BFAST: Blat-like Fast Accurate Search Tool.
    http://sourceforge.net/apps/mediawiki/bfast/index.php?title=Main_Page
    """
    default_version = "0.7.0a"
    version = env.get("tool_version", default_version)
    major_version_regex = "\d+\.\d+\.\d+"
    major_version = re.search(major_version_regex, version).group(0)
    url = "http://downloads.sourceforge.net/project/bfast/bfast/%s/bfast-%s.tar.gz"\
            % (major_version, version)
    _get_install(url, env, _configure_make)

@_if_not_installed("perm")
def install_perm(env):
    """Efficient mapping of short sequences accomplished with periodic full sensitive spaced seeds.
    https://code.google.com/p/perm/
    """
    default_version = "4"
    version = env.get("tool_version", default_version)
    url = "http://perm.googlecode.com/files/PerM%sSource.tar.gz" % version
    def gcc44_makefile_patch():
        gcc_cmd = "g++44"
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            result = env.safe_run("%s -v" % gcc_cmd)
        print result.return_code
        if result.return_code == 0:
            env.safe_sed("makefile", "g\+\+", gcc_cmd)
    _get_install(url, env, _make_copy("ls -1 perm", gcc44_makefile_patch))

@_if_not_installed("snap")
def install_snap(env):
    """Scalable Nucleotide Alignment Program
    http://snap.cs.berkeley.edu/
    """
    version = "0.15"
    url = "http://github.com/downloads/amplab/snap/" \
          "snap-%s-linux.tar.gz" % version
    _get_install(url, env, _make_copy("find . -perm -100 -type f", do_make=False))

def install_stampy(env):
    """Stampy: mapping of short reads from illumina sequencing machines onto a reference genome.
    http://www.well.ox.ac.uk/project-stampy
    """
    version = "1.0.21"
    #version = base_version
    #revision = "1654"
    #version = "{0}r{1}".format(base_version, revision)
    #url = "http://www.well.ox.ac.uk/bioinformatics/Software/" \
    #      "stampy-%s.tgz" % (version)
    # Ugh -- Stampy now uses a 'Stampy-latest' download target
    url = "http://www.well.ox.ac.uk/bioinformatics/Software/" \
          "Stampy-latest.tgz"
    def _clean_makefile(env):
        env.safe_sed("makefile", " -Wl", "")
    _get_install_local(url, env, _make_copy(),
                       dir_name="stampy-{0}".format(version),
                       post_unpack_fn=_clean_makefile)

@_if_not_installed("gmap")
def install_gmap(env):
    """GMAP and GSNAP: A Genomic Mapping and Alignment Program for mRNA EST and short reads.
    http://research-pub.gene.com/gmap/
    """
    version = "2012-11-09"
    url = "http://research-pub.gene.com/gmap/src/gmap-gsnap-%s.tar.gz" % version
    _get_install(url, env, _configure_make)

def _wget_with_cookies(ref_url, dl_url):
    env.safe_run("wget --cookies=on --keep-session-cookies --save-cookies=cookie.txt %s"
                 % (ref_url))
    env.safe_run("wget --referer=%s --cookies=on --load-cookies=cookie.txt "
                 "--keep-session-cookies --save-cookies=cookie.txt %s" %
                 (ref_url, dl_url))

@_if_not_installed("novoalign")
def install_novoalign(env):
    """Novoalign short read aligner using Needleman-Wunsch algorithm with affine gap penalties.
    http://www.novocraft.com/main/index.php
    """
    base_version = "V3.00.02"
    cs_version = "V1.03.02"
    _url = "http://www.novocraft.com/downloads/%s/" % base_version
    ref_url = "http://www.novocraft.com/main/downloadpage.php"
    base_url = "%s/novocraft%s.gcc.tar.gz" % (_url, base_version)
    cs_url = "%s/novoalignCS%s.gcc.tar.gz" % (_url, cs_version)
    install_dir = shared._get_bin_dir(env)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _wget_with_cookies(ref_url, base_url)
            env.safe_run("tar -xzvpf novocraft%s.gcc.tar.gz" % base_version)
            with cd("novocraft"):
                for fname in ["isnovoindex", "novo2maq", "novo2paf",
                              "novo2sam.pl", "novoalign", "novobarcode",
                              "novoindex", "novope2bed.pl", "novorun.pl",
                              "novoutil"]:
                    env.safe_sudo("mv %s %s" % (fname, install_dir))
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _wget_with_cookies(ref_url, cs_url)
            env.safe_run("tar -xzvpf novoalignCS%s.gcc.tar.gz" % cs_version)
            with cd("novoalignCS"):
                for fname in ["novoalignCS"]:
                    env.safe_sudo("mv %s %s" % (fname, install_dir))

@_if_not_installed("novosort")
def install_novosort(env):
    """Multithreaded sort and merge for BAM files.
    http://www.novocraft.com/wiki/tiki-index.php?page=Novosort
    """
    base_version = "V3.00.02"
    version = "V1.00.02"
    url = "http://www.novocraft.com/downloads/%s/novosort%s.gcc.tar.gz" % (base_version, version)
    ref_url = "http://www.novocraft.com/main/downloadpage.php"
    install_dir = shared._get_bin_dir(env)
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _wget_with_cookies(ref_url, url)
            env.safe_run("tar -xzvpf novosort%s.gcc.tar.gz" % version)
            with cd("novosort"):
                for fname in ["novosort"]:
                    env.safe_sudo("mv %s %s" % (fname, install_dir))

@_if_not_installed("lastz")
def install_lastz(env):
    """LASTZ sequence alignment program.
    http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html
    """
    default_version = "1.02.00"
    version = env.get("tool_version", default_version)
    url = "http://www.bx.psu.edu/miller_lab/dist/" \
          "lastz-%s.tar.gz" % version
    def _remove_werror(env):
        env.safe_sed("src/Makefile", " -Werror", "")
    _get_install(url, env, _make_copy("find . -perm -100 -name 'lastz'"),
                 post_unpack_fn=_remove_werror)

@_if_not_installed("MosaikAligner")
def install_mosaik(env):
    """MOSAIK: reference-guided aligner for next-generation sequencing technologies
    http://code.google.com/p/mosaik-aligner/
    """
    version = "2.1.73"
    url = "http://mosaik-aligner.googlecode.com/files/" \
          "MOSAIK-%s-binary.tar" % version
    _get_install(url, env, _make_copy("find . -perm -100 -type f", do_make=False))

# --- Utilities

def install_samtools(env):
    """SAM Tools provide various utilities for manipulating alignments in the SAM format.
    http://samtools.sourceforge.net/
    """
    default_version = "0.1.19"
    version = env.get("tool_version", default_version)
    if versioncheck.up_to_date(env, "samtools", version, stdout_flag="Version:"):
        env.logger.info("samtools version {0} is up to date; not installing"
                        .format(version))
        return
    url = "http://downloads.sourceforge.net/project/samtools/samtools/" \
          "%s/samtools-%s.tar.bz2" % (version, version)
    def _safe_ncurses_make(env):
        """Combine samtools, removing ncurses refs if not present on system.
        """
        with settings(warn_only=True):
            result = env.safe_run("make")
        # no ncurses, fix Makefile and rebuild
        if result.failed:
            env.safe_sed("Makefile", "-D_CURSES_LIB=1", "-D_CURSES_LIB=0")
            env.safe_sed("Makefile", "-lcurses", "# -lcurses")
            env.safe_run("make clean")
            env.safe_run("make")
        install_dir = shared._get_bin_dir(env)
        for fname in env.safe_run_output("ls -1 samtools bcftools/bcftools bcftools/vcfutils.pl misc/wgsim").split("\n"):
            env.safe_sudo("cp -f %s %s" % (fname.rstrip("\r"), install_dir))
    _get_install(url, env, _safe_ncurses_make)

def install_gemini(env):
    """A lightweight db framework for disease and population genetics.
    https://github.com/arq5x/gemini
    """
    version = "0.7.0"
    if versioncheck.up_to_date(env, "gemini -v", version, stdout_flag="gemini"):
        return
    elif not shared._executable_not_on_path("gemini -v"):
        env.safe_run("gemini update")
    else:
        iurl = "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py"
        data_dir = os.path.join(env.system_install,
                                "local" if env.system_install.find("/local") == -1 else "",
                                "share", "gemini")
        with _make_tmp_dir(ext="-gemini") as work_dir:
            with cd(work_dir):
                if env.safe_exists(os.path.basename(iurl)):
                    env.safe_run("rm -f %s" % os.path.basename(iurl))
                installer = shared._remote_fetch(env, iurl)
                env.safe_run("%s %s %s %s %s" %
                             (_python_cmd(env), installer, "" if env.use_sudo else "--nosudo",
                              env.system_install, data_dir))
                env.safe_run("rm -f gemini_install.py")

@_if_not_installed("vtools")
def install_varianttools(env):
    """Annotation, selection, and analysis of variants in the context of next-gen sequencing analysis.
    http://varianttools.sourceforge.net/
    """
    version = "1.0.6"
    url = "http://downloads.sourceforge.net/project/varianttools/" \
          "{ver}/variant_tools-{ver}-src.tar.gz".format(ver=version)
    _get_install(url, env, _python_make)

@_if_not_installed("dwgsim")
def install_dwgsim(env):
    """DWGSIM: simulating NGS data and evaluating mappings and variant calling.
    http://sourceforge.net/apps/mediawiki/dnaa/index.php?title=Main_Page
    """
    version = "0.1.10"
    samtools_version = "0.1.18"
    url = "http://downloads.sourceforge.net/project/dnaa/dwgsim/" \
          "dwgsim-{0}.tar.gz".format(version)
    samtools_url = "http://downloads.sourceforge.net/project/samtools/samtools/" \
                   "{ver}/samtools-{ver}.tar.bz2".format(ver=samtools_version)
    def _get_samtools(env):
        shared._remote_fetch(env, samtools_url)
        env.safe_run("tar jxf samtools-{0}.tar.bz2".format(samtools_version))
        env.safe_run("ln -s samtools-{0} samtools".format(samtools_version))
    _get_install(url, env, _make_copy("ls -1 dwgsim dwgsim_eval scripts/dwgsim_pileup_eval.pl"),
                 post_unpack_fn=_get_samtools)

@_if_not_installed("fastq_screen")
def install_fastq_screen(env):
    """A screening application for high througput sequence data.
    http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
    """
    version = "0.4"
    url = "http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/" \
          "fastq_screen_v%s.tar.gz" % version
    install_dir = shared._symlinked_shared_dir("fastqc_screen", version, env)
    executable = "fastq_screen"
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                out_file = shared._remote_fetch(env, url)
                env.safe_run("tar -xzvpf %s" % out_file)
                with cd("fastq_screen_v%s" % version):
                    env.safe_sudo("mv * %s" % install_dir)
                env.safe_sudo("ln -s %s/%s %s/bin/%s" % (install_dir, executable,
                                                         env.system_install, executable))

def install_bedtools(env):
    """A flexible suite of utilities for comparing genomic features.
    https://code.google.com/p/bedtools/
    """
    version = "2.17.0"
    if versioncheck.up_to_date(env, "bedtools --version", version, stdout_flag="bedtools"):
        return
    url = "https://bedtools.googlecode.com/files/" \
          "BEDTools.v%s.tar.gz" % version
    _get_install(url, env, _make_copy("ls -1 bin/*"))

_shrec_run = """
#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($RealBin);
use Getopt::Long;

my @java_args;
my @args;
foreach (@ARGV) {
  if (/^\-X/) {push @java_args,$_;}
  else {push @args,$_;}}
system("java -cp $RealBin @java_args Shrec @args");
"""

@_if_not_installed("shrec")
def install_shrec(env):
    """Shrec is a bioinformatics tool for error correction of HTS read data.
    http://sourceforge.net/projects/shrec-ec/
    """
    version = "2.2"
    url = "http://downloads.sourceforge.net/project/shrec-ec/SHREC%%20%s/bin.zip" % version
    install_dir = _symlinked_java_version_dir("shrec", version, env)
    if install_dir:
        shrec_script = "%s/shrec" % install_dir
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                out_file = shared._remote_fetch(env, url)
                env.safe_run("unzip %s" % out_file)
                env.safe_sudo("mv *.class %s" % install_dir)
                for line in _shrec_run.split("\n"):
                    if line.strip():
                        env.safe_append(shrec_script, line, use_sudo=env.use_sudo)
                env.safe_sudo("chmod a+rwx %s" % shrec_script)
                env.safe_sudo("ln -s %s %s/bin/shrec" % (shrec_script, env.system_install))

def install_echo(env):
    """ECHO: A reference-free short-read error correction algorithm
    http://uc-echo.sourceforge.net/
    """
    version = "1_12"
    url = "http://downloads.sourceforge.net/project/uc-echo/source%20release/" \
          "echo_v{0}.tgz".format(version)
    _get_install_local(url, env, _make_copy())

# -- Analysis

def install_picard(env):
    """Command-line utilities that manipulate BAM files with a Java API.
    http://picard.sourceforge.net/
    """
    version = "1.96"
    url = "http://downloads.sourceforge.net/project/picard/" \
          "picard-tools/%s/picard-tools-%s.zip" % (version, version)
    _java_install("picard", version, url, env)

def install_alientrimmer(env):
    """
    Adapter removal tool
    http://www.ncbi.nlm.nih.gov/pubmed/23912058
    """
    version = "0.3.2"
    url = ("ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/"
           "AlienTrimmer_%s.tar.gz" % version)
    _java_install("AlienTrimmer", version, url, env)

def install_rnaseqc(env):
    """Quality control metrics for RNA-seq data
    https://www.broadinstitute.org/cancer/cga/rna-seqc
    """
    version = "1.1.7"
    url = ("https://github.com/chapmanb/RNA-SeQC/releases/download/"
           "v%s/RNA-SeQC_v%s.jar" % (version, version))
    install_dir = _symlinked_java_version_dir("RNA-SeQC", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                out_file = shared._remote_fetch(env, url)
                env.safe_sudo("mv %s %s" % (out_file, install_dir))

def install_varscan(env):
    """Variant detection in massively parallel sequencing data
    http://varscan.sourceforge.net/
    """
    version = "2.3.7"
    url = "http://downloads.sourceforge.net/project/varscan/VarScan.v%s.jar" % version
    install_dir = _symlinked_java_version_dir("varscan", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                out_file = shared._remote_fetch(env, url)
                env.safe_sudo("mv %s %s" % (out_file, install_dir))

def install_mutect(env):
    version = "1.1.5"
    url = "https://github.com/broadinstitute/mutect/releases/download/" \
          "%s/muTect-%s-bin.zip" % (version, version)
    install_dir = _symlinked_java_version_dir("mutect", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                out_file = shared._remote_fetch(env, url)
                env.safe_run("unzip %s" % out_file)
                env.safe_sudo("mv *.jar version.txt LICENSE* %s" % install_dir)

@_if_not_installed("bam")
def install_bamutil(env):
    """Utilities for working with BAM files, from U of M Center for Statistical Genetics.
    http://genome.sph.umich.edu/wiki/BamUtil
    """
    version = "1.0.7"
    url = "http://genome.sph.umich.edu/w/images/5/5d/BamUtilLibStatGen.%s.tgz" % version
    _get_install(url, env, _make_copy("ls -1 bamUtil/bin/bam"),
                 dir_name="bamUtil_%s" % version)

@_if_not_installed("tabix")
def install_tabix(env):
    """Generic indexer for TAB-delimited genome position files
    http://samtools.sourceforge.net/tabix.shtml
    """
    version = "0.2.6"
    url = "http://downloads.sourceforge.net/project/samtools/tabix/tabix-%s.tar.bz2" % version
    _get_install(url, env, _make_copy("ls -1 tabix bgzip"))

@_if_not_installed("disambiguate.py")
def install_disambiguate(env):
    """a  tool for disambiguating reads aligning to multiple genomes
    https://github.com:mjafin/disambiguate
    """
    repository = "git clone https://github.com/mjafin/disambiguate.git"
    _get_install(repository, env, _python_make)

def install_grabix(env):
    """a wee tool for random access into BGZF files
    https://github.com/arq5x/grabix
    """
    version = "0.1.6"
    revision = "ba792bc872d38d3cb5a69b2de00e39a6ac367d69"
    try:
        uptodate = versioncheck.up_to_date(env, "grabix", version, stdout_flag="version:")
    # Old versions will not have any version information
    except IOError:
        uptodate = False
    if uptodate:
        return
    repository = "git clone https://github.com/arq5x/grabix.git"
    _get_install(repository, env, _make_copy("ls -1 grabix"),
                 revision=revision)

@_if_not_installed("pbgzip")
def install_pbgzip(env):
    """Parallel blocked bgzip -- compatible with bgzip but with thread support.
    https://github.com/nh13/samtools/tree/master/pbgzip
    """
    repository = "git clone https://github.com/chapmanb/samtools.git"
    revision = "2cce3ffa97"
    def _build(env):
        with cd("pbgzip"):
            env.safe_run("make")
            install_dir = shared._get_bin_dir(env)
            env.safe_sudo("cp -f pbgzip %s" % (install_dir))
    _get_install(repository, env, _build, revision=revision)

@_if_not_installed("bamtools")
def install_bamtools(env):
    """command-line toolkit for working with BAM data
    https://github.com/pezmaster31/bamtools
    """
    version = "3fe66b9"
    repository = "git clone --recursive https://github.com/pezmaster31/bamtools.git"
    def _cmake_bamtools(env):
        env.safe_run("mkdir build")
        with cd("build"):
            env.safe_run("cmake ..")
            env.safe_run("make")
        env.safe_sudo("cp bin/* %s" % shared._get_bin_dir(env))
        env.safe_sudo("cp lib/* %s" % shared._get_lib_dir(env))
    _get_install(repository, env, _cmake_bamtools,
                 revision=version)

@_if_not_installed("ogap")
def install_ogap(env):
    """gap opening realigner for BAM data streams
    https://github.com/ekg/ogap
    """
    version = "652c525"
    repository = "git clone --recursive https://github.com/ekg/ogap.git"
    _get_install(repository, env, _make_copy("ls ogap"),
                 revision=version)

def install_tophat(env):
    """TopHat is a fast splice junction mapper for RNA-Seq reads
    http://ccb.jhu.edu/software/tophat/index.shtml
    """
    default_version = "2.0.9"
    version = env.get("tool_version", default_version)
    if versioncheck.is_version(env, "tophat", version, args="--version", stdout_flag="TopHat"):
        env.logger.info("tophat version {0} is up to date; not installing"
            .format(version))
        return
    platform = "OSX" if env.distribution == "macosx" else "Linux"
    url = "http://ccb.jhu.edu/software/tophat/downloads/" \
          "tophat-%s.%s_x86_64.tar.gz" % (version, platform)

    _get_install(url, env,
                 _make_copy("find . -perm -100 -type f", do_make=False))

install_tophat2 = install_tophat

# --- Assembly

@_if_not_installed("ABYSS")
def install_abyss(env):
    """Assembly By Short Sequences - a de novo, parallel, paired-end sequence assembler.
    http://www.bcgsc.ca/platform/bioinfo/software/abyss
    """
    # XXX check for no sparehash on non-ubuntu systems
    default_version = "1.3.4"
    version = env.get("tool_version", default_version)
    url = "http://www.bcgsc.ca/downloads/abyss/abyss-%s.tar.gz" % version
    def _remove_werror_get_boost(env):
        env.safe_sed("configure", " -Werror", "")
        # http://osdir.com/ml/abyss-users-science/2011-10/msg00108.html
        url = "http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2"
        dl_file = shared._remote_fetch(env, url)
        env.safe_run("tar jxf %s" % dl_file)
        env.safe_run("ln -s boost_1_47_0/boost boost")
    _get_install(url, env, _configure_make, post_unpack_fn=_remove_werror_get_boost)

def install_transabyss(env):
    """Analyze ABySS multi-k-assembled shotgun transcriptome data.
    http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss
    """
    version = "1.4.4"
    url = "http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss/" \
          "releases/%s/trans-ABySS-v%s.tar.gz" % (version, version)
    _get_install_local(url, env, _make_copy(do_make=False))

@_if_not_installed("velvetg")
def install_velvet(env):
    """Sequence assembler for very short reads.
    http://www.ebi.ac.uk/~zerbino/velvet/
    """
    default_version = "1.2.08"
    version = env.get("tool_version", default_version)
    url = "http://www.ebi.ac.uk/~zerbino/velvet/velvet_%s.tgz" % version
    def _fix_library_order(env):
        """Fix library order problem in recent gcc versions
        http://biostar.stackexchange.com/questions/13713/
        error-installing-velvet-assembler-1-1-06-on-ubuntu-server
        """
        env.safe_sed("Makefile", "Z_LIB_FILES=-lz", "Z_LIB_FILES=-lz -lm")
    _get_install(url, env, _make_copy("find . -perm -100 -name 'velvet*'"),
                 post_unpack_fn=_fix_library_order)

@_if_not_installed("Ray")
def install_ray(env):
    """Ray -- Parallel genome assemblies for parallel DNA sequencing
    http://denovoassembler.sourceforge.net/
    """
    default_version = "2.2.0"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/denovoassembler/Ray-v%s.tar.bz2" % version
    def _ray_do_nothing(env):
        return
    _get_install(url, env, _make_copy("find . -name Ray"),
                 post_unpack_fn=_ray_do_nothing)

def install_trinity(env):
    """Efficient and robust de novo reconstruction of transcriptomes from RNA-seq data.
    http://trinityrnaseq.github.io/
    """
    version = "2.0.2"
    url = "https://github.com/trinityrnaseq/trinityrnaseq/archive/" \
          "v%s.tar.gz" % version
    dir_name = "trinityrnaseq-%s" % version
    _get_install_local(url, env, _make_copy(),
                       dir_name=dir_name)

def install_cortex_var(env):
    """De novo genome assembly and variation analysis from sequence data.
    http://cortexassembler.sourceforge.net/index_cortex_var.html
    """
    version = "1.0.5.21"
    url = "http://downloads.sourceforge.net/project/cortexassembler/cortex_var/" \
          "latest/CORTEX_release_v{0}.tgz".format(version)
    def _cortex_build(env):
        env.safe_sed("Makefile", "\-L/full/path/\S*",
                     "-L{0}/lib -L/usr/lib -L/usr/local/lib".format(env.system_install))
        env.safe_sed("Makefile", "^IDIR_GSL =.*$",
                     "IDIR_GSL={0}/include -I/usr/include -I/usr/local/include".format(env.system_install))
        env.safe_sed("Makefile", "^IDIR_GSL_ALSO =.*$",
                     "IDIR_GSL_ALSO={0}/include/gsl -I/usr/include/gsl -I/usr/local/include/gsl".format(
                         env.system_install))
        with cd("libs/gsl-1.15"):
            env.safe_run("make clean")
        with cd("libs/htslib"):
            env.safe_run("make clean")
            env.safe_run("make")
        for cols in ["1", "2", "3", "4", "5"]:
            for kmer in ["31", "63", "95"]:
                env.safe_run("make MAXK={0} NUM_COLS={1} cortex_var".format(kmer, cols))
        with cd("scripts/analyse_variants/needleman_wunsch"):
            env.safe_sed("Makefile", "string_buffer.c", "string_buffer.c -lz")
            # Fix incompatibilities with gzfile struct in zlib 1.2.6+
            for fix_gz in ["libs/string_buffer/string_buffer.c", "libs/bioinf/bioinf.c",
                           "libs/string_buffer/string_buffer.h", "libs/bioinf/bioinf.h"]:
                env.safe_sed(fix_gz, "gzFile \*", "gzFile ")
                env.safe_sed(fix_gz, "gzFile\*", "gzFile")
            env.safe_run("make")
    _get_install_local(url, env, _cortex_build)

def install_bcbio_variation(env):
    """Toolkit to analyze genomic variation data with comparison and ensemble approaches.
    https://github.com/chapmanb/bcbio.variation
    """
    version = "0.2.4"
    url = "https://github.com/chapmanb/bcbio.variation/releases/download/" \
          "v%s/bcbio.variation-%s-standalone.jar" % (version, version)
    install_dir = _symlinked_java_version_dir("bcbio_variation", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                jar_file = shared._remote_fetch(env, url)
                env.safe_sudo("mv %s %s" % (jar_file, install_dir))

# --- ChIP-seq

@_if_not_installed("macs14")
def install_macs(env):
    """Model-based Analysis for ChIP-Seq.
    http://liulab.dfci.harvard.edu/MACS/
    """
    default_version = "1.4.2"
    version = env.get("tool_version", default_version)
    url = "https://github.com/downloads/taoliu/MACS/" \
          "MACS-%s.tar.gz" % version
    _get_install(url, env, _python_make)

# --- Structural variation
@_if_not_installed("hydra")
def install_hydra(env):
    """Hydra detects structural variation breakpoints in both unique and duplicated genomic regions.
    https://code.google.com/p/hydra-sv/
    """
    version = "0.5.3"
    url = "http://hydra-sv.googlecode.com/files/Hydra.v{0}.tar.gz".format(version)
    def clean_libs(env):
        env.safe_run("make clean")
    _get_install(url, env, _make_copy("ls -1 bin/* scripts/*"),
                 post_unpack_fn=clean_libs)

def install_freec(env):
    """Control-FREEC: a tool for detection of copy number changes and allelic imbalances.
    http://bioinfo-out.curie.fr/projects/freec/
    """
    version = "6.4"
    if env.distribution in ["ubuntu", "debian"]:
        if env.is_64bit:
            url = "http://bioinfo-out.curie.fr/projects/freec/src/FREEC_Linux64.tar.gz"
        else:
            url = "http://bioinfo-out.curie.fr/projects/freec/src/FREEC_LINUX32.tar.gz"

        if not versioncheck.up_to_date(env, "freec", version, stdout_index=1):
            _get_install(url, env, _make_copy("find . -name 'freec'"), dir_name=".")

@_if_not_installed("CRISP.py")
def install_crisp(env):
    """Detect SNPs and short indels from pooled sequencing data.
    https://sites.google.com/site/vibansal/software/crisp/
    """
    version = "5"
    url = "https://sites.google.com/site/vibansal/software/crisp/" \
          "CRISP-linux-v{0}.tar.gz".format(version)
    def _make_executable():
        env.safe_run("chmod a+x *.py")
    _get_install(url, env, _make_copy("ls -1 CRISP.py crisp_to_vcf.py",
                                      premake_cmd=_make_executable,
                                      do_make=False))

@_if_not_installed("run_pipeline.pl")
def install_tassel(env):
    """TASSEL: evaluate traits associations, evolutionary patterns, and linkage disequilibrium.
    http://www.maizegenetics.net/index.php?option=com_content&task=view&id=89&/Itemid=119
    """
    version = "5"
    build_id = "1140d3fceb75"
    url = "https://bitbucket.org/tasseladmin/tassel-{0}-standalone/get/{1}.zip".format(version, build_id)
    executables = ["start_tassel.pl", "run_pipeline.pl"]
    install_dir = _symlinked_java_version_dir("tassel", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                dl_file = shared._remote_fetch(env, url)
                env.safe_run("unzip %s" % dl_file)
                with cd("tasseladmin-tassel-{0}-standalone-{1}".format(version, build_id)):
                    for x in executables:
                        env.safe_sed(x, "^my \$top.*;",
                                     "use FindBin qw($RealBin); my $top = $RealBin;")
                        env.safe_sudo("chmod a+rwx %s" % x)
                    env.safe_sudo("mv * %s" % install_dir)
                for x in executables:
                    env.safe_sudo("ln -s %s/%s %s/bin/%s" % (install_dir, x,
                                                             env.system_install, x))

@_if_not_installed("ustacks")
def install_stacks(env):
    """Stacks: build loci out of a set of short-read sequenced samples.
    http://creskolab.uoregon.edu/stacks/
    """
    version = "0.9999"
    url = "http://creskolab.uoregon.edu/stacks/source/" \
          "stacks-{0}.tar.gz".format(version)
    _get_install(url, env, _configure_make)

@_if_not_installed("seqlogo")
def install_weblogo(env):
    """Weblogo
    http://weblogo.berkeley.edu/
    """
    version = "2.8.2"
    url = "http://weblogo.berkeley.edu/release/weblogo.%s.tar.gz" % version
    _get_install(url, env, _make_copy("find . -perm -100 -type f", do_make=False))
    def _cp_pm(env):
        for perl_module in ["template.pm", "logo.pm", "template.eps"]:
            env.safe_sudo("cp %s %s/lib/perl5" % (perl_module, env.system_install))
    _get_install(url, env, _cp_pm(env))
