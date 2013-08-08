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
                    _get_lib_dir, _get_include_dir)
from cloudbio.custom import shared, versioncheck

from cloudbio import libraries
from cloudbio.flavor.config import get_config_file


@_if_not_installed("twoBitToFa")
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
             "twoBitToFa", "faCount"]
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
                    env.safe_run("wget %s%s" % (base_url, tool))
                    env.safe_sudo("cp -f %s %s" % (tool, install_dir))

# --- Alignment tools

@_if_not_installed("bowtie")
def install_bowtie(env):
    """The bowtie short read aligner.
    http://bowtie-bio.sourceforge.net/index.shtml
    """
    default_version = "1.0.0"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie/%s/" \
          "bowtie-%s-src.zip" % (version, version)
    _get_install(url, env, _make_copy("find -perm -100 -name 'bowtie*'"))

@_if_not_installed("bowtie2")
def install_bowtie2(env):
    """bowtie2 short read aligner, with gap support.
    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    """
    default_version = "2.1.0"
    version = env.get("tool_version", default_version)
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/%s/" \
          "bowtie2-%s-source.zip" % (version, version)
    _get_install(url, env, _make_copy("find -perm -100 -name 'bowtie2*'"))

def install_bwa(env):
    """BWA:  aligns short nucleotide sequences against a long reference sequence.
    http://bio-bwa.sourceforge.net/
    """
    default_version = "0.7.5a"
    version = env.get("tool_version", default_version)
    if versioncheck.up_to_date(env, "bwa", version, stdout_flag="Version:"):
        return
    url = "http://downloads.sourceforge.net/project/bio-bwa/bwa-%s.tar.bz2" % (
            version)
    def _fix_makefile():
        arch = env.safe_run_output("uname -m")
        # if not 64bit, remove the appropriate flag
        if arch.find("x86_64") == -1:
            env.safe_run("sed -i.bak -r -e 's/-O2 -m64/-O2/g' Makefile")
    _get_install(url, env, _make_copy("ls -1 bwa qualfa2fq.pl",
                                        _fix_makefile))

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
    _get_install(url, env, _make_copy("find -perm -100 -type f", do_make=False))

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
    _get_install(url, env, _make_copy("find -perm -100 -name 'lastz'"),
                 post_unpack_fn=_remove_werror)

@_if_not_installed("MosaikAligner")
def install_mosaik(env):
    """MOSAIK: reference-guided aligner for next-generation sequencing technologies
    http://code.google.com/p/mosaik-aligner/
    """
    version = "2.1.73"
    url = "http://mosaik-aligner.googlecode.com/files/" \
          "MOSAIK-%s-binary.tar" % version
    _get_install(url, env, _make_copy("find -perm -100 -type f", do_make=False))

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

@_if_not_installed("fastq_quality_boxplot_graph.sh")
def install_fastx_toolkit(env):
    """FASTX-Toolkit: collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
    http://hannonlab.cshl.edu/fastx_toolkit/
    """
    default_version = "0.0.13.2"
    version = env.get("tool_version", default_version)
    gtext_version = "0.6"
    url_base = "http://hannonlab.cshl.edu/fastx_toolkit/"
    fastx_url = "%sfastx_toolkit-%s.tar.bz2" % (url_base, version)
    gtext_url = "%slibgtextutils-%s.tar.bz2" % (url_base, gtext_version)
    def _remove_werror(env):
        env.safe_sed("configure", " -Werror", "")
    _get_install(gtext_url, env, _configure_make, post_unpack_fn=_remove_werror)
    _get_install(fastx_url, env, _configure_make, post_unpack_fn=_remove_werror)

@_if_not_installed("SolexaQA.pl")
def install_solexaqa(env):
    """SolexaQA creates visual representations of data quality from FASTQ files.
    http://solexaqa.sourceforge.net/
    """
    version = "1.4"
    url = "http://downloads.sourceforge.net/project/solexaqa/src/" \
            "SolexaQA_v.%s.pl.zip" % version
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            env.safe_run("wget %s" % url)
            env.safe_run("unzip %s" % os.path.basename(url))
            env.safe_sudo("mv SolexaQA.pl %s" % shared._get_bin_dir(env))

@_if_not_installed("gemini -v")
def install_gemini(env):
    """A lightweight db framework for disease and population genetics.
    https://github.com/arq5x/gemini
    """
    version = "github"
    installer = "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py"
    data_dir = os.path.join(env.system_install,
                            "local" if env.system_install.find("/local") == -1 else "",
                            "share", "gemini")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            env.safe_run("wget --no-check-certificate %s" % installer)
            env.safe_run("%s gemini_install.py %s %s %s" %
                         (_python_cmd(env), "" if env.use_sudo else "--nosudo",
                          env.system_install, data_dir))
            env.safe_run("rm -f gemini_install.py")

@_if_not_installed("vcftools")
def install_vcftools(env):
    """Work with VCF files, such as those generated by the 1000 Genomes Project.
    http://vcftools.sourceforge.net/
    """
    version = "0.1.9"
    url = "http://downloads.sourceforge.net/project/vcftools/vcftools_{v}.tar.gz".format(
        v=version)
    def _vcf_make(env):
        env.safe_sudo("make install PREFIX={dir}".format(dir=env.system_install))
        for perl_module in ["FaSlice.pm", "Vcf.pm", "VcfStats.pm"]:
            env.safe_sudo("cp perl/%s %s/lib/perl5" % (perl_module, env.system_install))
        env.safe_sudo("make clean")
    _get_install(url, env, _vcf_make)
    _get_install_local(url, env, _make_copy())

@_if_not_installed("vtools")
def install_varianttools(env):
    """Annotation, selection, and analysis of variants in the context of next-gen sequencing analysis.
    http://varianttools.sourceforge.net/
    """
    version = "1.0.6"
    url = "http://downloads.sourceforge.net/project/varianttools/" \
          "{ver}/variant_tools-{ver}-src.tar.gz".format(ver=version)
    _get_install(url, env, _python_make)

@_if_not_installed("pseq")
def install_plink_seq(env):
    """A toolset for working with human genetic variation data.
    http://atgu.mgh.harvard.edu/plinkseq/
    """
    version = "0.08"
    url = "http://atgu.mgh.harvard.edu/plinkseq/dist/" \
          "version-{v}/plinkseq-{v}-x86_64.tar.gz".format(v=version)
    def _plink_copy(env):
        for x in ["pseq"]:
            env.safe_sudo("cp {0} {1}/bin".format(x, env.system_install))
    _get_install(url, env, _plink_copy)

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
        env.safe_run("wget {0}".format(samtools_url))
        env.safe_run("tar jxf samtools-{0}.tar.bz2".format(samtools_version))
        env.safe_run("ln -s samtools-{0} samtools".format(samtools_version))
    _get_install(url, env, _make_copy("ls -1 dwgsim dwgsim_eval scripts/dwgsim_pileup_eval.pl"),
                 post_unpack_fn=_get_samtools)

@_if_not_installed("fastqc --version")
def install_fastqc(env):
    """A quality control tool for high throughput sequence data.
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """
    version = "0.10.1"
    url = "http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/" \
          "fastqc_v%s.zip" % version
    executable = "fastqc"
    install_dir = _symlinked_java_version_dir("fastqc", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget %s" % (url))
                env.safe_run("unzip %s" % os.path.basename(url))
                with cd("FastQC"):
                    env.safe_sudo("chmod a+rwx %s" % executable)
                    env.safe_sudo("mv * %s" % install_dir)
                env.safe_sudo("ln -s %s/%s %s/bin/%s" % (install_dir, executable,
                                                         env.system_install, executable))

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
                env.safe_run("wget %s" % (url))
                env.safe_run("tar -xzvpf %s" % os.path.basename(url))
                with cd("fastq_screen_v%s" % version):
                    env.safe_sudo("mv * %s" % install_dir)
                env.safe_sudo("ln -s %s/%s %s/bin/%s" % (install_dir, executable,
                                                         env.system_install, executable))

@_if_not_installed("bedtools")
def install_bedtools(env):
    """A flexible suite of utilities for comparing genomic features.
    https://code.google.com/p/bedtools/
    """
    version = "2.17.0"
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
                env.safe_run("wget %s" % (url))
                env.safe_run("unzip %s" % os.path.basename(url))
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
    version = "1.93"
    url = "http://downloads.sourceforge.net/project/picard/" \
          "picard-tools/%s/picard-tools-%s.zip" % (version, version)
    _java_install("picard", version, url, env)

def install_gatk(env):
    """GATK-lite: library for writing efficient analysis tools using next-generation sequencing data
    http://www.broadinstitute.org/gatk/
    """
    # Install main gatk executable
    version = "2.3-9-gdcdccbb"
    ext = ".tar.bz2"
    url = "ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/"\
          "GenomeAnalysisTKLite-%s%s" % (version, ext)
    _java_install("gatk", version, url, env)
    # Install R gsalib for report and pdf generation
    # XXX Currently have issues with gsalib R installation.
    # Need to make this into a proper R package and re-enable
    if False:
        with quiet():
            have_gsalib = env.safe_run("Rscript -e '\"gsalib\" %in% installed.packages()'")
        if have_gsalib and "FALSE" in have_gsalib:
            # install dependencies for gsalib
            rlib_config = get_config_file(env, "r-libs.yaml").base
            with open(rlib_config) as in_handle:
                config = yaml.load(in_handle)
            config["bioc"] = []
            config["update_packages"] = False
            config["cran"] = ["ggplot2", "gplots"]
            libraries.r_library_installer(config)
            # install gsalib
            git_repo = "git clone --depth 1 https://github.com/broadgsa/gatk.git"
            def install_gsalib(env):
                env.safe_sudo("ant gsalib")
            _get_install(git_repo, env, install_gsalib)

def install_varscan(env):
    """Variant detection in massively parallel sequencing data
    http://varscan.sourceforge.net/
    """
    version = "2.3.5"
    url = "http://downloads.sourceforge.net/project/varscan/VarScan.v%s.jar" % version
    install_dir = _symlinked_java_version_dir("varscan", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget --no-check-certificate %s" % url)
                env.safe_sudo("mv *.jar %s" % install_dir)

def install_mutect(env):
    version = "1.1.4"
    url = "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/" \
          "muTect-%s-bin.zip" % version
    install_dir = _symlinked_java_version_dir("mutect", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget --no-check-certificate %s" % url)
                env.safe_run("unzip %s" % os.path.basename(url))
                env.safe_sudo("mv *.jar version.txt LICENSE* %s" % install_dir)

def install_cram(env):
    """Highly efficient and tunable reference-based compression of sequence data.
    http://www.ebi.ac.uk/ena/about/cram_toolkit/
    """
    version = "2.0"
    url = "https://github.com/vadimzalunin/crammer/raw/master/" \
          "cramtools-%s.jar" % version
    install_dir = _symlinked_java_version_dir("cram", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget --no-check-certificate %s" % url)
                env.safe_sudo("mv *.jar %s" % install_dir)

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

@_if_not_installed("grabix")
def install_grabix(env):
    """a wee tool for random access into BGZF files
    https://github.com/arq5x/grabix
    """
    version = "fda4d2609"
    repository = "git clone https://github.com/arq5x/grabix.git"
    _get_install(repository, env, _make_copy("ls -1 grabix"),
                 revision=version)

def install_snpeff(env):
    """Variant annotation and effect prediction tool.
    http://snpeff.sourceforge.net/
    """
    version = "3_3"
    genomes = ["GRCh37.71", "hg19", "GRCm38.71"]
    #genomes_notinstalled = ["NCBIM37.66","athalianaTair10"]
    url = "http://downloads.sourceforge.net/project/snpeff/" \
          "snpEff_v%s_core.zip" % version
    genome_url_base = "http://downloads.sourceforge.net/project/snpeff/"\
                      "databases/v%s/snpEff_v%s_%s.zip"
    install_dir = _symlinked_java_version_dir("snpeff", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                dir_name = _fetch_and_unpack(url)
                with cd(dir_name):
                    env.safe_sudo("mv *.jar %s" % install_dir)
                    env.safe_run("sed -i.bak -r -e 's/^data_dir.*=.*/data_dir = %s\/data/' %s" %
                                 (install_dir.replace("/", "\/"), "snpEff.config"))
                    env.safe_run("chmod a+r *.config")
                    env.safe_sudo("mv *.config %s" % install_dir)
                    data_dir = os.path.join(install_dir, "data")
                    env.safe_sudo("mkdir %s" % data_dir)
                    for org in genomes:
                        if not env.safe_exists(os.path.join(data_dir, org)):
                            gurl = genome_url_base % (version, version, org)
                            _fetch_and_unpack(gurl, need_dir=False)
                            env.safe_sudo("mv data/%s %s" % (org, data_dir))

def install_vep(env):
    """Variant Effects Predictor (VEP) from Ensembl.
    http://ensembl.org/info/docs/variation/vep/index.html
    """
    version = "branch-ensembl-69"
    url = "http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-tools/scripts/" \
          "variant_effect_predictor.tar.gz?view=tar&root=ensembl" \
          "&pathrev={0}".format(version)
    cache_dbs = "24"
    def _vep_install(env):
        env.safe_sed("INSTALL.pl", 'my \$ok = <>', 'my $ok = "y"')
        env.safe_sed("INSTALL.pl", ", <>\)", ', "{0}")'.format(cache_dbs))
        env.safe_run("export FTP_PASSIVE=1 && perl INSTALL.pl")
    _get_install_local(url, env, _vep_install)

@_if_not_installed("freebayes")
def install_freebayes(env):
    """Bayesian haplotype-based polymorphism discovery and genotyping.
    https://github.com/ekg/freebayes
    """
    version = "296a0fa"
    repository = "git clone --recursive https://github.com/ekg/freebayes.git"
    def _fix_tabixpp_library_order(env):
        env.safe_sed("vcflib/tabixpp/Makefile", "-ltabix", "-ltabix -lz")
    _get_install(repository, env, _make_copy("ls -1 bin/*"),
                 post_unpack_fn=_fix_tabixpp_library_order,
                 revision=version)

@_if_not_installed("vcfallelicprimitives -h")
def install_vcflib(env):
    """Utilities for parsing and manipulating VCF files.
    https://github.com/ekg/vcflib
    """
    version = "06e664c"
    repository = "git clone --recursive https://github.com/ekg/vcflib.git"
    def _fix_tabixpp_library_order(env):
        env.safe_sed("tabixpp/Makefile", "-ltabix", "-ltabix -lz")
    _get_install(repository, env,
                 _make_copy("find -perm -100 -type f -name 'vcf*'"
                            " | grep -v '.sh$' | grep -v '.r$'"),
                 post_unpack_fn=_fix_tabixpp_library_order,
                 revision=version)

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

def _install_samtools_libs(env):
    repository = "svn co --non-interactive " \
                 "https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
    def _samtools_lib_install(env):
        lib_dir = _get_lib_dir(env)
        include_dir = os.path.join(env.system_install, "include", "bam")
        env.safe_run("make")
        env.safe_sudo("mv -f libbam* %s" % lib_dir)
        env.safe_sudo("mkdir -p %s" % include_dir)
        env.safe_sudo("mv -f *.h %s" % include_dir)
    check_dir = os.path.join(_get_include_dir(env), "bam")
    if not env.safe_exists(check_dir):
        _get_install(repository, env, _samtools_lib_install)

def _install_boost(env):
    version = "1.49.0"
    url = "http://downloads.sourceforge.net/project/boost/boost" \
          "/%s/boost_%s.tar.bz2" % (version, version.replace(".", "_"))
    check_version = "_".join(version.split(".")[:2])
    boost_dir = os.path.join(env.system_install, "boost")
    boost_version_file = os.path.join(boost_dir, "include", "boost", "version.hpp")
    def _boost_build(env):
        env.safe_run("./bootstrap.sh --prefix=%s --with-libraries=thread" % boost_dir)
        env.safe_run("./b2")
        env.safe_sudo("./b2 install")
    thread_lib = "libboost_thread.so.%s" % version
    final_thread_lib = os.path.join(env.system_install, "lib", thread_lib)
    if (not env.safe_exists(boost_version_file) or not env.safe_contains(boost_version_file, check_version)
          or not env.safe_exists(final_thread_lib)):
        _get_install(url, env, _boost_build)
        orig_lib = os.path.join(boost_dir, "lib", thread_lib)
        if not env.safe_exists(final_thread_lib):
            env.safe_sudo("ln -s %s %s" % (orig_lib, final_thread_lib))

def _cufflinks_configure_make(env):
    orig_eigen = "%s/include/eigen3" % env.system_install
    need_eigen = "%s/include/eigen3/include" % env.system_install
    if not env.safe_exists(need_eigen):
        env.safe_sudo("ln -s %s %s" % (orig_eigen, need_eigen))
    env.safe_run("./configure --disable-werror --prefix=%s --with-eigen=%s"
                 % (env.system_install, orig_eigen))
    #run("./configure --disable-werror --prefix=%s --with-eigen=%s" \
    #    " --with-boost=%s/boost" % (env.system_install, orig_eigen, env.system_install))
    env.safe_run("make")
    env.safe_sudo("make install")

@_if_not_installed("tophat")
def SRC_install_tophat(env):
    """TopHat is a fast splice junction mapper for RNA-Seq reads
    http://tophat.cbcb.umd.edu/
    """
    _install_samtools_libs(env)
    _install_boost(env)
    default_version = "2.0.7"
    version = env.get("tool_version", default_version)
    url = "http://tophat.cbcb.umd.edu/downloads/tophat-%s.tar.gz" % version
    _get_install(url, env, _cufflinks_configure_make)

@_if_not_installed("cufflinks")
def SRC_install_cufflinks(env):
    """Cufflinks assembles transcripts and tests for differential expression and regulation in RNA-Seq samples.
    http://cufflinks.cbcb.umd.edu/
    """
    _install_samtools_libs(env)
    _install_boost(env)
    default_version = "2.0.2"
    version = env.get("tool_version", default_version)
    url = "http://cufflinks.cbcb.umd.edu/downloads/cufflinks-%s.tar.gz" % version
    _get_install(url, env, _cufflinks_configure_make)

@_if_not_installed("tophat")
def install_tophat(env):
    """TopHat is a fast splice junction mapper for RNA-Seq reads
    http://tophat.cbcb.umd.edu/
    """
    default_version = "2.0.8b"
    version = env.get("tool_version", default_version)
    url = "http://tophat.cbcb.umd.edu/downloads/" \
          "tophat-%s.Linux_x86_64.tar.gz" % version
    _get_install(url, env, _make_copy("find -perm -100 -type f",
                                      do_make=False))

install_tophat2 = install_tophat

@_if_not_installed("cufflinks")
def install_cufflinks(env):
    """Cufflinks assembles transcripts and tests for differential expression and regulation in RNA-Seq samples.
    http://cufflinks.cbcb.umd.edu/
    """
    default_version = "2.1.1"
    version = env.get("tool_version", default_version)
    url = "http://cufflinks.cbcb.umd.edu/downloads/" \
          "cufflinks-%s.Linux_x86_64.tar.gz" % version
    _get_install(url, env, _make_copy("find -perm -100 -type f",
                                      do_make=False))

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
        env.safe_run("wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2")
        env.safe_run("tar jxf boost_1_47_0.tar.bz2")
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
    _get_install(url, env, _make_copy("find -perm -100 -name 'velvet*'"),
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
    _get_install(url, env, _make_copy("find -name Ray"),
                 post_unpack_fn=_ray_do_nothing)

def install_trinity(env):
    """Efficient and robust de novo reconstruction of transcriptomes from RNA-seq data.
    http://trinityrnaseq.sourceforge.net/
    """
    version = "r2012-10-05"
    url = "http://downloads.sourceforge.net/project/trinityrnaseq/" \
          "trinityrnaseq_%s.tgz" % version
    def _remove_werror(env):
        env.safe_sed("trinity-plugins/jellyfish/Makefile.in", " -Werror", "")
    _get_install_local(url, env, _make_copy(),
                       post_unpack_fn=_remove_werror)

def install_cortex_var(env):
    """De novo genome assembly and variation analysis from sequence data.
    http://cortexassembler.sourceforge.net/index_cortex_var.html
    """
    version = "1.0.5.20"
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
    version = "0.0.9"
    url = "https://s3.amazonaws.com/bcbio.variation/" \
          "bcbio.variation-%s-standalone.jar" % version
    install_dir = _symlinked_java_version_dir("bcbio_variation", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget %s" % url)
                env.safe_sudo("mv *.jar %s" % install_dir)

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

@_if_not_installed("lumpy")
def install_lumpy(env):
    """a general probabilistic framework for structural variant discovery
    https://github.com/arq5x/lumpy-sv
    """
    version = "fca4706573"
    repository = "git clone https://github.com/arq5x/lumpy-sv.git"
    _get_install(repository, env, _make_copy("ls -1 bin/*"), revision=version)

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
    version = "4.0"
    url = "http://www.maizegenetics.net/tassel/tassel{0}_standalone.zip".format(version)
    executables = ["start_tassel.pl", "run_pipeline.pl"]
    install_dir = _symlinked_java_version_dir("tassel", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget %s" % (url))
                env.safe_run("unzip %s" % os.path.basename(url))
                with cd("tassel{0}_standalone".format(version)):
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

@_if_not_installed("sambamba")
def install_sambamba(env):
    """Library for working with SAM/BAM formats written in D programming language
    https://github.com/lomereiter/sambamba/wiki
    """
    version = "0.2.9"
    url = "https://github.com/downloads/lomereiter/sambamba/" \
          "sambamba-{0}_amd64.deb".format(version)
    if env.distribution in ["ubuntu", "debian"] and env.is_64bit:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget {0}".format(url))
                env.safe_sudo("sudo dpkg -i {0}".format(
                        os.path.basename(url)))
