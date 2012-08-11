"""Install next gen sequencing analysis tools not currently packaged.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_installed, _make_tmp_dir,
                    _get_install, _get_install_local, _make_copy, _configure_make,
                    _java_install,
                    _symlinked_java_version_dir, _fetch_and_unpack, _python_make)

@_if_not_installed("faToTwoBit")
def install_ucsc_tools(env):
    """Useful executables from UCSC.

    todo: install from source to handle 32bit and get more programs
    http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
    """
    tools = ["liftOver", "faToTwoBit", "bedToBigBed",
             "bigBedInfo", "bigBedSummary", "bigBedToBed",
             "bigWigInfo", "bigWigSummary", "bigWigToBedGraph", "bigWigToWig",
             "fetchChromSizes", "wigToBigWig", "faSize", "twoBitInfo",
             "faCount"]
    url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
    install_dir = os.path.join(env.system_install, "bin")
    for tool in tools:
        with cd(install_dir):
            if not exists(tool):
                env.safe_sudo("wget %s%s" % (url, tool))
                env.safe_sudo("chmod a+rwx %s" % tool)

# --- Alignment tools

@_if_not_installed("bowtie")
def install_bowtie(env):
    """The bowtie short read aligner.
    http://bowtie-bio.sourceforge.net/index.shtml
    """
    version = "0.12.7"
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie/%s/" \
          "bowtie-%s-src.zip" % (version, version)
    _get_install(url, env, _make_copy("find -perm -100 -name 'bowtie*'"))

@_if_not_installed("bowtie2")
def install_bowtie2(env):
    """bowtie2 short read aligner, with gap support.
    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    """
    version = "2.0.0-beta6"
    url = "http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/%s/" \
          "bowtie2-%s-source.zip" % (version, version)
    _get_install(url, env, _make_copy("find -perm -100 -name 'bowtie2*'"))

@_if_not_installed("bwa")
def install_bwa(env):
    """BWA:  aligns short nucleotide sequences against a long reference sequence.
    http://bio-bwa.sourceforge.net/
    """
    version = "0.5.9"
    url = "http://downloads.sourceforge.net/project/bio-bwa/bwa-%s.tar.bz2" % (
            version)
    def _fix_makefile():
        arch = run("uname -m")
        # if not 64bit, remove the appropriate flag
        if arch.find("x86_64") == -1:
            run("sed -i.bak -r -e 's/-O2 -m64/-O2/g' Makefile")
    _get_install(url, env, _make_copy("ls -1 bwa solid2fastq.pl qualfa2fq.pl",
                                        _fix_makefile))

@_if_not_installed("bfast")
def install_bfast(env):
    """BFAST: Blat-like Fast Accurate Search Tool.
    http://sourceforge.net/apps/mediawiki/bfast/index.php?title=Main_Page
    """
    version = "0.7.0"
    vext = "a"
    url = "http://downloads.sourceforge.net/project/bfast/bfast/%s/bfast-%s%s.tar.gz"\
            % (version, version, vext)
    _get_install(url, env, _configure_make)

@_if_not_installed("perm")
def install_perm(env):
    """Efficient mapping of short sequences accomplished with periodic full sensitive spaced seeds.
    https://code.google.com/p/perm/
    """
    version = "3.6"
    url = "http://perm.googlecode.com/files/PerM_%s_Source.tar.gz" % version
    def gcc44_makefile_patch():
        gcc_cmd = "g++44"
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            result = run("%s -v" % gcc_cmd)
        print result.return_code
        if result.return_code == 0:
            sed("makefile", "g\+\+", gcc_cmd)
    _get_install(url, env, _make_copy("ls -1 perm", gcc44_makefile_patch))

def install_stampy(env):
    """Stampy: mapping of short reads from illumina sequencing machines onto a reference genome.
    http://www.well.ox.ac.uk/project-stampy
    """
    base_version = "1.0.15"
    revision = "1360"
    version = "{0}r{1}".format(base_version, revision)
    url = "http://www.well.ox.ac.uk/~gerton/software/Stampy/" \
          "stampy-{0}.tgz".format(version)
    _get_install_local(url, env, _make_copy(), dir_name="stampy-{0}".format(base_version))

@_if_not_installed("gmap")
def install_gmap(env):
    """GMAP and GSNAP: A Genomic Mapping and Alignment Program for mRNA EST and short reads.
    http://research-pub.gene.com/gmap/
    """
    version = "2011-11-12"
    url = "http://research-pub.gene.com/gmap/src/gmap-gsnap-%s.tar.gz" % version
    _get_install(url, env, _configure_make)

def _wget_with_cookies(ref_url, dl_url):
    run("wget --cookies=on --keep-session-cookies --save-cookies=cookie.txt %s"
            % (ref_url))
    run("wget --referer=%s --cookies=on --load-cookies=cookie.txt "
        "--keep-session-cookies --save-cookies=cookie.txt %s" %
        (ref_url, dl_url))

@_if_not_installed("novoalign")
def install_novoalign(env):
    """Novoalign short read aligner using Needleman-Wunsch algorithm with affine gap penalties.
    http://www.novocraft.com/main/index.php
    """
    base_version = "V2.08.01"
    cs_version = "V1.02.01"
    _url = "http://www.novocraft.com/downloads/%s/" % base_version
    ref_url = "http://www.novocraft.com/main/downloadpage.php"
    base_url = "%s/novocraft%s.gcc.tar.gz" % (_url, base_version)
    cs_url = "%s/novoalignCS%s.gcc.tar.gz" % (_url, cs_version)
    install_dir = os.path.join(env.system_install, "bin")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _wget_with_cookies(ref_url, base_url)
            run("tar -xzvpf novocraft%s.gcc.tar.gz" % base_version)
            with cd("novocraft"):
                for fname in ["isnovoindex", "novo2maq", "novo2paf",
                        "novo2sam.pl", "novoalign", "novobarcode",
                        "novoindex", "novope2bed.pl", "novorun.pl",
                        "novoutil"]:
                    env.safe_sudo("mv %s %s" % (fname, install_dir))
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            _wget_with_cookies(ref_url, cs_url)
            run("tar -xzvpf novoalignCS%s.gcc.tar.gz" % cs_version)
            with cd("novoalignCS"):
                for fname in ["novoalignCS"]:
                    env.safe_sudo("mv %s %s" % (fname, install_dir))

@_if_not_installed("lastz")
def install_lastz(env):
    """LASTZ sequence alignment program.
    http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html
    """
    version = "1.02.00"
    url = "http://www.bx.psu.edu/miller_lab/dist/" \
          "lastz-%s.tar.gz" % version
    def _remove_werror(env):
        sed("src/Makefile", " -Werror", "")
    _get_install(url, env, _make_copy("find -perm -100 -name 'lastz'"),
                 post_unpack_fn=_remove_werror)

@_if_not_installed("MosaikAligner")
def install_mosaik(env):
    """MOSAIK: reference-guided aligner for next-generation sequencing technologies
    http://code.google.com/p/mosaik-aligner/
    """
    version = "github"
    repository = "git clone git://github.com/wanpinglee/MOSAIK.git"
    def _chdir_src(work_cmd):
        def do_work(env):
            with cd("src"):
                work_cmd(env)
        return do_work
    _get_install(repository, env, _chdir_src(_make_copy("ls -1 ../bin/*")))

# --- Utilities

@_if_not_installed("samtools")
def install_samtools(env):
    """SAM Tools provide various utilities for manipulating alignments in the SAM format.
    http://samtools.sourceforge.net/
    """
    version = "0.1.18"
    url = "http://downloads.sourceforge.net/project/samtools/samtools/" \
          "%s/samtools-%s.tar.bz2" % (version, version)
    _get_install(url, env, _make_copy("find -perm -100 -type f"))

@_if_not_installed("fastq_quality_boxplot_graph.sh")
def install_fastx_toolkit(env):
    """FASTX-Toolkit: collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
    http://hannonlab.cshl.edu/fastx_toolkit/
    """
    version = "0.0.13"
    gtext_version = "0.6"
    url_base = "http://hannonlab.cshl.edu/fastx_toolkit/"
    fastx_url = "%sfastx_toolkit-%s.tar.bz2" % (url_base, version)
    gtext_url = "%slibgtextutils-%s.tar.bz2" % (url_base, gtext_version)
    def _remove_werror(env):
        sed("configure", " -Werror", "")
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
            run("wget %s" % url)
            run("unzip %s" % os.path.basename(url))
            env.safe_sudo("mv SolexaQA.pl %s" % os.path.join(env.system_install, "bin"))

@_if_not_installed("vcftools")
def install_vcftools(env):
    """Work with VCF files, such as those generated by the 1000 Genomes Project.
    http://vcftools.sourceforge.net/
    """
    version = "0.1.7"
    url = "http://downloads.sourceforge.net/project/vcftools/vcftools_{v}.tar.gz".format(
        v=version)
    def _vcf_make(env):
        env.safe_sudo("make install PREFIX={dir}".format(dir=env.system_install))
        for perl_module in ["FaSlice.pm", "Vcf.pm", "VcfStats.pm"]:
            env.safe_sudo("cp perl/%s %s/lib/perl5" % (perl_module, env.system_install))
        env.safe_sudo("make clean")
    _get_install(url, env, _vcf_make)

@_if_not_installed("vtools")
def install_varianttools(env):
    """Annotation, selection, and analysis of variants in the context of next-gen sequencing analysis.
    http://varianttools.sourceforge.net/
    """
    version = "1.0.1"
    version_ext = "a"
    url = "http://downloads.sourceforge.net/project/varianttools/" \
          "{ver}/variant_tools-{ver}{ext}.tar.gz".format(ver=version, ext=version_ext)
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
    version = "0.1.8"
    samtools_version = "0.1.18"
    url = "http://downloads.sourceforge.net/project/dnaa/dwgsim/" \
          "dwgsim-{0}.tar.gz".format(version)
    samtools_url = "http://downloads.sourceforge.net/project/samtools/samtools/" \
                   "{ver}/samtools-{ver}.tar.bz2".format(ver=samtools_version)
    def _get_samtools(env):
        run("wget {0}".format(samtools_url))
        run("tar jxf samtools-{0}.tar.bz2".format(samtools_version))
        run("ln -s samtools-{0} samtools".format(samtools_version))
    _get_install(url, env, _make_copy("ls -1 dwgsim dwgsim_eval scripts/dwgsim_pileup_eval.pl"),
                 post_unpack_fn=_get_samtools)

@_if_not_installed("fastqc")
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
                run("wget %s" % (url))
                run("unzip %s" % os.path.basename(url))
                with cd("FastQC"):
                    env.safe_sudo("chmod a+rwx %s" % executable)
                    env.safe_sudo("mv * %s" % install_dir)
                env.safe_sudo("ln -s %s/%s %s/bin/%s" % (install_dir, executable,
                                                         env.system_install, executable))

@_if_not_installed("bedtools")
def install_bedtools(env):
    """A flexible suite of utilities for comparing genomic features.
    https://code.google.com/p/bedtools/
    """
    version = "github"
    repository = "git clone git://github.com/arq5x/bedtools.git"
    _get_install(repository, env, _make_copy("ls -1 bin/*"))

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
                run("wget %s" % (url))
                run("unzip %s" % os.path.basename(url))
                env.safe_sudo("mv *.class %s" % install_dir)
                for line in _shrec_run.split("\n"):
                    if line.strip():
                        append(shrec_script, line, use_sudo=env.use_sudo)
                env.safe_sudo("chmod a+rwx %s" % shrec_script)
                env.safe_sudo("ln -s %s %s/bin/shrec" % (shrec_script, env.system_install))

def install_echo(env):
    """ECHO: A reference-free short-read error correction algorithm
    http://uc-echo.sourceforge.net/
    """
    version = "1_11"
    url = "http://downloads.sourceforge.net/project/uc-echo/source%20release/" \
          "echo_v{0}.tgz".format(version)
    _get_install_local(url, env, _make_copy())

# -- Analysis

def install_picard(env):
    """Command-line utilities that manipulate BAM files with a Java API.
    http://picard.sourceforge.net/
    """
    version = "1.74"
    url = "http://downloads.sourceforge.net/project/picard/" \
          "picard-tools/%s/picard-tools-%s.zip" % (version, version)
    _java_install("picard", version, url, env)

def install_gatk(env):
    """GATK-lite: library for writing efficient analysis tools using next-generation sequencing data
    http://www.broadinstitute.org/gatk/
    """
    version = "2.0-34-g6d0be9b"
    ext = ".tar.bz2"
    url = "ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/"\
          "GenomeAnalysisTKLite-%s%s" % (version, ext)
    _java_install("gatk", version, url, env)

def install_gatk_queue(env):
    """Command-line scripting framework for defining multi-stage genomic analysis pipelines.
    http://www.broadinstitute.org/gsa/wiki/index.php/GATK-Queue
    """
    version = "2.0-34-g6d0be9b"
    ext = ".tar.bz2"
    url = "ftp://ftp.broadinstitute.org/pub/gsa/Queue/"\
          "QueueLite-%s%s" % (version, ext)
    _java_install("gatk_queue", version, url, env)

def install_snpeff(env):
    """Variant annotation and effect prediction tool.
    http://snpeff.sourceforge.net/
    """
    version = "2_0_5"
    genomes = ["GRCh37.64", "NCBIM37.64", "athalianaTair10"]
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
                    run("sed -i.bak -r -e 's/data_dir = \.\/data\//data_dir = %s\/data/' %s" %
                        (install_dir.replace("/", "\/"), "snpEff.config"))
                    run("chmod a+r *.config")
                    env.safe_sudo("mv *.config %s" % install_dir)
                    data_dir = os.path.join(install_dir, "data")
                    env.safe_sudo("mkdir %s" % data_dir)
                    for org in genomes:
                        if not exists(os.path.join(data_dir, org)):
                            gurl = genome_url_base % (version, version, org)
                            _fetch_and_unpack(gurl, need_dir=False)
                            env.safe_sudo("mv data/%s %s" % (org, data_dir))

@_if_not_installed("freebayes")
def install_freebayes(env):
    """Bayesian haplotype-based polymorphism discovery and genotyping.
    https://github.com/ekg/freebayes
    """
    version = "github"
    repository = "git clone --recursive git://github.com/ekg/freebayes.git"
    def _fix_library_order(env):
        sed("vcflib/tabixpp/Makefile", "-ltabix", "-ltabix -lz")
    _get_install(repository, env, _make_copy("ls -1 bin/*"),
                 post_unpack_fn=_fix_library_order)

def _install_samtools_libs(env):
    repository = "svn co --non-interactive " \
                 "https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
    def _samtools_lib_install(env):
        lib_dir = os.path.join(env.system_install, "lib")
        include_dir = os.path.join(env.system_install, "include", "bam")
        run("make")
        env.safe_sudo("mv -f libbam* %s" % lib_dir)
        env.safe_sudo("mkdir -p %s" % include_dir)
        env.safe_sudo("mv -f *.h %s" % include_dir)
    check_dir = os.path.join(env.system_install, "include", "bam")
    if not exists(check_dir):
        _get_install(repository, env, _samtools_lib_install)

def _install_boost(env):
    version = "1.49.0"
    url = "http://downloads.sourceforge.net/project/boost/boost" \
          "/%s/boost_%s.tar.bz2" % (version, version.replace(".", "_"))
    check_version = "_".join(version.split(".")[:2])
    boost_dir = os.path.join(env.system_install, "boost")
    boost_version_file = os.path.join(boost_dir, "include", "boost", "version.hpp")
    def _boost_build(env):
        run("./bootstrap.sh --prefix=%s --with-libraries=thread" % boost_dir)
        run("./b2")
        env.safe_sudo("./b2 install")
    if not exists(boost_version_file) or not contains(boost_version_file, check_version):
        _get_install(url, env, _boost_build)
        thread_lib = "libboost_thread.so.%s" % version
        final_lib = os.path.join(env.system_install, "lib", thread_lib)
        orig_lib = os.path.join(boost_dir, "lib", thread_lib)
        if not exists(final_lib):
            env.safe_sudo("ln -s %s %s" % (orig_lib, final_lib))

def _cufflinks_configure_make(env):
    orig_eigen = "%s/include/eigen3" % env.system_install
    need_eigen = "%s/include/eigen3/include" % env.system_install
    if not exists(need_eigen):
        env.safe_sudo("ln -s %s %s" % (orig_eigen, need_eigen))
    run("./configure --disable-werror --prefix=%s --with-eigen=%s" \
        " --with-boost=%s/boost" % (env.system_install, orig_eigen, env.system_install))
    run("make")
    env.safe_sudo("make install")

@_if_not_installed("tophat")
def install_tophat(env):
    """TopHat is a fast splice junction mapper for RNA-Seq reads
    http://tophat.cbcb.umd.edu/
    """
    _install_samtools_libs(env)
    _install_boost(env)
    version = "2.0.0"
    url = "http://tophat.cbcb.umd.edu/downloads/tophat-%s.tar.gz" % version
    _get_install(url, env, _cufflinks_configure_make)

@_if_not_installed("cufflinks")
def install_cufflinks(env):
    """Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples.
    http://cufflinks.cbcb.umd.edu/
    """
    _install_samtools_libs(env)
    _install_boost(env)
    version = "2.0.0"
    url = "http://cufflinks.cbcb.umd.edu/downloads/cufflinks-%s.tar.gz" % version
    _get_install(url, env, _cufflinks_configure_make)

# --- Assembly

@_if_not_installed("ABYSS")
def install_abyss(env):
    """Assembly By Short Sequences - a de novo, parallel, paired-end sequence assembler.
    http://www.bcgsc.ca/platform/bioinfo/software/abyss
    """
    # XXX check for no sparehash on non-ubuntu systems
    version = "1.3.3"
    url = "http://www.bcgsc.ca/downloads/abyss/abyss-%s.tar.gz" % version
    def _remove_werror_get_boost(env):
        sed("configure", " -Werror", "")
        # http://osdir.com/ml/abyss-users-science/2011-10/msg00108.html
        run("wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2")
        run("tar jxf boost_1_47_0.tar.bz2")
        run("ln -s boost_1_47_0/boost boost")
    _get_install(url, env, _configure_make, post_unpack_fn=_remove_werror_get_boost)

def install_transabyss(env):
    """Analyze ABySS multi-k-assembled shotgun transcriptome data.
    http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss
    """
    version = "1.3.2"
    ext = "_20120516"
    url = "http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss/" \
          "releases/%s/trans-ABySS-v%s%s.tar.gz" % (version, version, ext)
    _get_install_local(url, env, _make_copy(do_make=False))

@_if_not_installed("velvetg")
def install_velvet(env):
    """Sequence assembler for very short reads.
    http://www.ebi.ac.uk/~zerbino/velvet/
    """
    version = "1.2.05"
    url = "http://www.ebi.ac.uk/~zerbino/velvet/velvet_%s.tgz" % version
    def _fix_library_order(env):
        """Fix library order problem in recent gcc versions
        http://biostar.stackexchange.com/questions/13713/
        error-installing-velvet-assembler-1-1-06-on-ubuntu-server
        """
        sed("Makefile", "Z_LIB_FILES=-lz", "Z_LIB_FILES=-lz -lm")
    _get_install(url, env, _make_copy("find -perm -100 -name 'velvet*'"),
                 post_unpack_fn=_fix_library_order)

def install_trinity(env):
    """Efficient and robust de novo reconstruction of transcriptomes from RNA-seq data.
    http://trinityrnaseq.sourceforge.net/
    """
    version = "r2012-05-18"
    url = "http://downloads.sourceforge.net/project/trinityrnaseq/" \
          "trinityrnaseq_%s.tar.gz" % version
    _get_install_local(url, env, _make_copy())

# --- ChIP-seq

@_if_not_installed("macs14")
def install_macs(env):
    """Model-based Analysis for ChIP-Seq.
    http://liulab.dfci.harvard.edu/MACS/
    """
    version = "1.4.2"
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
        run("make clean")
    _get_install(url, env, _make_copy("ls -1 bin/* scripts/*"),
                 post_unpack_fn=clean_libs)

@_if_not_installed("CRISP.py")
def install_crisp(env):
    """Detect SNPs and short indels from pooled sequencing data.
    https://sites.google.com/site/vibansal/software/crisp/
    """
    version = "5"
    url = "https://sites.google.com/site/vibansal/software/crisp/" \
          "CRISP-linux-v{0}.tar.gz".format(version)
    def _make_executable():
        run("chmod a+x *.py")
    _get_install(url, env, _make_copy("ls -1 CRISP.py crisp_to_vcf.py",
                                      premake_cmd=_make_executable,
                                      do_make=False))

@_if_not_installed("start_tassel.pl")
def install_tassel(env):
    """TASSEL: evaluate traits associations, evolutionary patterns, and linkage disequilibrium.
    http://www.maizegenetics.net/index.php?option=com_content&task=view&id=89&Itemid=119
    """
    version = "3.0"
    url = "http://www.maizegenetics.net/tassel/tassel{0}_standalone.zip".format(version)
    executables = ["start_tassel.pl", "run_pipeline.pl"]
    install_dir = _symlinked_java_version_dir("tassel", version, env)
    if install_dir:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget %s" % (url))
                run("unzip %s" % os.path.basename(url))
                with cd("tassel{0}_standalone".format(version)):
                    for x in executables:
                        sed(x, "^my \$top.*;",
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
    version = "0.998"
    url = "http://creskolab.uoregon.edu/stacks/source/" \
          "stacks-{0}.tar.gz".format(version)
    _get_install(url, env, _configure_make)

@_if_not_installed("sambamba")
def install_sambamba(env):
    """Library for working with SAM/BAM formats written in D programming language
    https://github.com/lomereiter/sambamba/wiki
    """
    version = "0.1.0"
    url = "http://cloud.github.com/downloads/lomereiter/sambamba/"\
          "sambamba-{0}_amd64.deb".format(version)
    if env.distribution in ["ubuntu", "debian"] and env.is_64bit:
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget {0}".format(url))
                env.safe_sudo("sudo dpkg -i {0}".format(
                        os.path.basename(url)))
        
    
