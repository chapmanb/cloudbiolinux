"""Download and install structured genome data and aligner index files.

Downloads prepared FASTA, indexes for aligners like BWA, Bowtie and novoalign
and other genome data in automated pipelines. Specify the genomes and aligners
to use in an input biodata.yaml configuration file.

The main targets are fabric functions:

  - install_data -- Install biological data from scratch, including indexing genomes.
  - install_data_s3 -- Install biological data, downloading pre-computed indexes from S3.
  - upload_s3 -- Upload created indexes to biodata S3 bucket.

"""
import os
import operator
import socket
import subprocess
from math import log

from fabric.api import *
from fabric.contrib.files import *
from fabric.context_managers import path
try:
    import yaml
except ImportError:
    yaml = None
try:
    import boto
except ImportError:
    boto = None

from cloudbio.biodata import galaxy, ggd
from cloudbio.biodata.dbsnp import download_dbsnp
from cloudbio.biodata.rnaseq import download_transcripts
from cloudbio.custom import shared
from cloudbio.fabutils import quiet
import multiprocessing as mp

# -- Configuration for genomes to download and prepare

class _DownloadHelper:
    def __init__(self):
        self.config = {}

    def ucsc_name(self):
        return None

    def _exists(self, fname, seq_dir):
        """Check if a file exists in either download or final destination.
        """
        return env.safe_exists(fname) or env.safe_exists(os.path.join(seq_dir, fname))

class UCSCGenome(_DownloadHelper):
    def __init__(self, genome_name, dl_name=None):
        _DownloadHelper.__init__(self)
        self.data_source = "UCSC"
        self._name = genome_name
        self.dl_name = dl_name if dl_name is not None else genome_name
        self._url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips" % \
                genome_name

    def ucsc_name(self):
        return self._name

    def _karyotype_sort(self, xs):
        """Sort reads in karyotypic order to work with GATK's defaults.
        """
        def karyotype_keyfn(x):
            base = os.path.splitext(os.path.basename(x))[0]
            if base.startswith("chr"):
                base = base[3:]
            parts = base.split("_")
            try:
                parts[0] = int(parts[0])
            except ValueError:
                pass
            # unplaced at the very end
            if parts[0] == "Un":
                parts.insert(0, "z")
            # mitochondrial special case -- after X/Y
            elif parts[0] in ["M", "MT"]:
                parts.insert(0, "x")
            # sort random and extra chromosomes after M
            elif len(parts) > 1:
                parts.insert(0, "y")
            return parts
        return sorted(xs, key=karyotype_keyfn)

    def _split_multifasta(self, fasta_file):
        chrom = ""
        file_handle = None
        file_names = []
        out_dir = os.path.dirname(fasta_file)
        with open(fasta_file) as in_handle:
            for line in in_handle:
                if line.startswith(">"):
                    chrom = line.split(">")[1].strip()
                    file_handle.close() if file_handle else None
                    file_names.append(chrom + ".fa")
                    file_handle = open(os.path.join(out_dir, chrom + ".fa"), "w")
                    file_handle.write(line)
                else:
                    file_handle.write(line)
        file_handle.close()
        return file_names

    def download(self, seq_dir):
        zipped_file = None
        genome_file = "%s.fa" % self._name
        if not self._exists(genome_file, seq_dir):
            prep_dir = "seq_prep"
            env.safe_run("mkdir -p %s" % prep_dir)
            with cd(prep_dir):
                zipped_file = self._download_zip(seq_dir)
                if zipped_file.endswith(".tar.gz"):
                    env.safe_run("tar -xzpf %s" % zipped_file)
                elif zipped_file.endswith(".zip"):
                    env.safe_run("unzip %s" % zipped_file)
                elif zipped_file.endswith(".gz"):
                    if not env.safe_exists("out.fa"):
                        env.safe_run("gunzip -c %s > out.fa" % zipped_file)
                else:
                    raise ValueError("Do not know how to handle: %s" % zipped_file)
                tmp_file = genome_file.replace(".fa", ".txt")
                result = env.safe_run_output("find `pwd` -name '*.fa'")
                result = [x.strip() for x in result.split("\n")]
                if len(result) == 1:
                    orig_result = result[0]
                    result = self._split_multifasta(result[0])
                    env.safe_run("rm %s" % orig_result)
                result = self._karyotype_sort(result)
                env.safe_run("rm -f inputs.txt")
                for fname in result:
                    with quiet():
                        env.safe_run("echo '%s' >> inputs.txt" % fname)
                env.safe_run("cat `cat inputs.txt` > %s" % (tmp_file))
                for fname in result:
                    with quiet():
                        env.safe_run("rm -f %s" % fname)
                env.safe_run("mv %s %s" % (tmp_file, genome_file))
                zipped_file = os.path.join(prep_dir, zipped_file)
                genome_file = os.path.join(prep_dir, genome_file)
        return genome_file, [zipped_file]

    def _download_zip(self, seq_dir):
        for zipped_file in ["chromFa.tar.gz", "%s.fa.gz" % self._name,
                            "chromFa.zip"]:
            if not self._exists(zipped_file, seq_dir):
                result = shared._remote_fetch(env, "%s/%s" % (self._url, zipped_file), allow_fail=True)
                if result:
                    break
            else:
                break
        return zipped_file

class NCBIRest(_DownloadHelper):
    """Retrieve files using the TogoWS REST server pointed at NCBI.
    """
    def __init__(self, name, refs, dl_name=None):
        _DownloadHelper.__init__(self)
        self.data_source = "NCBI"
        self._name = name
        self._refs = refs
        self.dl_name = dl_name if dl_name is not None else name
        self._base_url = "http://togows.dbcls.jp/entry/ncbi-nucleotide/%s.fasta"

    def download(self, seq_dir):
        genome_file = "%s.fa" % self._name
        if not self._exists(genome_file, seq_dir):
            for ref in self._refs:
                shared._remote_fetch(env, self._base_url % ref)
                env.safe_run("ls -l")
                env.safe_sed('%s.fasta' % ref, '^>.*$', '>%s' % ref, '1')
            tmp_file = genome_file.replace(".fa", ".txt")
            env.safe_run("cat *.fasta > %s" % tmp_file)
            env.safe_run("rm -f *.fasta")
            env.safe_run("rm -f *.bak")
            env.safe_run("mv %s %s" % (tmp_file, genome_file))
        return genome_file, []

class VectorBase(_DownloadHelper):
    """Retrieve genomes from VectorBase) """

    def __init__(self, name, genus, species, strain, release, assembly_types):
        _DownloadHelper.__init__(self)
        self._name = name
        self.data_source = "VectorBase"
        self._base_url = ("http://www.vectorbase.org/sites/default/files/ftp/"
                          "downloads/")
        _base_file = ("{genus}-{species}-{strain}_{assembly}"
                      "_{release}.fa.gz")
        self._to_get = []
        for assembly in assembly_types:
            self._to_get.append(_base_file.format(**locals()))

    def download(self, seq_dir):
        print os.getcwd()
        genome_file = "%s.fa" % self._name
        for fn in self._to_get:
            url = self._base_url + fn
            if not self._exists(fn, seq_dir):
                shared._remote_fetch(env, url)
                env.safe_run("gunzip -c %s >> %s" % (fn, genome_file))
        return genome_file, []

class EnsemblGenome(_DownloadHelper):
    """Retrieve genome FASTA files from Ensembl.

    ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/
    arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.22.dna.toplevel.fa.gz
    ftp://ftp.ensembl.org/pub/release-75/fasta/
    caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.75.dna.toplevel.fa.gz
    ftp://ftp.ensemblgenomes.org/pub/bacteria/release-23/bacteria/fasta/
    bacteria_17_collection/pseudomonas_aeruginosa_ucbpp_pa14/dna/
    Pseudomonas_aeruginosa_ucbpp_pa14.GCA_000014625.1.23.dna.toplevel.fa.gz
    """
    def __init__(self, ensembl_section, release, organism, name, subsection=None):
        _DownloadHelper.__init__(self)
        self.data_source = "Ensembl"
        if ensembl_section == "standard":
            url = "ftp://ftp.ensembl.org/pub/"
        else:
            url = "ftp://ftp.ensemblgenomes.org/pub/%s/" % ensembl_section
        url += "release-%s/fasta/" % release
        if subsection:
            url += "%s/" % subsection
        url += "%s/dna/" % organism.lower()
        self._url = url
        if ensembl_section == "standard":
            self._get_file = "%s.%s.dna.toplevel.fa.gz" % (organism, name)
        else:
            self._get_file = "%s.%s.%s.dna.toplevel.fa.gz" % (organism, name, release)
        self._name = name
        self.dl_name = name

    def download(self, seq_dir):
        genome_file = "%s.fa" % self._name
        if not self._exists(self._get_file, seq_dir):
            shared._remote_fetch(env, "%s%s" % (self._url, self._get_file))
        if not self._exists(genome_file, seq_dir):
            env.safe_run("gunzip -c %s > %s" % (self._get_file, genome_file))
        return genome_file, [self._get_file]

class BroadGenome(_DownloadHelper):
    """Retrieve genomes organized and sorted by Broad for use with GATK.

    Uses the UCSC-name compatible versions of the GATK bundles.
    """
    def __init__(self, name, bundle_version, target_fasta, dl_name=None):
        _DownloadHelper.__init__(self)
        self.data_source = "UCSC"
        self._name = name
        self.dl_name = dl_name if dl_name is not None else name
        self._target = target_fasta
        self._ftp_url = "ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/bundle/" + \
                        "{ver}/{org}/".format(ver=bundle_version, org=self.dl_name)

    def download(self, seq_dir):
        org_file = "%s.fa" % self._name
        if not self._exists(org_file, seq_dir):
            shared._remote_fetch(env, "%s%s.gz" % (self._ftp_url, self._target))
            env.safe_run("gunzip %s.gz" % self._target)
            env.safe_run("mv %s %s" % (self._target, org_file))
        return org_file, []

class GGDGenome:
    """Genome with download specified via a GGD recipe.
    """
    def __init__(self, name):
        self._name = name

BROAD_BUNDLE_VERSION = "2.8"
DBSNP_VERSION = "138"

GENOMES_SUPPORTED = [
           ("phiX174", "phix", NCBIRest("phix", ["NC_001422.1"])),
           ("Scerevisiae", "sacCer3", UCSCGenome("sacCer3")),
           ("Mmusculus", "mm10", UCSCGenome("mm10")),
           ("Mmusculus", "mm9", UCSCGenome("mm9")),
           ("Mmusculus", "mm8", UCSCGenome("mm8")),
           ("Hsapiens", "hg18", BroadGenome("hg18", BROAD_BUNDLE_VERSION,
                                            "Homo_sapiens_assembly18.fasta")),
           ("Hsapiens", "hg19", BroadGenome("hg19", BROAD_BUNDLE_VERSION,
                                            "ucsc.hg19.fasta")),
           ("Hsapiens", "GRCh37", BroadGenome("GRCh37", BROAD_BUNDLE_VERSION,
                                              "human_g1k_v37.fasta", "b37")),
           ("Hsapiens", "hg38", GGDGenome("hg38")),
           ("Hsapiens", "hg38-noalt", GGDGenome("hg38-noalt")),
           ("Rnorvegicus", "rn5", UCSCGenome("rn5")),
           ("Rnorvegicus", "rn4", UCSCGenome("rn4")),
           ("Xtropicalis", "xenTro3", UCSCGenome("xenTro3")),
           ("Athaliana", "TAIR10", EnsemblGenome("plants", "26",
                                                 "Arabidopsis_thaliana", "TAIR10")),
           ("Dmelanogaster", "dm3", UCSCGenome("dm3")),
           ("Celegans", "WBcel235", EnsemblGenome("standard", "80",
                                                  "Caenorhabditis_elegans", "WBcel235")),
           ("Mtuberculosis_H37Rv", "mycoTube_H37RV", NCBIRest("mycoTube_H37RV",
               ["NC_000962"])),
           ("Msmegmatis", "92", NCBIRest("92", ["NC_008596.1"])),
           ("Paeruginosa_UCBPP-PA14", "pseudomonas_aeruginosa_ucbpp_pa14",
            EnsemblGenome("bacteria", "26", "Pseudomonas_aeruginosa_ucbpp_pa14",
                          "GCA_000014625.1", "bacteria_17_collection")),
           ("Ecoli", "eschColi_K12", NCBIRest("eschColi_K12", ["U00096.2"])),
           ("Amellifera_Honeybee", "apiMel3", UCSCGenome("apiMel3")),
           ("Cfamiliaris_Dog", "canFam3", UCSCGenome("canFam3")),
           ("Cfamiliaris_Dog", "canFam2", UCSCGenome("canFam2")),
           ("Drerio_Zebrafish", "Zv9", EnsemblGenome("standard", "80", "Danio_rerio", "Zv9")),
           ("Ecaballus_Horse", "equCab2", UCSCGenome("equCab2")),
           ("Fcatus_Cat", "felCat3", UCSCGenome("felCat3")),
           ("Ggallus_Chicken", "galGal3", UCSCGenome("galGal3")),
           ("Tguttata_Zebra_finch", "taeGut1", UCSCGenome("taeGut1")),
           ("Aalbimanus", "AalbS1", VectorBase("AalbS1", "Anopheles",
                                               "albimanus", "STECLA",
                                               "AalbS1", ["SCAFFOLDS"])),
           ("Agambiae", "AgamP3", VectorBase("AgamP3", "Anopheles",
                                               "gambiae", "PEST",
                                               "AgamP3", ["CHROMOSOMES"])),]


GENOME_INDEXES_SUPPORTED = ["bowtie", "bowtie2", "bwa", "maq", "novoalign", "novoalign-cs",
                            "ucsc", "mosaik", "snap", "star"]
DEFAULT_GENOME_INDEXES = ["seq"]

# -- Fabric instructions

def _check_version():
    version = env.version
    if int(version.split(".")[0]) < 1:
        raise NotImplementedError("Please install fabric version 1 or better")

def install_data(config_source, approaches=None):
    """Main entry point for installing useful biological data.
    """
    PREP_FNS = {"s3": _download_s3_index,
                "ggd": _install_with_ggd,
                "raw": _prep_raw_index}
    if approaches is None: approaches = ["raw"]
    ready_approaches = []
    for approach in approaches:
        ready_approaches.append((approach, PREP_FNS[approach]))
    _check_version()
    # Append a potentially custom system install path to PATH so tools are found
    with path(os.path.join(env.system_install, 'bin')):
        genomes, genome_indexes, config = _get_genomes(config_source)
        genome_indexes = [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes] + genome_indexes
        _make_genome_directories(env, genomes)
        download_transcripts(genomes, env)
        _prep_genomes(env, genomes, genome_indexes, ready_approaches)
        _install_additional_data(genomes, genome_indexes, config)

def install_data_s3(config_source):
    """Install data using pre-existing genomes present on Amazon s3.
    """
    _check_version()
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    _make_genome_directories(env, genomes)
    download_transcripts(genomes, env)
    _download_genomes(genomes, genome_indexes)
    _install_additional_data(genomes, genome_indexes, config)

def install_data_rsync(config_source):
    """Install data using pre-existing genomes from Galaxy rsync servers.
    """
    _check_version()
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    # Galaxy stores FASTAs in ucsc format and generates on the fly
    if "ucsc" not in genome_indexes:
        genome_indexes.append("ucsc")
    genome_dir = _make_genome_dir()
    galaxy.rsync_genomes(genome_dir, genomes, genome_indexes)

def upload_s3(config_source):
    """Upload prepared genome files by identifier to Amazon s3 buckets.
    """
    if boto is None:
        raise ImportError("install boto to upload to Amazon s3")
    if env.host != "localhost" and not env.host.startswith(socket.gethostname()):
        raise ValueError("Need to run S3 upload on a local machine")
    _check_version()
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    _data_ngs_genomes(genomes, genome_indexes)
    _upload_genomes(genomes, genome_indexes)


def _install_additional_data(genomes, genome_indexes, config):
    download_dbsnp(genomes, BROAD_BUNDLE_VERSION, DBSNP_VERSION)
    for custom in (config.get("custom") or []):
        _prep_custom_genome(custom, genomes, genome_indexes, env)
    if config.get("install_liftover", False):
        lift_over_genomes = [g.ucsc_name() for (_, _, g) in genomes if g.ucsc_name()]
        _data_liftover(lift_over_genomes)
    if config.get("install_uniref", False):
        _data_uniref()

def _get_genomes(config_source):
    if isinstance(config_source, dict):
        config = config_source
    else:
        if yaml is None:
            raise ImportError("install yaml to read configuration from %s" % config_source)
        with open(config_source) as in_handle:
            config = yaml.load(in_handle)
    genomes = []
    genomes_config = config["genomes"] or []
    env.logger.info("List of genomes to get (from the config file at '{0}'): {1}"
                    .format(config_source, ', '.join(g.get('name', g["dbkey"]) for g in genomes_config)))
    for g in genomes_config:
        ginfo = None
        for info in GENOMES_SUPPORTED:
            if info[1] == g["dbkey"]:
                ginfo = info
                break
        assert ginfo is not None, "Did not find download info for %s" % g["dbkey"]
        name, gid, manager = ginfo
        manager.config = g
        genomes.append((name, gid, manager))
    indexes = config["genome_indexes"] or []
    if "seq" in indexes:
        indexes.remove("seq")
        indexes.insert(0, "seq")
    return genomes, indexes, config

# ## Decorators and context managers

def _if_installed(pname):
    """Run if the given program name is installed.
    """
    def argcatcher(func):
        def decorator(*args, **kwargs):
            if not shared._executable_not_on_path(pname):
                return func(*args, **kwargs)
        return decorator
    return argcatcher

# ## Generic preparation functions

def _make_genome_dir():
    genome_dir = os.path.join(env.data_files, "genomes")
    if not env.safe_exists(genome_dir):
        with settings(warn_only=True):
            result = env.safe_run_output("mkdir -p %s" % genome_dir)
    else:
        result = None
    if result is not None and result.failed:
        env.safe_sudo("mkdir -p %s" % genome_dir)
        env.safe_sudo("chown -R %s %s" % (env.user, genome_dir))
    return genome_dir


def _make_genome_directories(env, genomes):
    genome_dir = _make_genome_dir()
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not env.safe_exists(org_dir):
            env.safe_run('mkdir -p %s' % org_dir)

def _prep_genomes(env, genomes, genome_indexes, retrieve_fns):
    """Prepare genomes with the given indexes, supporting multiple retrieval methods.
    """
    genome_dir = _make_genome_dir()
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not env.safe_exists(org_dir):
            env.safe_run('mkdir -p %s' % org_dir)
        for idx in genome_indexes + manager.config.get("annotations", []):
            with cd(org_dir):
                if not env.safe_exists(idx):
                    finished = False
                    for method, retrieve_fn in retrieve_fns:
                        try:
                            retrieve_fn(env, manager, gid, idx)
                            finished = True
                            break
                        except KeyboardInterrupt:
                            raise
                        except:
                            # Fail on incorrect GGD recipes
                            if idx in manager.config.get("annotations", []) and method == "ggd":
                                raise
                            else:
                                env.logger.info("Genome preparation method {0} failed, trying next".format(method))
                    if not finished:
                        raise IOError("Could not prepare index {0} for {1} by any method".format(idx, gid))
        ref_file = os.path.join(org_dir, "seq", "%s.fa" % gid)
        if not env.safe_exists(ref_file):
            ref_file = os.path.join(org_dir, "seq", "%s.fa" % manager._name)
        assert env.safe_exists(ref_file), ref_file
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(org_dir, ref_file, gid, cur_indexes, manager.config)

# ## Genomes index for next-gen sequencing tools

def _get_ref_seq(env, manager):
    """Check for or retrieve the reference sequence.
    """
    seq_dir = os.path.join(env.cwd, "seq")
    ref_file = os.path.join(seq_dir, "%s.fa" % manager._name)
    if not env.safe_exists(ref_file):
        ref_file, base_zips = manager.download(seq_dir)
        ref_file = _move_seq_files(ref_file, base_zips, seq_dir)
    return ref_file

def _prep_raw_index(env, manager, gid, idx):
    """Prepare genome from raw downloads and indexes.
    """
    env.logger.info("Preparing genome {0} with index {1}".format(gid, idx))
    ref_file = _get_ref_seq(env, manager)
    get_index_fn(idx)(ref_file)

def _data_ngs_genomes(genomes, genome_indexes):
    """Download and create index files for next generation genomes.
    """
    genome_dir = _make_genome_dir()
    for organism, genome, manager in genomes:
        cur_dir = os.path.join(genome_dir, organism, genome)
        env.logger.info("Processing genome {0} and putting it to {1}"\
            .format(organism, cur_dir))
        if not env.safe_exists(cur_dir):
            env.safe_run('mkdir -p %s' % cur_dir)
        with cd(cur_dir):
            if hasattr(env, "remove_old_genomes") and env.remove_old_genomes:
                _clean_genome_directory()
            seq_dir = 'seq'
            ref_file, base_zips = manager.download(seq_dir)
            ref_file = _move_seq_files(ref_file, base_zips, seq_dir)
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(cur_dir, ref_file, genome, cur_indexes, manager.config)

def _index_to_galaxy(work_dir, ref_file, gid, genome_indexes, config):
    """Index sequence files and update associated Galaxy loc files.
    """
    indexes = {}
    with cd(work_dir):
        for idx in genome_indexes:
            index_file = get_index_fn(idx)(ref_file)
            if index_file:
                indexes[idx] = os.path.join(work_dir, index_file)
    galaxy.prep_locs(gid, indexes, config)

class CustomMaskManager:
    """Create a custom genome based on masking an existing genome.
    """
    def __init__(self, custom, config):
        assert custom.has_key("mask")
        self._custom = custom
        self.config = config

    def download(self, seq_dir):
        base_seq = os.path.join(os.pardir, self._custom["base"],
                                "seq", "{0}.fa".format(self._custom["base"]))
        assert env.safe_exists(base_seq)
        mask_file = os.path.basename(self._custom["mask"])
        ready_mask = apply("{0}-complement{1}".format, os.path.splitext(mask_file))
        out_fasta = "{0}.fa".format(self._custom["dbkey"])
        if not env.safe_exists(os.path.join(seq_dir, out_fasta)):
            if not env.safe_exists(mask_file):
                shared._remote_fetch(env, self._custom["mask"])
            if not env.safe_exists(ready_mask):
                env.safe_run("bedtools complement -i {i} -g {g}.fai > {o}".format(
                    i=mask_file, g=base_seq, o=ready_mask))
            if not env.safe_exists(out_fasta):
                env.safe_run("bedtools maskfasta -fi {fi} -bed {bed} -fo {fo}".format(
                    fi=base_seq, bed=ready_mask, fo=out_fasta))
        return out_fasta, [mask_file, ready_mask]

def _prep_custom_genome(custom, genomes, genome_indexes, env):
    """Prepare a custom genome derived from existing genome.
    Allows creation of masked genomes for specific purposes.
    """
    cur_org = None
    cur_manager = None
    for org, gid, manager in genomes:
        if gid == custom["base"]:
            cur_org = org
            cur_manager = manager
            break
    assert cur_org is not None
    _data_ngs_genomes([[cur_org, custom["dbkey"],
                        CustomMaskManager(custom, cur_manager.config)]],
                      genome_indexes)

def _clean_genome_directory():
    """Remove any existing sequence information in the current directory.
    """
    for dirname in GENOME_INDEXES_SUPPORTED + DEFAULT_GENOME_INDEXES:
        if env.safe_exists(dirname):
            env.safe_run("rm -rf %s" % dirname)

def _move_seq_files(ref_file, base_zips, seq_dir):
    if not env.safe_exists(seq_dir):
        env.safe_run('mkdir %s' % seq_dir)
    for move_file in [ref_file] + base_zips:
        if env.safe_exists(move_file):
            env.safe_run("mv %s %s" % (move_file, seq_dir))
    path, fname = os.path.split(ref_file)
    moved_ref = os.path.join(path, seq_dir, fname)
    assert env.safe_exists(moved_ref), moved_ref
    return moved_ref

# ## Indexing for specific aligners

def _index_w_command(dir_name, command, ref_file, pre=None, post=None, ext=None):
    """Low level function to do the indexing and paths with an index command.
    """
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    if ext is not None: index_name += ext
    full_ref_path = os.path.join(os.pardir, ref_file)
    if not env.safe_exists(dir_name):
        env.safe_run("mkdir %s" % dir_name)
        with cd(dir_name):
            if pre:
                full_ref_path = pre(full_ref_path)
            env.safe_run(command.format(ref_file=full_ref_path, index_name=index_name))
            if post:
                post(full_ref_path)
    return os.path.join(dir_name, index_name)

@_if_installed("faToTwoBit")
def _index_twobit(ref_file):
    """Index reference files using 2bit for random access.
    """
    dir_name = "ucsc"
    cmd = "faToTwoBit {ref_file} {index_name}"
    return _index_w_command(dir_name, cmd, ref_file)

def _index_bowtie(ref_file):
    dir_name = "bowtie"
    cmd = "bowtie-build -f {ref_file} {index_name}"
    return _index_w_command(dir_name, cmd, ref_file)

def _index_bowtie2(ref_file):
    dir_name = "bowtie2"
    cmd = "bowtie2-build {ref_file} {index_name}"
    out_suffix = _index_w_command(dir_name, cmd, ref_file)
    bowtie_link = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.path.pardir,
                                                out_suffix + ".fa"))
    relative_ref_file = os.path.relpath(ref_file, os.path.dirname(bowtie_link))
    if not env.safe_exists(bowtie_link):
        env.safe_run("ln -sf %s %s" % (relative_ref_file, bowtie_link))
    return out_suffix

def _index_bwa(ref_file):
    dir_name = "bwa"
    local_ref = os.path.split(ref_file)[-1]
    if not env.safe_exists(dir_name):
        env.safe_run("mkdir %s" % dir_name)
        with cd(dir_name):
            env.safe_run("ln -sf %s" % os.path.join(os.pardir, ref_file))
            with settings(warn_only=True):
                result = env.safe_run("bwa index -a bwtsw %s" % local_ref)
            # work around a bug in bwa indexing for small files
            if result.failed:
                env.safe_run("bwa index %s" % local_ref)
            env.safe_run("rm -f %s" % local_ref)
    return os.path.join(dir_name, local_ref)

def _index_maq(ref_file):
    dir_name = "maq"
    cmd = "maq fasta2bfa {ref_file} {index_name}"
    def link_local(ref_file):
        local = os.path.basename(ref_file)
        env.safe_run("ln -sf {0} {1}".format(ref_file, local))
        return local
    def rm_local(local_file):
        env.safe_run("rm -f {0}".format(local_file))
    return _index_w_command(dir_name, cmd, ref_file, pre=link_local, post=rm_local)

@_if_installed("novoindex")
def _index_novoalign(ref_file):
    dir_name = "novoalign"
    cmd = "novoindex {index_name} {ref_file}"
    return _index_w_command(dir_name, cmd, ref_file)

@_if_installed("novoalignCS")
def _index_novoalign_cs(ref_file):
    dir_name = "novoalign_cs"
    cmd = "novoindex -c {index_name} {ref_file}"
    return _index_w_command(dir_name, cmd, ref_file)

def _index_sam(ref_file):
    (ref_dir, local_file) = os.path.split(ref_file)
    with cd(ref_dir):
        if not env.safe_exists("%s.fai" % local_file):
            env.safe_run("samtools faidx %s" % local_file)
    galaxy.index_picard(ref_file)
    return ref_file

def _index_star(ref_file):
    (ref_dir, local_file) = os.path.split(ref_file)
    gtf_file = os.path.join(ref_dir, os.pardir, "rnaseq", "ref-transcripts.gtf")
    if not os.path.exists(gtf_file):
        print "%s not found, skipping creating the STAR index." % (gtf_file)
        return None
    GenomeLength = os.path.getsize(ref_file)
    Nbases = int(round(min(14, log(GenomeLength, 2)/2 - 2), 0))
    dir_name = os.path.normpath(os.path.join(ref_dir, os.pardir, "star"))
    cpu = mp.cpu_count()
    cmd = ("STAR --genomeDir %s --genomeFastaFiles {ref_file} "
           "--runThreadN %s "
           "--runMode genomeGenerate --sjdbOverhang 99 --sjdbGTFfile %s --genomeSAindexNbases %s" % (dir_name, str(cpu), gtf_file, Nbases))
    return _index_w_command(dir_name, cmd, ref_file)

def _index_snap(ref_file):
    """Snap indexing is computationally expensive. Ask for all cores and need 64Gb of memory.
    """
    dir_name = "snap"
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    org_arg = "-hg19" if index_name in ["hg19", "GRCh37"] else ""
    cmd = "snap index {ref_file} {dir_name} -bSpace {org_arg}"
    if not env.safe_exists(os.path.join(dir_name, "GenomeIndex")):
        env.safe_run(cmd.format(**locals()))
    return dir_name

@_if_installed("MosaikJump")
def _index_mosaik(ref_file):
    hash_size = 15
    dir_name = "mosaik"
    cmd = "MosaikBuild -fr {ref_file} -oa {index_name}"
    def create_jumpdb(ref_file):
        jmp_base = os.path.splitext(os.path.basename(ref_file))[0]
        dat_file = "{0}.dat".format(jmp_base)
        if not env.safe_exists("{0}_keys.jmp".format(jmp_base)):
            cmd = "export MOSAIK_TMP=`pwd` && MosaikJump -hs {hash_size} -ia {ref_file} -out {index_name}".format(
                hash_size=hash_size, ref_file=dat_file, index_name=jmp_base)
            env.safe_run(cmd)
    return _index_w_command(dir_name, cmd, ref_file,
                            post=create_jumpdb, ext=".dat")

# -- Retrieve using GGD recipes

def _install_with_ggd(env, manager, gid, recipe):
    assert env.hosts == ["localhost"], "GGD recipes only work for local runs"
    recipe_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                               os.pardir, os.pardir, "ggd-recipes"))
    recipe_file = os.path.join(recipe_dir, gid, "%s.yaml" % recipe)
    if os.path.exists(recipe_file):
        ggd.install_recipe(env.cwd, recipe_file)
    else:
        raise NotImplementedError("GGD recipe not available for %s %s" % (gid, recipe))

# -- Genome upload and download to Amazon s3 buckets

def _download_s3_index(env, manager, gid, idx):
    env.logger.info("Downloading genome from s3: {0} {1}".format(gid, idx))
    url = "https://s3.amazonaws.com/biodata/genomes/%s-%s.tar.xz" % (gid, idx)
    out_file = shared._remote_fetch(env, url)
    env.safe_run("xz -dc %s | tar -xvpf -" % out_file)
    env.safe_run("rm -f %s" % out_file)

def _download_genomes(genomes, genome_indexes):
    """Download a group of genomes from Amazon s3 bucket.
    """
    genome_dir = _make_genome_dir()
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not env.safe_exists(org_dir):
            env.safe_run('mkdir -p %s' % org_dir)
        for idx in genome_indexes:
            with cd(org_dir):
                if not env.safe_exists(idx):
                    _download_s3_index(env, manager, gid, idx)
        ref_file = os.path.join(org_dir, "seq", "%s.fa" % gid)
        if not env.safe_exists(ref_file):
            ref_file = os.path.join(org_dir, "seq", "%s.fa" % manager._name)
        assert env.safe_exists(ref_file), ref_file
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(org_dir, ref_file, gid, cur_indexes, manager.config)

def _upload_genomes(genomes, genome_indexes):
    """Upload our configured genomes to Amazon s3 bucket.
    """
    conn = boto.connect_s3()
    bucket = conn.create_bucket("biodata")
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, _) in genomes:
        cur_dir = os.path.join(genome_dir, orgname, gid)
        _clean_directory(cur_dir, gid)
        for idx in genome_indexes:
            idx_dir = os.path.join(cur_dir, idx)
            tarball = _tar_directory(idx_dir, "%s-%s" % (gid, idx))
            _upload_to_s3(tarball, bucket)
    bucket.make_public()

def _upload_to_s3(tarball, bucket):
    """Upload the genome tarball to s3.
    """
    upload_script = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir,
                                 "utils", "s3_multipart_upload.py")
    s3_key_name = os.path.join("genomes", os.path.basename(tarball))
    if not bucket.get_key(s3_key_name):
        gb_size = int(run("du -sm %s" % tarball).split()[0]) / 1000.0
        print "Uploading %s %.1fGb" % (s3_key_name, gb_size)
        cl = ["python", upload_script, tarball, bucket.name, s3_key_name, "--public"]
        subprocess.check_call(cl)

def _tar_directory(dir, tar_name):
    """Create a tarball of the directory.
    """
    base_dir, tar_dir = os.path.split(dir)
    tarball = os.path.join(base_dir, "%s.tar.xz" % tar_name)
    if not env.safe_exists(tarball):
        with cd(base_dir):
            env.safe_run("tar -cvpf - %s | xz -zc - > %s" % (tar_dir,
                                                             os.path.basename(tarball)))
    return tarball

def _clean_directory(dir, gid):
    """Clean duplicate files from directories before tar and upload.
    """
    # get rid of softlinks
    bowtie_ln = os.path.join(dir, "bowtie", "%s.fa" % gid)
    maq_ln = os.path.join(dir, "maq", "%s.fa" % gid)
    for to_remove in [bowtie_ln, maq_ln]:
        if env.safe_exists(to_remove):
            env.safe_run("rm -f %s" % to_remove)
    # remove any downloaded original sequence files
    remove_exts = ["*.gz", "*.zip"]
    with cd(os.path.join(dir, "seq")):
        for rext in remove_exts:
            fnames = env.safe_run("find . -name '%s'" % rext)
            for fname in (f.strip() for f in fnames.split("\n") if f.strip()):
                env.safe_run("rm -f %s" % fname)

# == Liftover files

def _data_liftover(lift_over_genomes):
    """Download chain files for running liftOver.

    Does not install liftOver binaries automatically.
    """
    lo_dir = os.path.join(env.data_files, "liftOver")
    if not env.safe_exists(lo_dir):
        env.safe_run("mkdir %s" % lo_dir)
    lo_base_url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/liftOver/%s"
    lo_base_file = "%sTo%s.over.chain.gz"
    for g1 in lift_over_genomes:
        for g2 in [g for g in lift_over_genomes if g != g1]:
            g2u = g2[0].upper() + g2[1:]
            cur_file = lo_base_file % (g1, g2u)
            non_zip = os.path.splitext(cur_file)[0]
            worked = False
            with cd(lo_dir):
                if not env.safe_exists(non_zip):
                    result = shared._remote_fetch(env, "%s" % (lo_base_url % (g1, cur_file)), allow_fail=True)
                    # Lift over back and forths don't always exist
                    # Only move forward if we found the file
                    if result:
                        worked = True
                        env.safe_run("gunzip %s" % result)
            if worked:
                ref_parts = [g1, g2, os.path.join(lo_dir, non_zip)]
                galaxy.update_loc_file("liftOver.loc", ref_parts)

# == UniRef
def _data_uniref():
    """Retrieve and index UniRef databases for protein searches.

    http://www.ebi.ac.uk/uniref/

    These are currently indexed for FASTA searches. Are other indexes desired?
    Should this be separated out and organized by program like genome data?
    This should also check the release note and automatically download and
    replace older versions.
    """
    site = "ftp://ftp.uniprot.org"
    base_url = site + "/pub/databases/uniprot/" \
               "current_release/uniref/%s/%s"
    for uniref_db in ["uniref50", "uniref90", "uniref100"]:
        work_dir = os.path.join(env.data_files, "uniref", uniref_db)
        if not env.safe_exists(work_dir):
            env.safe_run("mkdir -p %s" % work_dir)
        base_work_url = base_url % (uniref_db, uniref_db)
        fasta_url = base_work_url + ".fasta.gz"
        base_file = os.path.splitext(os.path.basename(fasta_url))[0]
        with cd(work_dir):
            if not env.safe_exists(base_file):
                out_file = shared._remote_fetch(env, fasta_url)
                env.safe_run("gunzip %s" % out_file)
                shared._remote_fetch(env, base_work_url + ".release_note")
        _index_blast_db(work_dir, base_file, "prot")

def _index_blast_db(work_dir, base_file, db_type):
    """Index a database using blast+ for similary searching.
    """
    type_to_ext = dict(prot = ("phr", "pal"), nucl = ("nhr", "nal"))
    db_name = os.path.splitext(base_file)[0]
    with cd(work_dir):
        if not reduce(operator.or_,
            (env.safe_exists("%s.%s" % (db_name, ext)) for ext in type_to_ext[db_type])):
            env.safe_run("makeblastdb -in %s -dbtype %s -out %s" %
                         (base_file, db_type, db_name))


def get_index_fn(index):
    """
    return the index function for an index, if it is missing return a function
    that is a no-op
    """
    return INDEX_FNS.get(index, lambda x: None)

INDEX_FNS = {
    "seq" : _index_sam,
    "bwa" : _index_bwa,
    "bowtie": _index_bowtie,
    "bowtie2": _index_bowtie2,
    "maq": _index_maq,
    "mosaik": _index_mosaik,
    "novoalign": _index_novoalign,
    "novoalign_cs": _index_novoalign_cs,
    "ucsc": _index_twobit,
    "star": _index_star,
    "snap": _index_snap,
    }
