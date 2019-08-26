"""Download and install structured genome data and aligner index files.

Downloads prepared FASTA, indexes for aligners like BWA, Bowtie and novoalign
and other genome data in automated pipelines. Specify the genomes and aligners
to use in an input biodata.yaml configuration file.

The main targets are fabric functions:

  - install_data -- Install biological data from scratch, including indexing genomes.
  - install_data_s3 -- Install biological data, downloading pre-computed indexes from S3.
  - upload_s3 -- Upload created indexes to biodata S3 bucket.
"""
from __future__ import print_function
import collections
import os
import operator
import socket
import subprocess
import sys
import traceback
from math import log

try:
    import yaml
except ImportError:
    yaml = None
try:
    import boto
except ImportError:
    boto = None

from cloudbio.biodata import galaxy, ggd, rnaseq
from cloudbio.custom import shared

# -- Configuration for genomes to download and prepare

class _DownloadHelper:
    def __init__(self):
        self.config = {}

    def ucsc_name(self):
        return None

    def _exists(self, fname, seq_dir):
        """Check if a file exists in either download or final destination.
        """
        return os.path.exists(fname) or os.path.exists(os.path.join(seq_dir, fname))

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
                base = int(parts[0])
            except ValueError:
                base = sys.maxsize
            # unplaced at the very end
            if parts[0] == "Un":
                parts.insert(0, "z")
            # mitochondrial special case -- after X/Y
            elif parts[0] in ["M", "MT"]:
                parts.insert(0, "x")
            # sort random and extra chromosomes after M
            elif len(parts) > 1:
                parts.insert(0, "y")
            # standard integers, sort first
            else:
                parts.insert(0, "a")
            return [base] + parts
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
            subprocess.check_call("mkdir -p %s" % prep_dir, shell=True)
            with shared.chdir(prep_dir):
                zipped_file = self._download_zip(seq_dir)
                if zipped_file.endswith(".tar.gz"):
                    subprocess.check_call("tar -xzpf %s" % zipped_file, shell=True)
                elif zipped_file.endswith(".zip"):
                    subprocess.check_call("unzip %s" % zipped_file, shell=True)
                elif zipped_file.endswith(".gz"):
                    if not os.path.exists("out.fa"):
                        subprocess.check_call("gunzip -c %s > out.fa" % zipped_file, shell=True)
                else:
                    raise ValueError("Do not know how to handle: %s" % zipped_file)
                tmp_file = genome_file.replace(".fa", ".txt")
                result = subprocess.check_output("find `pwd` -name '*.fa'", shell=True).decode()
                result = [x.strip() for x in result.split("\n")]
                if len(result) == 1:
                    orig_result = result[0]
                    result = self._split_multifasta(result[0])
                    subprocess.check_call("rm %s" % orig_result, shell=True)
                result = self._karyotype_sort(result)
                subprocess.check_call("rm -f inputs.txt", shell=True)
                for fname in result:
                    subprocess.check_output("echo '%s' >> inputs.txt" % fname, shell=True).decode()
                subprocess.check_call("cat `cat inputs.txt` > %s" % (tmp_file), shell=True)
                for fname in result:
                    subprocess.check_output("rm -f %s" % fname, shell=True).decode()
                subprocess.check_call("mv %s %s" % (tmp_file, genome_file), shell=True)
                zipped_file = os.path.join(prep_dir, zipped_file)
                genome_file = os.path.join(prep_dir, genome_file)
        return genome_file, [zipped_file]

    def _download_zip(self, seq_dir):
        for zipped_file in ["chromFa.tar.gz", "%s.fa.gz" % self._name,
                            "chromFa.zip"]:
            if not self._exists(zipped_file, seq_dir):
                result = shared._remote_fetch(None, "%s/%s" % (self._url, zipped_file), allow_fail=True)
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
                shared._remote_fetch(None, self._base_url % ref)
                subprocess.check_call("ls -l", shell=True)
                subprocess.check_call(r"sed -i 's/^>.*$/>%s/' %s.fasta" % (ref, ref), shell=True)
            tmp_file = genome_file.replace(".fa", ".txt")
            subprocess.check_call("cat *.fasta > %s" % tmp_file, shell=True)
            subprocess.check_call("rm -f *.fasta", shell=True)
            subprocess.check_call("rm -f *.bak", shell=True)
            subprocess.check_call("mv %s %s" % (tmp_file, genome_file), shell=True)
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
        genome_file = "%s.fa" % self._name
        for fn in self._to_get:
            url = self._base_url + fn
            if not self._exists(fn, seq_dir):
                shared._remote_fetch(None, url)
                subprocess.check_call("gunzip -c %s >> %s" % (fn, genome_file), shell=True)
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
            shared._remote_fetch(None, "%s%s" % (self._url, self._get_file))
        if not self._exists(genome_file, seq_dir):
            subprocess.check_call("gunzip -c %s > %s" % (self._get_file, genome_file), shell=True)
        return genome_file, [self._get_file]

class BroadGenome(_DownloadHelper):
    """Retrieve genomes organized and sorted by Broad for use with GATK.

    Uses the UCSC-name compatible versions of the GATK bundles.
    """
    def __init__(self, name, target_fasta, dl_name=None):
        _DownloadHelper.__init__(self)
        self.data_source = "UCSC"
        self._name = name
        self.dl_name = dl_name if dl_name is not None else name
        self._target = target_fasta
        self._ftp_url = "ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/bundle/" + \
                        "{org}/".format(org=self.dl_name)

    def download(self, seq_dir):
        org_file = "%s.fa" % self._name
        if not self._exists(org_file, seq_dir):
            shared._remote_fetch(None, "%s%s.gz" % (self._ftp_url, self._target))
            subprocess.check_call("gunzip %s.gz" % self._target, shell=True)
            subprocess.check_call("mv %s %s" % (self._target, org_file), shell=True)
        return org_file, []

class GGDGenome:
    """Genome with download specified via a GGD recipe.
    """
    def __init__(self, name):
        self._name = name

GENOMES_SUPPORTED = [
           ("phiX174", "phix", NCBIRest("phix", ["NC_001422.1"])),
           ("Scerevisiae", "sacCer3", UCSCGenome("sacCer3")),
           ("Mmusculus", "mm10", UCSCGenome("mm10")),
           ("Mmusculus", "mm9", UCSCGenome("mm9")),
           ("Mmusculus", "mm8", UCSCGenome("mm8")),
           ("Hsapiens", "hg18", BroadGenome("hg18", "Homo_sapiens_assembly18.fasta")),
           ("Hsapiens", "hg19", BroadGenome("hg19", "ucsc.hg19.fasta")),
           ("Hsapiens", "GRCh37", BroadGenome("GRCh37", "human_g1k_v37.fasta", "b37")),
           ("Hsapiens", "hg38", GGDGenome("hg38")),
           ("Hsapiens", "hg38-noalt", GGDGenome("hg38-noalt")),
           ("Rnorvegicus", "rn6", GGDGenome("rn6")),
           ("Rnorvegicus", "rn5", UCSCGenome("rn5")),
           ("Rnorvegicus", "rn4", UCSCGenome("rn4")),
           ("Xtropicalis", "xenTro3", UCSCGenome("xenTro3")),
           ("Athaliana", "TAIR10", EnsemblGenome("plants", "26",
                                                 "Arabidopsis_thaliana", "TAIR10")),
           ("Dmelanogaster", "dm3", UCSCGenome("dm3")),
           ("Dmelanogaster", "BDGP6", GGDGenome("BDGP6")),
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
           ("Drerio_Zebrafish", "GRCz10", EnsemblGenome("standard", "81", "Danio_rerio", "GRCz10")),
           ("Drerio_Zebrafish", "GRCz11", EnsemblGenome("standard", "92", "Danio_rerio", "GRCz11")),
           ("Sscrofa", "Sscrofa11.1", EnsemblGenome("standard", "92", "Sus_scrofa", "Sscrofa11.1")),
           ("Ecaballus_Horse", "equCab2", UCSCGenome("equCab2")),
           ("Fcatus_Cat", "felCat3", UCSCGenome("felCat3")),
           ("Ggallus_Chicken", "galGal4", UCSCGenome("galGal4")),
           ("Tguttata_Zebra_finch", "taeGut1", UCSCGenome("taeGut1")),
           ("Aalbimanus", "AalbS1", VectorBase("AalbS1", "Anopheles",
                                               "albimanus", "STECLA",
                                               "AalbS1", ["SCAFFOLDS"])),
           ("Agambiae", "AgamP3", VectorBase("AgamP3", "Anopheles",
                                             "gambiae", "PEST",
                                             "AgamP3", ["CHROMOSOMES"]))]
GENOME_INDEXES_SUPPORTED = ["bowtie", "bowtie2", "bwa", "maq", "minimap2", "novoalign",
                            "novoalign-cs", "ucsc", "mosaik", "snap", "star",
                            "rtg", "hisat2", "bbmap", "bismark"]
DEFAULT_GENOME_INDEXES = ["seq"]

# -- Fabric instructions

def _check_version(env):
    version = env.version
    if int(version.split(".")[0]) < 1:
        raise NotImplementedError("Please install fabric version 1 or better")

def install_data(config_source, approaches=None):
    """Main entry point for installing useful biological data, back compatible.
    """
    from fabric.api import env
    _check_version(env)
    install_data_local(config_source, env.system_install, env.data_files,
                       env.galaxy_home, env.tool_data_table_conf_file, env.cores, approaches)

def install_data_local(config_source, system_installdir, data_filedir,
                       galaxy_home=None, tool_data_table_conf_file=None,
                       cores=None, approaches=None):
    """Local installation of biological data, avoiding fabric usage.
    """
    if not cores:
        cores = 1
    PREP_FNS = {"s3": _download_s3_index,
                "ggd": _install_with_ggd,
                "raw": _prep_raw_index}
    if approaches is None: approaches = ["ggd", "s3", "raw"]
    ready_approaches = []
    Env = collections.namedtuple("Env", "system_install, galaxy_home, tool_data_table_conf_file, cores")
    env = Env(system_installdir, galaxy_home, tool_data_table_conf_file, cores)
    for approach in approaches:
        ready_approaches.append((approach, PREP_FNS[approach]))
    # Append a potentially custom system install path to PATH so tools are found
    os.environ["PATH"] = "%s/bin:%s" % (os.path.join(system_installdir), os.environ["PATH"])
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes = [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes] + genome_indexes
    _make_genome_directories(genomes, data_filedir)
    rnaseq.cleanup(genomes, data_filedir)
    _prep_genomes(env, genomes, genome_indexes, ready_approaches, data_filedir)
    rnaseq.finalize(genomes, data_filedir)

def install_data_s3(config_source):
    """Install data using pre-existing genomes present on Amazon s3.
    """
    from fabric.api import env
    _check_version(env)
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    _make_genome_directories(genomes, env.data_files)
    rnaseq.cleanup(genomes, env.data_files)
    _download_genomes(env, genomes, genome_indexes)
    rnaseq.finalize(genomes, env.data_files)
    _install_additional_data(env, genomes, genome_indexes, config)

def install_data_rsync(config_source):
    """Install data using pre-existing genomes from Galaxy rsync servers.
    """
    from fabric.api import env
    _check_version(env)
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    # Galaxy stores FASTAs in ucsc format and generates on the fly
    if "ucsc" not in genome_indexes:
        genome_indexes.append("ucsc")
    genome_dir = _make_genome_dir(env.data_files)
    galaxy.rsync_genomes(genome_dir, genomes, genome_indexes)

def upload_s3(config_source):
    """Upload prepared genome files by identifier to Amazon s3 buckets.
    """
    from fabric.api import env
    if boto is None:
        raise ImportError("install boto to upload to Amazon s3")
    if env.host != "localhost" and not env.host.startswith(socket.gethostname()):
        raise ValueError("Need to run S3 upload on a local machine")
    _check_version(env)
    genomes, genome_indexes, config = _get_genomes(config_source)
    genome_indexes += [x for x in DEFAULT_GENOME_INDEXES if x not in genome_indexes]
    _data_ngs_genomes(env, genomes, genome_indexes)
    _upload_genomes(env, genomes, genome_indexes)


def _install_additional_data(env, genomes, genome_indexes, config):
    for custom in (config.get("custom") or []):
        _prep_custom_genome(custom, genomes, genome_indexes, env)
    if config.get("install_liftover", False):
        lift_over_genomes = [g.ucsc_name() for (_, _, g) in genomes if g.ucsc_name()]
        _data_liftover(env, lift_over_genomes)
    if config.get("install_uniref", False):
        _data_uniref(env)

def _get_genomes(config_source):
    if isinstance(config_source, dict):
        config = config_source
    else:
        if yaml is None:
            raise ImportError("install yaml to read configuration from %s" % config_source)
        with open(config_source) as in_handle:
            config = yaml.safe_load(in_handle)
    genomes = []
    genomes_config = config["genomes"] or []
    print("List of genomes to get (from the config file at '{0}'): {1}"
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
            envs = [x for x in args if hasattr(x, "system_install")]
            env = envs[0] if envs else None
            if shared.which(pname, env):
                return func(*args, **kwargs)
        return decorator
    return argcatcher

# ## Generic preparation functions

def _make_genome_dir(data_filedir):
    genome_dir = os.path.join(data_filedir, "genomes")
    subprocess.check_output("mkdir -p %s" % genome_dir, shell=True).decode()
    return genome_dir

def _make_genome_directories(genomes, data_filedir):
    genome_dir = _make_genome_dir(data_filedir)
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not os.path.exists(org_dir):
            subprocess.check_call('mkdir -p %s' % org_dir, shell=True)

def _prep_genomes(env, genomes, genome_indexes, retrieve_fns, data_filedir):
    """Prepare genomes with the given indexes, supporting multiple retrieval methods.
    """
    genome_dir = _make_genome_dir(data_filedir)
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not os.path.exists(org_dir):
            subprocess.check_call('mkdir -p %s' % org_dir, shell=True)
        ggd_recipes = manager.config.get("annotations", []) + manager.config.get("validation", [])
        ggd_recipes += [x for x in manager.config.get("indexes", []) if x in genome_indexes]
        for idx in genome_indexes + ggd_recipes:
            with shared.chdir(org_dir):
                if idx in ggd_recipes or not os.path.exists(idx):
                    finished = False
                    last_exc = None
                    for method, retrieve_fn in retrieve_fns:
                        try:
                            retrieve_fn(env, manager, gid, idx)
                            finished = True
                            break
                        except KeyboardInterrupt:
                            raise
                        except BaseException as e:
                            # Fail on incorrect GGD recipes
                            if idx in ggd_recipes and method == "ggd":
                                raise
                            else:
                                last_exc = traceback.format_exc()
                                print("Moving on to next genome prep method after trying {0}\n{1}".format(
                                      method, str(e)))
                    if not finished:
                        raise IOError("Could not prepare index {0} for {1} by any method\n{2}"
                                      .format(idx, gid, last_exc))
        ref_file = os.path.join(org_dir, "seq", "%s.fa" % gid)
        if not os.path.exists(ref_file):
            ref_file = os.path.join(org_dir, "seq", "%s.fa" % manager._name)
        assert os.path.exists(ref_file), ref_file
        _index_to_galaxy(env, org_dir, ref_file, gid, genome_indexes, manager.config)

# ## Genomes index for next-gen sequencing tools

def _get_ref_seq(manager):
    """Check for or retrieve the reference sequence.
    """
    seq_dir = os.path.join(os.getcwd(), "seq")
    ref_file = os.path.join(seq_dir, "%s.fa" % manager._name)
    if not os.path.exists(ref_file):
        ref_file, base_zips = manager.download(seq_dir)
        ref_file = _move_seq_files(ref_file, base_zips, seq_dir)
    return ref_file

def _prep_raw_index(env, manager, gid, idx):
    """Prepare genome from raw downloads and indexes.
    """
    print("Preparing genome {0} with index {1}".format(gid, idx))
    ref_file = _get_ref_seq(manager)
    get_index_fn(idx)(env, ref_file)

def _data_ngs_genomes(env, genomes, genome_indexes):
    """Download and create index files for next generation genomes.
    """
    genome_dir = _make_genome_dir(env.data_files)
    for organism, genome, manager in genomes:
        cur_dir = os.path.join(genome_dir, organism, genome)
        print("Processing genome {0} and putting it to {1}".format(organism, cur_dir))
        if not os.path.exists(cur_dir):
            subprocess.check_call('mkdir -p %s' % cur_dir, shell=True)
        with shared.chdir(cur_dir):
            if hasattr(env, "remove_old_genomes") and env.remove_old_genomes:
                _clean_genome_directory()
            seq_dir = 'seq'
            ref_file, base_zips = manager.download(seq_dir)
            ref_file = _move_seq_files(ref_file, base_zips, seq_dir)
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(env, cur_dir, ref_file, genome, cur_indexes, manager.config)

def _index_to_galaxy(env, work_dir, ref_file, gid, genome_indexes, config):
    """Index sequence files and update associated Galaxy loc files.
    """
    indexes = {}
    with shared.chdir(work_dir):
        for idx in genome_indexes:
            index_file = get_index_fn(idx)(env, ref_file)
            if index_file:
                indexes[idx] = os.path.join(work_dir, index_file)
    galaxy.prep_locs(env, gid, indexes, config)

class CustomMaskManager:
    """Create a custom genome based on masking an existing genome.
    """
    def __init__(self, custom, config):
        assert "mask" in custom
        self._custom = custom
        self.config = config

    def download(self, seq_dir):
        base_seq = os.path.join(os.pardir, self._custom["base"],
                                "seq", "{0}.fa".format(self._custom["base"]))
        assert os.path.exists(base_seq)
        mask_file = os.path.basename(self._custom["mask"])
        ready_mask = apply("{0}-complement{1}".format, os.path.splitext(mask_file))
        out_fasta = "{0}.fa".format(self._custom["dbkey"])
        if not os.path.exists(os.path.join(seq_dir, out_fasta)):
            if not os.path.exists(mask_file):
                shared._remote_fetch(None, self._custom["mask"])
            if not os.path.exists(ready_mask):
                subprocess.check_call("bedtools complement -i {i} -g {g}.fai > {o}".format(
                    i=mask_file, g=base_seq, o=ready_mask), shell=True)
            if not os.path.exists(out_fasta):
                subprocess.check_call("bedtools maskfasta -fi {fi} -bed {bed} -fo {fo}".format(
                    fi=base_seq, bed=ready_mask, fo=out_fasta), shell=True)
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
    _data_ngs_genomes(env, [[cur_org, custom["dbkey"],
                             CustomMaskManager(custom, cur_manager.config)]],
                      genome_indexes)

def _clean_genome_directory():
    """Remove any existing sequence information in the current directory.
    """
    for dirname in GENOME_INDEXES_SUPPORTED + DEFAULT_GENOME_INDEXES:
        if os.path.exists(dirname):
            subprocess.check_call("rm -rf %s" % dirname, shell=True)

def _move_seq_files(ref_file, base_zips, seq_dir):
    if not os.path.exists(seq_dir):
        subprocess.check_call('mkdir %s' % seq_dir, shell=True)
    for move_file in [ref_file] + base_zips:
        if os.path.exists(move_file):
            subprocess.check_call("mv %s %s" % (move_file, seq_dir), shell=True)
    path, fname = os.path.split(ref_file)
    moved_ref = os.path.join(path, seq_dir, fname)
    assert os.path.exists(moved_ref), moved_ref
    return moved_ref

# ## Indexing for specific aligners

def _index_w_command(env, dir_name, command, ref_file, pre=None, post=None, ext=None):
    """Low level function to do the indexing and paths with an index command.
    """
    path_export = _get_path_export(env)
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    if ext is not None: index_name += ext
    full_ref_path = os.path.join(os.pardir, ref_file)
    if not os.path.exists(dir_name):
        subprocess.check_call("mkdir %s" % dir_name, shell=True)
        with shared.chdir(dir_name):
            if pre:
                full_ref_path = pre(full_ref_path)
            subprocess.check_call(path_export + command.format(ref_file=full_ref_path, index_name=index_name),
                                  shell=True)
            if post:
                post(full_ref_path)
    return os.path.join(dir_name, index_name)

@_if_installed("faToTwoBit")
def _index_twobit(env, ref_file):
    """Index reference files using 2bit for random access.
    """
    dir_name = "ucsc"
    cmd = "faToTwoBit {ref_file} {index_name}"
    return _index_w_command(env, dir_name, cmd, ref_file)

def _index_bowtie(env, ref_file):
    dir_name = "bowtie"
    cmd = "bowtie-build -f {ref_file} {index_name}"
    return _index_w_command(env, dir_name, cmd, ref_file)

def _index_bowtie2(env, ref_file):
    dir_name = "bowtie2"
    cmd = "bowtie2-build {ref_file} {index_name}"
    out_suffix = _index_w_command(env, dir_name, cmd, ref_file)
    bowtie_link = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.path.pardir,
                                                out_suffix + ".fa"))
    relative_ref_file = os.path.relpath(ref_file, os.path.dirname(bowtie_link))
    if not os.path.exists(bowtie_link):
        subprocess.check_call("ln -sf %s %s" % (relative_ref_file, bowtie_link), shell=True)
    return out_suffix

def _index_bwa(env, ref_file):
    dir_name = "bwa"
    local_ref = os.path.split(ref_file)[-1]
    if not os.path.exists(os.path.join(dir_name, "%s.bwt" % local_ref)):
        subprocess.check_call("mkdir -p %s" % dir_name, shell=True)
        with shared.chdir(dir_name):
            subprocess.check_call("ln -sf %s" % os.path.join(os.pardir, ref_file), shell=True)
            try:
                subprocess.check_call("bwa index -a bwtsw %s" % local_ref, shell=True)
            except subprocess.CalledProcessError:
                # work around a bug in bwa indexing for small files
                subprocess.check_call("bwa index %s" % local_ref, shell=True)
            subprocess.check_call("rm -f %s" % local_ref, shell=True)
    return os.path.join(dir_name, local_ref)

def _index_bbmap(env, ref_file):
    dir_name = "bbmap"
    try:
        cores = env.cores
    except:
        cores = 1
    if not os.path.exists(os.path.join(dir_name, "ref", "genome", "1", "summary.txt")):
        subprocess.check_call("mkdir -p %s" % dir_name, shell=True)
        subprocess.check_call("bbmap.sh -Xms%sg -Xmx%sg path=%s ref=%s" %
                              (cores, 3 * int(cores), dir_name, ref_file), shell=True)
    return dir_name

def _index_bismark(env, ref_file):
    """Bismark indexing happens in GGD recipes. This links to the proper output directory.
    """
    dir_name = "bismark/Bisulfite_Genome"
    return dir_name

def _index_maq(env, ref_file):
    dir_name = "maq"
    cmd = "maq fasta2bfa {ref_file} {index_name}"
    def link_local(ref_file):
        local = os.path.basename(ref_file)
        subprocess.check_call("ln -sf {0} {1}".format(ref_file, local), shell=True)
        return local
    def rm_local(local_file):
        subprocess.check_call("rm -f {0}".format(local_file), shell=True)
    return _index_w_command(env, dir_name, cmd, ref_file, pre=link_local, post=rm_local)

def _index_minimap2(env, ref_file):
    dir_name = "minimap2"
    indexes = []
    for preset in ["sr"]:
        index_name = "%s-%s.mmi" % (os.path.splitext(os.path.basename(ref_file))[0], preset)
        cmd = "minimap2 -x %s -d %s {ref_file}" % (preset, index_name)
        out_basename = _index_w_command(env, dir_name, cmd, ref_file)
        indexes.append(os.path.join(os.path.dirname(out_basename), index_name))
    return indexes[0]

@_if_installed("novoindex")
def _index_novoalign(env, ref_file):
    dir_name = "novoalign"
    cmd = "novoindex {index_name} {ref_file}"
    return _index_w_command(env, dir_name, cmd, ref_file)

@_if_installed("novoalignCS")
def _index_novoalign_cs(env, ref_file):
    dir_name = "novoalign_cs"
    cmd = "novoindex -c {index_name} {ref_file}"
    return _index_w_command(env, dir_name, cmd, ref_file)

def _index_sam(env, ref_file):
    (ref_dir, local_file) = os.path.split(ref_file)
    with shared.chdir(ref_dir):
        if not os.path.exists("%s.fai" % local_file):
            subprocess.check_call("samtools faidx %s" % local_file, shell=True)
    galaxy.index_picard(ref_file)
    return ref_file

@_if_installed("STAR")
def _index_star(env, ref_file):
    (ref_dir, local_file) = os.path.split(ref_file)
    GenomeLength = os.path.getsize(ref_file)
    Nbases = int(round(min(14, log(GenomeLength, 2) / 2 - 2), 0))
    dir_name = os.path.normpath(os.path.join(ref_dir, os.pardir, "star"))
    # if there is a large number of contigs, scale nbits down
    # https://github.com/alexdobin/STAR/issues/103#issuecomment-173009628
    # if there is a small genome, scale nbits down
    # https://groups.google.com/forum/#!topic/rna-star/9g8Uoe1Igho
    cmd = 'grep ">" {ref_file} | wc -l'.format(ref_file=ref_file)
    nrefs = float(subprocess.check_output(cmd, shell=True).decode())
    nbits = int(round(min(14, log(GenomeLength / nrefs, 2), log(GenomeLength, 2) / 2 - 1)))
    # first we estimate the number of bits we need to hold the genome and allocate
    # double that plus some padding to build the index
    mem = ((GenomeLength + 1) / nbits + 1) * nbits
    mem = (mem + 10000) * 2
    mem = mem + mem / 3
    mem = max(mem, 30000000000)
    try:
        cpu = env.cores
    except:
        cpu = 1
    cmd = ("STAR --genomeDir %s --genomeFastaFiles {ref_file} "
           "--runThreadN %s "
           "--limitGenomeGenerateRAM %s "
           "--genomeChrBinNbits %s "
           "--runMode genomeGenerate "
           "--genomeSAindexNbases %s" % (dir_name, str(cpu), str(mem), Nbases,
                                         nbits))
    if not os.path.exists(os.path.join(dir_name, "SA")):
        _index_w_command(env, dir_name, cmd, ref_file)
    return dir_name

@_if_installed("hisat2-build")
def _index_hisat2(env, ref_file):
    path_export = _get_path_export(env)
    build = os.path.splitext(os.path.basename(ref_file))[0]
    (ref_dir, local_file) = os.path.split(ref_file)
    gtf_file = os.path.join(ref_dir, os.pardir, "rnaseq", "ref-transcripts.gtf")
    dir_name = os.path.normpath(os.path.join(ref_dir, os.pardir, "hisat2"))
    index_prefix = os.path.join(dir_name, build)
    if os.path.exists(os.path.join(index_prefix + ".1.ht2")):
        return dir_name
    if not os.path.exists(dir_name):
        subprocess.check_call('mkdir -p %s' % dir_name, shell=True)
    try:
        cpu = env.cores
    except:
        cpu = 1
    cmd = "{path_export}hisat2-build -p {cpu} "

    exons_file = index_prefix + ".exons"
    splicesites_file = index_prefix + ".splicesites"
    if os.path.exists(gtf_file):
        if not os.path.exists(exons_file):
            with open(exons_file, "w") as out_handle:
                exons_cmd = ["hisat2_extract_exons.py", gtf_file]
                subprocess.check_call(path_export + " ".join(exons_cmd), stdout=out_handle, shell=True)
        if not os.path.exists(splicesites_file):
            with open(splicesites_file, "w") as out_handle:
                splicesites_cmd = ["hisat2_extract_splice_sites.py", gtf_file]
                subprocess.check_call(path_export + " ".join(splicesites_cmd), stdout=out_handle, shell=True)

        if os.stat(exons_file).st_size > 0 and os.stat(splicesites_file).st_size > 0:
            cmd += "--exon {exons_file} --ss {splicesites_file} "
    cmd += "{ref_file} {index_prefix} "
    if not os.path.exists(os.path.join(index_prefix + ".1.ht2")):
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return dir_name

def _index_snap(env, ref_file):
    """Snap indexing is computationally expensive. Requests all cores and 64Gb of memory.
    """
    dir_name = "snap"
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    org_arg = "-hg19" if index_name in ["hg19", "GRCh37"] else ""
    cmd = "snap-aligner index {ref_file} {dir_name} -bSpace {org_arg}"
    if not os.path.exists(os.path.join(dir_name, "GenomeIndex")):
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return dir_name

def _get_path_export(env):
    """Ensure PATH points to local install directory.
    """
    path_export = ""
    if hasattr(env, "system_install") and env.system_install:
        local_bin = os.path.join(env.system_install, 'bin')
        if os.path.exists(local_bin):
            path_export = "export PATH=%s:$PATH && " % local_bin
    return path_export

def _index_rtg(env, ref_file):
    """Perform indexing for use with Real Time Genomics tools.

    https://github.com/RealTimeGenomics/rtg-tools
    """
    path_export = _get_path_export(env)
    dir_name = "rtg"
    index_name = "%s.sdf" % os.path.splitext(os.path.basename(ref_file))[0]
    if not os.path.exists(os.path.join(dir_name, index_name, "done")):
        cmd = ("{path_export}export RTG_JAVA_OPTS='-Xms1g' && export RTG_MEM=2g && "
               "rtg format -o {dir_name}/{index_name} {ref_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return dir_name

@_if_installed("MosaikJump")
def _index_mosaik(env, ref_file):
    hash_size = 15
    dir_name = "mosaik"
    cmd = "MosaikBuild -fr {ref_file} -oa {index_name}"
    def create_jumpdb(ref_file):
        jmp_base = os.path.splitext(os.path.basename(ref_file))[0]
        dat_file = "{0}.dat".format(jmp_base)
        if not os.path.exists("{0}_keys.jmp".format(jmp_base)):
            cmd = "export MOSAIK_TMP=`pwd` && MosaikJump -hs {hash_size} -ia {ref_file} -out {index_name}".format(
                hash_size=hash_size, ref_file=dat_file, index_name=jmp_base)
            subprocess.check_call(cmd, shell=True)
    return _index_w_command(env, dir_name, cmd, ref_file,
                            post=create_jumpdb, ext=".dat")

# -- Retrieve using GGD recipes

def _install_with_ggd(env, manager, gid, recipe):
    recipe_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                               os.pardir, os.pardir, "ggd-recipes"))
    recipe_file = os.path.join(recipe_dir, gid, "%s.yaml" % recipe)
    if os.path.exists(recipe_file):
        ggd.install_recipe(os.getcwd(), env.system_install, recipe_file, gid)
    else:
        raise NotImplementedError("GGD recipe not available for %s %s" % (gid, recipe))

# -- Genome upload and download to Amazon s3 buckets

def _download_s3_index(env, manager, gid, idx):
    print("Downloading genome from s3: {0} {1}".format(gid, idx))
    url = "https://s3.amazonaws.com/biodata/genomes/%s-%s.tar.xz" % (gid, idx)
    if gid in ["GRCh37", "hg19", "mm10"] and idx in ["bowtie2", "bwa", "novoalign"]:
        out_file = shared._remote_fetch(env, url, samedir=True)
        subprocess.check_call("xz -dc %s | tar -xvpf -" % out_file, shell=True)
        subprocess.check_call("rm -f %s" % out_file, shell=True)
    else:
        raise NotImplementedError("No pre-computed indices for %s %s" % (gid, idx))

def _download_genomes(env, genomes, genome_indexes):
    """Download a group of genomes from Amazon s3 bucket.
    """
    genome_dir = _make_genome_dir(env.data_files)
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname, gid)
        if not os.path.exists(org_dir):
            subprocess.check_call('mkdir -p %s' % org_dir, shell=True)
        for idx in genome_indexes:
            with shared.chdir(org_dir):
                if not os.path.exists(idx):
                    _download_s3_index(env, manager, gid, idx)
        ref_file = os.path.join(org_dir, "seq", "%s.fa" % gid)
        if not os.path.exists(ref_file):
            ref_file = os.path.join(org_dir, "seq", "%s.fa" % manager._name)
        assert os.path.exists(ref_file), ref_file
        cur_indexes = manager.config.get("indexes", genome_indexes)
        _index_to_galaxy(env, org_dir, ref_file, gid, cur_indexes, manager.config)

def _upload_genomes(env, genomes, genome_indexes):
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
        gb_size = int(subprocess.check_output("du -sm %s" % tarball, shell=True).decode().split()[0]) / 1000.0
        print("Uploading %s %.1fGb" % (s3_key_name, gb_size))
        cl = ["python", upload_script, tarball, bucket.name, s3_key_name, "--public"]
        subprocess.check_call(cl)

def _tar_directory(dir, tar_name):
    """Create a tarball of the directory.
    """
    base_dir, tar_dir = os.path.split(dir)
    tarball = os.path.join(base_dir, "%s.tar.xz" % tar_name)
    if not os.path.exists(tarball):
        with shared.chdir(base_dir):
            subprocess.check_call("tar -cvpf - %s | xz -zc - > %s" %
                                  (tar_dir, os.path.basename(tarball)), shell=True)
    return tarball

def _clean_directory(dir, gid):
    """Clean duplicate files from directories before tar and upload.
    """
    # get rid of softlinks
    bowtie_ln = os.path.join(dir, "bowtie", "%s.fa" % gid)
    maq_ln = os.path.join(dir, "maq", "%s.fa" % gid)
    for to_remove in [bowtie_ln, maq_ln]:
        if os.path.exists(to_remove):
            subprocess.check_call("rm -f %s" % to_remove, shell=True)
    # remove any downloaded original sequence files
    remove_exts = ["*.gz", "*.zip"]
    with shared.chdir(os.path.join(dir, "seq")):
        for rext in remove_exts:
            fnames = subprocess.check_output("find . -name '%s'" % rext, shell=True).decode()
            for fname in (f.strip() for f in fnames.split("\n") if f.strip()):
                subprocess.check_call("rm -f %s" % fname, shell=True)

# == Liftover files

def _data_liftover(env, lift_over_genomes):
    """Download chain files for running liftOver.

    Does not install liftOver binaries automatically.
    """
    lo_dir = os.path.join(env.data_files, "liftOver")
    if not os.path.exists(lo_dir):
        subprocess.check_call("mkdir %s" % lo_dir, shell=True)
    lo_base_url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/liftOver/%s"
    lo_base_file = "%sTo%s.over.chain.gz"
    for g1 in lift_over_genomes:
        for g2 in [g for g in lift_over_genomes if g != g1]:
            g2u = g2[0].upper() + g2[1:]
            cur_file = lo_base_file % (g1, g2u)
            non_zip = os.path.splitext(cur_file)[0]
            worked = False
            with shared.chdir(lo_dir):
                if not os.path.exists(non_zip):
                    result = shared._remote_fetch(env, "%s" % (lo_base_url % (g1, cur_file)), allow_fail=True)
                    # Lift over back and forths don't always exist
                    # Only move forward if we found the file
                    if result:
                        worked = True
                        subprocess.check_call("gunzip %s" % result, shell=True)
            if worked:
                ref_parts = [g1, g2, os.path.join(lo_dir, non_zip)]
                galaxy.update_loc_file(env, "liftOver.loc", ref_parts)

# == UniRef
def _data_uniref(env):
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
        if not os.path.exists(work_dir):
            subprocess.check_call("mkdir -p %s" % work_dir, shell=True)
        base_work_url = base_url % (uniref_db, uniref_db)
        fasta_url = base_work_url + ".fasta.gz"
        base_file = os.path.splitext(os.path.basename(fasta_url))[0]
        with shared.chdir(work_dir):
            if not os.path.exists(base_file):
                out_file = shared._remote_fetch(env, fasta_url)
                subprocess.check_call("gunzip %s" % out_file, shell=True)
                shared._remote_fetch(env, base_work_url + ".release_note")
        _index_blast_db(work_dir, base_file, "prot")

def _index_blast_db(work_dir, base_file, db_type):
    """Index a database using blast+ for similary searching.
    """
    type_to_ext = dict(prot=("phr", "pal"), nucl=("nhr", "nal"))
    db_name = os.path.splitext(base_file)[0]
    with shared.chdir(work_dir):
        if not reduce(operator.or_,
                      (os.path.exists("%s.%s" % (db_name, ext)) for ext in type_to_ext[db_type])):
            subprocess.check_call("makeblastdb -in %s -dbtype %s -out %s" %
                                  (base_file, db_type, db_name), shell=True)


def get_index_fn(index):
    """
    return the index function for an index, if it is missing return a function
    that is a no-op
    """
    def noop(env, ref_file):
        pass
    return INDEX_FNS.get(index, noop)

INDEX_FNS = {
    "seq": _index_sam,
    "bbmap": _index_bbmap,
    "bismark": _index_bismark,
    "bwa": _index_bwa,
    "bowtie": _index_bowtie,
    "bowtie2": _index_bowtie2,
    "maq": _index_maq,
    "mosaik": _index_mosaik,
    "minimap2": _index_minimap2,
    "novoalign": _index_novoalign,
    "novoalign_cs": _index_novoalign_cs,
    "ucsc": _index_twobit,
    "twobit": _index_twobit,
    "star": _index_star,
    "snap": _index_snap,
    "rtg": _index_rtg,
    "hisat2": _index_hisat2
    }
