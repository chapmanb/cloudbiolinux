"""Retrieve indexed genomes using Galaxy's rsync server resources.

http://wiki.galaxyproject.org/Admin/Data%20Integration
"""
from xml.etree import ElementTree

from fabric.api import *
from fabric.contrib.files import *

# ## Compatibility definitions

server = "rsync://datacache.g2.bx.psu.edu"

index_map = {"bowtie": "bowtie_index",
             "bowtie2": "bowtie2_index",
             "bwa": "bwa_index",
             "novoalign": None,
             "ucsc": "seq",
             "seq": "sam_index"}

org_remap = {"phix": "phiX",
             "GRCh37": "hg_g1k_v37",
             "araTha_tair9": "Arabidopsis_thaliana_TAIR9",
             "araTha_tair10": "Arabidopsis_thaliana_TAIR10",
             "WS210": "ce10",
             "WS220": "ce10"}

galaxy_subdirs = ["", "/microbes"]

# ## Galaxy location files

class LocCols(object):
    # Hold all possible .loc file column fields making sure the local
    # variable names match column names in Galaxy's tool_data_table_conf.xml
    def __init__(self, config, dbkey, file_path):
        self.dbkey = dbkey
        self.path = file_path
        self.value = config.get("value", dbkey)
        self.name = config.get("name", dbkey)
        self.species = config.get('species', '')
        self.index = config.get('index', 'index')
        self.formats = config.get('index', 'fastqsanger')
        self.dbkey1 = config.get('index', dbkey)
        self.dbkey2 = config.get('index', dbkey)

def _get_tool_conf(tool_name):
    """
    Parse the tool_data_table_conf.xml from installed_files subfolder and extract
    values for the 'columns' tag and 'path' parameter for the 'file' tag, returning
    those as a dict.
    """
    tool_conf = {}
    tdtc = ElementTree.parse(env.tool_data_table_conf_file)
    tables = tdtc.getiterator('table')
    for t in tables:
        if tool_name in t.attrib.get('name', ''):
            tool_conf['columns'] = t.find('columns').text.replace(' ', '').split(',')
            tool_conf['file'] = t.find('file').attrib.get('path', '')
    return tool_conf

def _build_galaxy_loc_line(dbkey, file_path, config, prefix, tool_name):
    """Prepare genome information to write to a Galaxy *.loc config file.
    """
    if tool_name:
        str_parts = []
        tool_conf = _get_tool_conf(tool_name)
        loc_cols = LocCols(config, dbkey, file_path)
        # Compose the .loc file line as str_parts list by looking for column values
        # from the retrieved tool_conf (as defined in tool_data_table_conf.xml).
        # Any column values required but missing in the tool_conf are
        # supplemented by the defaults defined in LocCols class
        for col in tool_conf.get('columns', []):
            str_parts.append(config.get(col, getattr(loc_cols, col)))
    else:
        str_parts = [dbkey, file_path]
    if prefix:
        str_parts.insert(0, prefix)
    return str_parts

def update_loc_file(ref_file, line_parts):
    """Add a reference to the given genome to the base index file.
    """
    if getattr(env, "galaxy_home", None) is not None:
        tools_dir = os.path.join(env.galaxy_home, "tool-data")
        if not exists(tools_dir):
            run("mkdir -p %s" % tools_dir)
            put(env.tool_data_table_conf_file,
                os.path.join(env.galaxy_home, "tool_data_table_conf.xml"))
        add_str = "\t".join(line_parts)
        with cd(tools_dir):
            if not exists(ref_file):
                run("touch %s" % ref_file)
            if not contains(ref_file, add_str):
                append(ref_file, add_str)

def prep_locs(gid, indexes, config):
    """Prepare Galaxy location files for all available indexes.
    """
    for ref_index_file, cur_index, prefix, tool_name in [
          ("sam_fa_indices.loc", indexes.get("seq", None), "", 'sam_fa_indexes'),
          ("picard_index.loc", indexes.get("seq", None), "", "picard_indexes"),
          ("gatk_sorted_picard_index.loc", indexes.get("seq", None), "", "gatk_picard_indexes"),
          ("alignseq.loc", indexes.get("ucsc", None), "seq", None),
          ("twobit.loc", indexes.get("ucsc", None), "", None),
          ("bowtie_indices.loc", indexes.get("bowtie", None), "", 'bowtie_indexes'),
          ("bowtie2_indices.loc", indexes.get("bowtie2", None), "", 'bowtie2_indexes'),
          ("mosaik_index.loc", indexes.get("mosaik", None), "", "mosaik_indexes"),
          ("bwa_index.loc", indexes.get("bwa", None), "", 'bwa_indexes')]:
        if cur_index:
            str_parts = _build_galaxy_loc_line(gid, cur_index, config, prefix, tool_name)
            update_loc_file(ref_index_file, str_parts)

# ## Finalize downloads

def index_picard(ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    index_file = "%s.dict" % os.path.splitext(ref_file)[0]
    dirs_to_try = ["%s/share/java/picard" % env.system_install,
                   getattr(env, "picard_home", None)]
    picard_jar = None
    for dname in dirs_to_try:
        if dname:
            test_jar = os.path.join(dname, "CreateSequenceDictionary.jar")
            if exists(test_jar):
                picard_jar = test_jar
                break
    if picard_jar and not exists(index_file):
        run("java -jar {jar} REFERENCE={ref} OUTPUT={out}".format(
            jar=picard_jar, ref=ref_file, out=index_file))
    return index_file

def _finalize_index_seq(fname):
    """Convert UCSC 2bit file into fasta file.
    """
    out_fasta = fname + ".fa"
    if not exists(out_fasta):
        run("twoBitToFa {base}.2bit {out}".format(
            base=fname, out=out_fasta))

finalize_fns = {"ucsc": _finalize_index_seq,
                "seq": index_picard}

def _finalize_index(idx, fname):
    """Perform final processing on an rsync'ed index file if necessary.
    """
    finalize_fn = finalize_fns.get(idx)
    if finalize_fn:
        finalize_fn(fname)

# ## Retrieve data from Galaxy

def rsync_genomes(genome_dir, genomes, genome_indexes):
    """Top level entry point to retrieve rsync'ed indexes from Galaxy.
    """
    for gid in (x[1] for x in genomes):
        galaxy_gid = org_remap.get(gid, gid)
        indexes = _get_galaxy_genomes(galaxy_gid, genome_dir, genomes, genome_indexes)
        _finalize_index("ucsc", indexes["ucsc"])
        for idx, fname in indexes.iteritems():
            _finalize_index(idx, fname)
        prep_locs(galaxy_gid, indexes, {})

def _get_galaxy_genomes(gid, genome_dir, genomes, genome_indexes):
    """Retrieve the provided genomes and indexes from Galaxy rsync.
    """
    out = {}
    org_dir = os.path.join(genome_dir, gid)
    if not exists(org_dir):
        run('mkdir -p %s' % org_dir)
    for idx in genome_indexes:
        galaxy_index_name = index_map.get(idx)
        index_file = None
        if galaxy_index_name:
            index_file = _rsync_genome_index(gid, galaxy_index_name, org_dir)
        if index_file:
            out[idx] = index_file
        else:
            print "Galaxy does not support {0} for {1}".format(idx, gid)
    return out

def _rsync_genome_index(gid, idx, org_dir):
    """Retrieve index for a genome from rsync server, returning path to files.
    """
    idx_dir = os.path.join(org_dir, idx)
    if not exists(idx_dir):
        org_rsync = None
        for subdir in galaxy_subdirs:
            test_rsync = "{server}/indexes{subdir}/{gid}/{idx}/".format(
                server=server, subdir=subdir, gid=gid, idx=idx)
            with quiet():
                check_dir = run("rsync --list-only {server}".format(server=test_rsync))
            if check_dir.succeeded:
                org_rsync = test_rsync
                break
        if org_rsync is None:
            raise ValueError("Could not find genome %s on Galaxy rsync" % gid)
        with quiet():
            check_dir = run("rsync --list-only {server}".format(server=org_rsync))
        if check_dir.succeeded:
            if not exists(idx_dir):
                run('mkdir -p %s' % idx_dir)
            with cd(idx_dir):
                run("rsync -avzP {server} {idx_dir}".format(server=org_rsync,
                                                            idx_dir=idx_dir))
    if exists(idx_dir):
        with quiet():
            has_fa_ext = run("ls {idx_dir}/{gid}.fa*".format(idx_dir=idx_dir,
                                                             gid=gid))
        ext = ".fa" if (has_fa_ext.succeeded and idx not in ["seq"]) else ""
        return os.path.join(idx_dir, gid + ext)
