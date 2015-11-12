#!/usr/bin/env python
"""Map AZ300 gene name list to transcript regions.

This requires a 3 pass approach to get gene names to coordinates:

- it first tries against ensemble BED coordinates in a transcript file with
  the supplied name
- it then remaps the name with gene.info to alternative symbols and tries to
  find those in the transcript file.
- Finally, it uses coordinates from gene.info if those exist.
"""
import os
import sys

import requests

ref_dir = "/human"
in_file = sys.argv[1]

def read_targets(in_file):
    targets = set([])
    with open(in_file) as in_handle:
        for line in in_handle:
            cur_symbol = (line.strip())
            targets.add(cur_symbol)
    return targets

def get_gene_info(cur_symbol):
    chroms = [str(x) for x in range(1, 23)] + ["X", "Y"]
    fields = "symbol,genomic_pos_hg19"
    url = "http://mygene.info/v2/query?q=%s&species=human&fields=%s" % (cur_symbol, fields)
    info = requests.get(url).json()
    hits = info["hits"]
    symbols = [x["symbol"] for x in hits]
    pos = []
    for ps in [x["genomic_pos_hg19"] for x in hits if "genomic_pos_hg19" in x]:
        if not isinstance(ps, (list, tuple)):
            ps = [ps]
        for p in ps:
            if p["chr"] in chroms:
                pos.append(p)
    return symbols, pos

def find_missing_targets(missing, in_file, genome):
    missing_file = "%s-%s-missingsymbols.txt" % (os.path.splitext(in_file)[0], genome)
    out = set([])
    if os.path.exists(missing_file):
        with open(missing_file) as in_handle:
            for line in in_handle:
                out.add(line.strip())
    else:
        with open(missing_file, "w") as out_handle:
            for cur_symbol in missing:
                symbols, pos = get_gene_info(cur_symbol)
                if cur_symbol not in symbols and len(symbols) > 0:
                    cur_symbol = symbols[0]
                out_handle.write("%s\n" % cur_symbol)
                out.add(cur_symbol)
    return out

def write_from_transcript_file(targets, ref_dir, genome, out_handle):
    ref_file = os.path.join(ref_dir, genome, "rnaseq", "ref-transcripts.bed")
    found = set([])
    with open(ref_file) as in_handle:
        for line in in_handle:
            name = line.split()[3]
            if name in targets:
                found.add(name)
                out_handle.write(line)
    return targets - found

def write_from_remap_names(targets, ref_dir, genome, out_handle, in_file):
    ref_file = os.path.join(ref_dir, genome, "rnaseq", "ref-transcripts.bed")
    targets = find_missing_targets(targets, in_file, genome)
    found = set([])
    with open(ref_file) as in_handle:
        for line in in_handle:
            name = line.split()[3]
            if name in targets:
                found.add(name)
                out_handle.write(line)
    return targets - found

def write_from_gene_info(targets, genome, out_handle):
    missing = []
    for target in sorted(list(read_targets(in_file))):
        symbols, pos = get_gene_info(target)
        if pos:
            assert isinstance(pos, (list, tuple))
            if symbols:
                target = symbols[0]
            for p in pos:
                chrom = "%s%s" % ("chr" if genome == "hg19" else "", p["chr"])
                out_handle.write("\t".join([chrom, str(p["start"]), str(p["end"]), target]) + "\n")
        else:
            missing.append(target)
    return missing

for genome in ["hg19", "GRCh37"]:
    out_file = "%s-%s.bed" % (os.path.splitext(in_file)[0], genome)
    with open(out_file, "w") as out_handle:
        targets = read_targets(in_file)
        print "total", len(targets)
        targets = write_from_transcript_file(targets, ref_dir, genome, out_handle)
        print("after first name pass", len(targets))
        targets = write_from_remap_names(targets, ref_dir, genome, out_handle, in_file)
        print("after rename name pass", len(targets))
        targets = write_from_gene_info(targets, genome, out_handle)
        print("after coordinate retrieval", targets)