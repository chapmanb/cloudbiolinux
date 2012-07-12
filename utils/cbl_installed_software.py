#!/usr/bin/env python
"""Provide dump of software and libraries installed on CloudBioLinux image.

This provides an output YAML file with package details, providing a complete
dump of installed software and packages. The YAML output feeds into a BioGems
style webpage that provides a more human friendly view of installed packages.
The version information provides a reproducible dump of software on a system.
"""
import os
import collections
import urllib2
import subprocess

import yaml
import yolk.yolklib, yolk.metadata

def main():
    out_dir = os.path.join(os.getcwd(), "manifest")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    write_debian_pkg_info(out_dir)
    write_python_pkg_info(out_dir)
    write_r_pkg_info(out_dir)

# ## R packages

def get_r_pkg_info():
    r_command = "options(width=10000); subset(installed.packages(fields=c('Title', 'URL')), select=c('Version', 'Title','URL'))"
    out = subprocess.check_output(["Rscript", "-e",r_command])
    pkg_raw_list = []
    for line in out.split("\n")[1:]:
	pkg_raw_list.append(filter(None, [entry.strip(' ') for entry in line.split('"')]))
    for pkg in pkg_raw_list:
	if len(pkg)>0:
            yield {"name": pkg[0], "version": pkg[1],
                   "description": pkg[2],
                   "homepage_uri": (pkg[3],'')[pkg[3]=='NA']}

def write_r_pkg_info(out_dir):
    out_file = os.path.join(out_dir, "r-packages.yaml")
    if not os.path.exists(out_file):
        out = {}
        for pkg in get_r_pkg_info():
            out[pkg["name"]] = pkg
        with open(out_file, "w") as out_handle:
            out_handle.write(yaml.dump(out, default_flow_style=False, allow_unicode=False))
    return out_file

# ## Python packages

def get_python_pkg_info():
    for dist in yolk.yolklib.Distributions().get_packages("all"):
        md = yolk.metadata.get_metadata(dist)
        yield {"name": md["Name"], "version": md["Version"],
               "description": md.get("Summary", ""),
               "homepage_uri": md.get("Home-page", "")}

def _resolve_latest_pkg(pkgs):
    if len(pkgs) == 1:
        return pkgs[0]
    else:
        latest_version = yolk.yolklib.Distributions().get_highest_installed(pkgs[0]["name"])
        return [x for x in pkgs if x["version"] == latest_version][0]

def write_python_pkg_info(out_dir):
    out_file = os.path.join(out_dir, "python-packages.yaml")
    if not os.path.exists(out_file):
        pkgs_by_name = collections.defaultdict(list)
        for pkg in get_python_pkg_info():
            pkgs_by_name[pkg["name"]].append(pkg)
        out = {}
        for name in sorted(pkgs_by_name.keys()):
            out[name] = _resolve_latest_pkg(pkgs_by_name[name])
        with open(out_file, "w") as out_handle:
            out_handle.write(yaml.dump(out, default_flow_style=False, allow_unicode=False))
    return out_file

# ## Debian packages

def _parse_pkg_info(pkg):
    details = {}
    for line in subprocess.check_output(["apt-cache", "show", pkg["name"]]).split("\n"):
        if ":" in line:
            name, info = line.split(":", 1)
            details[name] = info
    return {"homepage_uri": details.get("Homepage", "").strip(),
            "section": details.get("Section", "").strip()}

def _get_pkg_popcon():
    """Retrieve popularity information for debian packages.
    """
    url = "http://popcon.debian.org/by_vote"
    popcon = {}
    for line in (l for l in urllib2.urlopen(url) if not l.startswith(("#", "--"))):
        parts = line.split()
        popcon[parts[1]] = int(parts[3])
    return popcon

def get_debian_pkg_info():
    pkg_popcon = _get_pkg_popcon()
    for pkg_line in [l for l in subprocess.check_output(["dpkg", "-l"]).split("\n")
                     if l.startswith("ii")]:
        info = pkg_line.split()[1:]
        pkg = {"name": info[0], "version": info[1], "description": " ".join(info[2:]),
               "downloads": pkg_popcon.get(info[0], 0)}
        pkg.update(_parse_pkg_info(pkg))
        yield pkg

def write_debian_pkg_info(out_dir):
    out_file = os.path.join(out_dir, "debian-packages.yaml")
    if not os.path.exists(out_file):
        out = {}
        for pkg in get_debian_pkg_info():
            out[pkg["name"]] = pkg
        with open(out_file, "w") as out_handle:
            out_handle.write(yaml.dump(out, default_flow_style=False, allow_unicode=False))
    return out_file

if __name__ == "__main__":
    main()
