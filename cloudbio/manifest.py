"""Provide dump of software and libraries installed on CloudBioLinux image.

This provides an output YAML file with package details, providing a complete
dump of installed software and packages. The YAML output feeds into a BioGems
style webpage that provides a more human friendly view of installed packages.
The version information provides a reproducible dump of software on a system.
"""
import os
import collections
import inspect
import urllib2
import subprocess
import sys

import yaml
try:
    import yolk.yolklib
    import yolk.metadata
except ImportError:
    yolk = None

def create(out_dir, tooldir="/usr/local", fetch_remote=False):
    """Create a manifest in the output directory with installed packages.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    write_debian_pkg_info(out_dir, fetch_remote)
    write_python_pkg_info(out_dir)
    write_r_pkg_info(out_dir)
    write_brew_pkg_info(out_dir, tooldir)
    write_custom_pkg_info(out_dir, tooldir)

# ## Custom packages

def _get_custom_pkg_info(name, fn):
    """Retrieve information about the installed package from the install function.
    """
    vals = dict((k, v) for k, v in inspect.getmembers(fn))
    code = inspect.getsourcelines(fn)
    if vals["__name__"] == "decorator":
        fn = [x for x in fn.func_closure if not isinstance(x.cell_contents, str)][0].cell_contents
        vals = dict((k, v) for k, v in inspect.getmembers(fn))
        code = inspect.getsourcelines(fn)
    version = ""
    for line in (l.strip() for l in code[0]):
        if line.find("version") >= 0 and line.find(" =") > 0:
            version = line.split()[-1].replace('"', '').replace("'", "")
        if version:
            break
    doc = vals.get("func_doc", "")
    descr, homepage = "", ""
    if doc is not None:
        descr = doc.split("\n")[0]
        for line in doc.split("\n"):
            if line.strip().startswith("http"):
                homepage = line.strip()
    return {"name": name.replace("install_", ""),
            "description": descr,
            "homepage_uri": homepage,
            "version": version}

def write_custom_pkg_info(out_dir, tooldir):
    custom_names = ["bio_general", "bio_nextgen", "cloudman", "distributed",
                    "java", "python", "phylogeny", "system"]
    out_file = os.path.join(out_dir, "custom-packages.yaml")
    if not os.path.exists(out_file):
        out = {}
        for modname in custom_names:
            mod = getattr(__import__("cloudbio.custom", globals(), locals(),
                                     [modname], -1),
                          modname)
            for prog in [x for x in dir(mod) if x.startswith("install")]:
                pkg = _get_custom_pkg_info(prog, getattr(mod, prog))
                out[pkg["name"]] = pkg
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

# ## Homebrew/Linuxbrew packages

def write_brew_pkg_info(out_dir, tooldir):
    """Extract information for packages installed by homebrew/linuxbrew.
    """
    out_file = os.path.join(out_dir, "brew-packages.yaml")
    if not os.path.exists(out_file):
        brew_cmd = os.path.join(tooldir, "bin", "brew") if tooldir else None
        if not brew_cmd or not os.path.exists(brew_cmd):
            brew_cmd = "brew"
        try:
            vout = subprocess.check_output([brew_cmd, "list", "--versions"])
        except OSError:  # brew not installed/used
            vout = ""
        out = {}
        for vstr in vout.split("\n"):
            if vstr.strip():
                parts = vstr.rstrip().split()
                name = parts[0]
                v = parts[-1]
                out[name] = {"name": name, "version": v}
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

# ## R packages

def get_r_pkg_info():
    r_command = ("options(width=10000); subset(installed.packages(fields=c('Title', 'URL')), "
                 "select=c('Version', 'Title','URL'))")
    try:
        out = subprocess.check_output(["Rscript", "-e", r_command])
    except (subprocess.CalledProcessError, OSError):
        out = ""
    pkg_raw_list = []
    for line in out.split("\n")[1:]:
        pkg_raw_list.append(filter(None, [entry.strip(' ') for entry in line.split('"')]))
    for pkg in pkg_raw_list:
        if len(pkg) > 2:
            yield {"name": pkg[0], "version": pkg[1],
                   "description": pkg[2],
                   "homepage_uri": (pkg[3], '')[pkg[3] == 'NA'] if len(pkg) > 3 else ""}

def write_r_pkg_info(out_dir):
    out_file = os.path.join(out_dir, "r-packages.yaml")
    if not os.path.exists(out_file):
        out = {}
        for pkg in get_r_pkg_info():
            out[pkg["name"]] = pkg
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

# ## Python packages

def get_python_pkg_info():
    if yolk:
        for dist in yolk.yolklib.Distributions().get_packages("all"):
            md = yolk.metadata.get_metadata(dist)
            yield {"name": md["Name"].lower(), "version": md["Version"],
                   "description": md.get("Summary", ""),
                   "homepage_uri": md.get("Home-page", "")}
    else:
        base_dir = os.path.dirname(sys.executable)
        if os.path.exists(os.path.join(base_dir, "conda")):
            for line in subprocess.check_output([os.path.join(base_dir, "conda"), "list"]).split("\n"):
                if line.strip() and not line.startswith("#"):
                    name, version = line.split()[:2]
                    yield {"name": name.lower(), "version": version}
        else:
            for line in subprocess.check_output([os.path.join(base_dir, "pip"), "list"]).split("\n"):
                if line.strip() and not line.startswith("#"):
                    name, version = line.split()[:2]
                    yield {"name": name.lower(), "version": version[1:-1]}

def _resolve_latest_pkg(pkgs):
    if len(pkgs) == 1 or not yolk:
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
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

# ## Debian packages

def _get_pkg_popcon():
    """Retrieve popularity information for debian packages.
    """
    url = "http://popcon.debian.org/by_vote"
    popcon = {}
    for line in (l for l in urllib2.urlopen(url) if not l.startswith(("#", "--"))):
        parts = line.split()
        popcon[parts[1]] = int(parts[3])
    return popcon

def get_debian_pkg_info(fetch_remote=False):
    pkg_popcon = _get_pkg_popcon() if fetch_remote else {}
    cmd = ("dpkg-query --show --showformat "
           "'${Status}\t${Package}\t${Version}\t${Section}\t${Homepage}\t${binary:Summary}\n'")
    for pkg_line in [l for l in subprocess.check_output(cmd, shell=True).split("\n")
                     if l.startswith("install ok")]:
        parts = pkg_line.rstrip("\n").split("\t")
        if len(parts) > 5:
            pkg = {"name": parts[1], "version": parts[2],
                   "section": parts[3], "homepage_uri": parts[4],
                   "description": parts[5]}
            if pkg_popcon.get(pkg["name"]):
                pkg["downloads"] = pkg_popcon.get(pkg["name"], 0)
            yield pkg

def write_debian_pkg_info(out_dir, fetch_remote=False):
    base_sections = set(["gnome", "admin", "utils", "web", "games",
                         "sound", "devel", "kde", "x11", "net", "text",
                         "graphics", "misc", "editors", "fonts", "doc",
                         "mail", "otherosfs", "video", "kernel",
                         "libs", "libdevel", "comm", "metapackages", "tex",
                         "introspection"])
    for s in list(base_sections):
        base_sections.add("universe/%s" % s)
        base_sections.add("partner/%s" % s)
    out_file = os.path.join(out_dir, "debian-packages.yaml")
    out_base_file = os.path.join(out_dir, "debian-base-packages.yaml")
    try:
        subprocess.check_call(["dpkg", "--help"], stdout=subprocess.PIPE)
        has_dpkg = True
    except (subprocess.CalledProcessError, OSError):
        has_dpkg = False
    if has_dpkg and (not os.path.exists(out_file) or not os.path.exists(out_base_file)):
        out = {}
        out_base = {}
        for pkg in get_debian_pkg_info(fetch_remote):
            if pkg.get("section") in base_sections:
                out_base[pkg["name"]] = pkg
            else:
                out[pkg["name"]] = pkg
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
        with open(out_base_file, "w") as out_handle:
            yaml.safe_dump(out_base, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file
