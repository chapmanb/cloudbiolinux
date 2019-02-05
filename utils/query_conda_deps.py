#!/usr/bin/env python
"""Query conda dependencies, summarizing for environment stratification.

Helps identify grouping of packages for different sub-environments.
"""
import sys

import requests
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    channels = config["channels"]
    channels.reverse()
    for p in sorted(config["bio_nextgen"]):
        p = p.split(";")[0]
        env = "default"
        for c in channels:
            deps = get_dependencies(c, p)
            if deps is not None:
                if "python" in deps and deps.get("python"):
                    if not any([s[1].startswith("3.") for s in deps["python"]]):
                        env = "python2"
                break
        print(p, env)

def get_dependencies(channel, package):
    url = "https://api.anaconda.org/package/%s/%s/files" % (channel, package)
    out = requests.get(url).json()
    if isinstance(out, dict) and out.get("error"):
        return None
    else:
        deps = {"python": set([])}
        for f in out:
            for d in f["dependencies"]["depends"]:
                for k in deps.keys():
                    if d["name"] == k:
                        for s in d["specs"]:
                            deps[k].add(tuple(s))
        return deps

if __name__ == "__main__":
    main(sys.argv[1])
