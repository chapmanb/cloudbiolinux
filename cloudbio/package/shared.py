"""Shared functionality useful for multiple package managers.
"""
from __future__ import print_function
import yaml

def _yaml_to_packages(yaml_file, to_install=None, subs_yaml_file=None, namesort=True, env=None):
    """Read a list of packages from a nested YAML configuration file.
    """
    print("Reading packages from %s" % yaml_file)
    with open(yaml_file) as in_handle:
        full_data = yaml.safe_load(in_handle)
        if full_data is None:
            full_data = {}
    if subs_yaml_file is not None:
        with open(subs_yaml_file) as in_handle:
            subs = yaml.safe_load(in_handle)
    else:
        subs = {}
    # filter the data based on what we have configured to install
    data = [(k, v) for (k, v) in full_data.items()
            if (to_install is None or k in to_install) and k not in ["channels"]]
    data.sort()
    packages = []
    pkg_to_group = dict()
    while len(data) > 0:
        cur_key, cur_info = data.pop(0)
        if cur_info:
            if isinstance(cur_info, (list, tuple)):
                packages.extend(_filter_subs_packages(cur_info, subs, namesort))
                for p in cur_info:
                    pkg_to_group[p] = cur_key
            elif isinstance(cur_info, dict):
                for key, val in cur_info.items():
                    # if we are okay, propagate with the top level key
                    if env and key == 'needs_64bit':
                        if env.is_64bit:
                            data.insert(0, (cur_key, val))
                    elif env and key.startswith(env.distribution):
                        if key.endswith(env.dist_name):
                            data.insert(0, (cur_key, val))
                    else:
                        data.insert(0, (cur_key, val))
            else:
                raise ValueError(cur_info)
    return packages, pkg_to_group

def _filter_subs_packages(initial, subs, namesort=True):
    """Rename and filter package list with subsitutions; for similar systems.
    """
    final = []
    for p in initial:
        try:
            new_p = subs[p]
        except KeyError:
            new_p = p
        if new_p:
            final.append(new_p)
    if namesort:
        final.sort()
    return final
