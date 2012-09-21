"""Handle alternative configuration file locations for flavor customizations.
"""
import os
import collections

def _find_fname(env, fname):
    for dirname in [env.get("flavor_dir", None), env.config_dir]:
        print dirname, fname
        if dirname:
            full_path = os.path.join(dirname, fname)
            if os.path.exists(full_path):
                return full_path
    return None

def get_config_file(env, fname):
    """Retrieve YAML configuration file from the default config directory or flavor directory.

    This combines all options for getting distribution or flavor specific customizations.
    """
    base, ext = os.path.splitext(fname)
    distribution_fname = "{0}-{1}{2}".format(base, env.get("distribution", "notspecified"), ext)
    Config = collections.namedtuple("Config", "base dist")
    out = Config(base=_find_fname(env, fname), dist=_find_fname(env, distribution_fname))
    return out
    
