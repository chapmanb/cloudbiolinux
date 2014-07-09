import inspect
import os
import yaml


def parse_settings(name="deploy/settings.yaml"):
    return _read_yaml(_path_from_root(name))


def _path_from_root(name):
    root_path = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), "..", "..")
    file_path = os.path.join(root_path, name)
    return file_path


def _read_yaml(yaml_file):
    with open(yaml_file) as in_handle:
        return yaml.load(in_handle)
