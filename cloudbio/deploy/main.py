from argparse import ArgumentParser
import yaml

from cloudbio.deploy import deploy

DESC = "Creates an on-demand cloud instance, sets up applications, and transfer files to it."


def main():
    args = parse_args()
    options = parse_settings(args.settings)
    options["files"] = args.files
    options["compressed_files"] = args.compressed_files
    options["actions"] = args.actions
    options["runtime_properties"] = args.runtime_properties
    deploy(options)


def parse_args():
    parser = ArgumentParser(DESC)
    parser.add_argument("--settings", dest="settings", default="settings.yaml")
    parser.add_argument('--action', dest="actions", action="append", default=[])
    parser.add_argument('--runtime_property', dest="runtime_properties", action="append", default=[])
    parser.add_argument('--compressed_file', dest="compressed_files", action="append", default=[], help="file to transfer to new instance and decompress")
    parser.add_argument('--file', dest="files", action="append", default=[], help="file to transfer to new instance")
    args = parser.parse_args()
    if len(args.actions) == 0:
        args.actions = ["transfer"]
    return args


def parse_settings(name):
    return _read_yaml(name)


def _read_yaml(yaml_file):
    with open(yaml_file) as in_handle:
        return yaml.load(in_handle)


if __name__ == "__main__":
    main()
