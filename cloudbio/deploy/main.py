from argparse import ArgumentParser
import yaml

from cloudbio.deploy import deploy

DESC = "Creates an on-demand cloud instance, sets up applications, and transfer files to it."

## Properties that may be specified as args or in settings file,
## argument takes precedence.
ARG_PROPERTIES = [
  # VM launcher options
  "files",
  "compressed_files",
  "actions",
  "runtime_properties",
  "vm_provider",
  "hostname",

  # CloudBioLinux options
  "target",
  "flavor",
  "package",

  # CloudMan options
  "target_bucket",
]


def main():
    args = parse_args()
    options = parse_settings(args.settings)

    for property in ARG_PROPERTIES:
        _copy_arg_to_options(options, args, property)

    deploy(options)


def _copy_arg_to_options(options, args, property):
    arg_property = getattr(args, property)
    if arg_property or not property in options:
        options[property] = arg_property


def parse_args():
    parser = ArgumentParser(DESC)
    parser.add_argument("--settings", dest="settings", default="settings.yaml")
    parser.add_argument('--action', dest="actions", action="append", default=[])
    parser.add_argument('--runtime_property', dest="runtime_properties", action="append", default=[])
    parser.add_argument('--compressed_file', dest="compressed_files", action="append", default=[], help="file to transfer to new instance and decompress")
    parser.add_argument('--file', dest="files", action="append", default=[], help="file to transfer to new instance")
    parser.add_argument("--vm_provider", dest="vm_provider", default=None, help="libcloud driver to use (or vagrant) (e.g. aws, openstack)")
    parser.add_argument("--hostname", dest="hostname", default=None, help="Newly created nodes are created with this specified hostname.")

    # CloudBioLinux options
    parser.add_argument("--target", dest="target", default=None, help="Specify a CloudBioLinux target, used with action install_biolinux action")
    parser.add_argument("--flavor", dest="flavor", default=None, help="Specify a CloudBioLinux flavor, used with action install_biolinux action")
    parser.add_argument("--package", dest="package", default=None, help="Specify a CloudBioLinux package, used with action install_custom")

    # CloudMan related options
    parser.add_argument("--target_bucket", dest="target_bucket", default=None, help="Specify a target bucket for CloudMan bucket related actions.")

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
