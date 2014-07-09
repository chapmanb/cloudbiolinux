from cloudbio.galaxy.tools import _install_application


def install_tool(options):
    version = options.get("galaxy_tool_version")
    name = options.get("galaxy_tool_name")
    install_dir = options.get("galaxy_tool_dir", None)
    _install_application(name, version, tool_install_dir=install_dir)


configure_actions = {
    "install_galaxy_tool": install_tool,
}
