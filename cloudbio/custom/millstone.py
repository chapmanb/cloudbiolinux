"""Install instructions for non-packaged programs required by Millstone.
"""

from fabric.api import cd

from cloudbio.custom.shared import _make_tmp_dir


def install_unafold(env):
    """Required by optmage.
    """
    # Since unafold is distributed as an .rpm, we need the program alien to
    # convert it into a .deb that can be installed on this system.
    env.safe_sudo("apt-get install -y alien")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            env.safe_run("wget http://dinamelt.rit.albany.edu/download/unafold-3.8-1.x86_64.rpm")
            env.safe_sudo("alien -i unafold-3.8-1.x86_64.rpm")
