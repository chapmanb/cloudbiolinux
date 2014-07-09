import os
import tempfile


from cloudbio.custom.shared import _make_tmp_dir
from fabric.api import sudo, put, cd

r_packages_template = """
r <- getOption("repos");
r["CRAN"] <- "http://watson.nci.nih.gov/cran_mirror";
options(repos=r);
install.packages( c( %s ), dependencies = TRUE);
source("http://bioconductor.org/biocLite.R");
biocLite( c( %s ) );
"""


def _install_r_packages(tools_conf):
    f = tempfile.NamedTemporaryFile()
    r_packages = tools_conf["r_packages"]
    bioconductor_packages = tools_conf["bioconductor_packages"]
    if not r_packages and not bioconductor_packages:
        return
    r_cmd = r_packages_template % (_concat_strings(r_packages), _concat_strings(bioconductor_packages))
    f.write(r_cmd)
    f.flush()
    with _make_tmp_dir() as work_dir:
        put(f.name, os.path.join(work_dir, 'install_packages.r'))
        with cd(work_dir):
            sudo("R --vanilla --slave < install_packages.r")
    f.close()


def _concat_strings(strings):
    if strings:
        return ", ".join(map(lambda x: '"%s"' % x, strings))
    else:
        return ""
