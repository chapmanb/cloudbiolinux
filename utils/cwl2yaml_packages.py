#!/usr/bin/env python
# Parses CWL's "SoftwareRequirement" hints section and dumps a cbl compatible yaml file.
# The purpose with this script is to create smaller composable docker containers for bcbio-nextgen.
#
# Usage: cwl2yaml_packages.py test_bcbio_cwl/run_info-cwl-workflow/steps/process_alignment.cwl > cloudbiolinux/contrib/flavor/cwl_dockers/packages-bcbio-alignment.yaml
import os
import sys
import yaml

CWL_STEPS=sys.argv[1]
cwl_pkgs=yaml.safe_load(open(CWL_STEPS,'r'))
cbl_yml=dict()
cbl_pkgs=[]

# take the filename as the flavor/dockerfile name
cbl_flavor="bcbio-"+os.path.splitext(os.path.basename(sys.argv[1]))[0]

for pkg in cwl_pkgs['hints'][1]['packages']:
    cbl_pkgs.append(pkg['package'])

cbl_yml['channels']=['bioconda', 'conda-forge']
cbl_yml[cbl_flavor]=cbl_pkgs

#print cbl_yml

print yaml.safe_dump(cbl_yml, default_flow_style=False, indent=4)