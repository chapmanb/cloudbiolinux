#!/bin/bash

project_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $project_directory
if [ ! -e .venv-deploy ];
then
    ./setup.sh
fi

export PYTHONPATH=..:$PYTHONPATH
tools/with_venv.sh python $project_directory/../cloudbio/deploy/main.py "$@"
