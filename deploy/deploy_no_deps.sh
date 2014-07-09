#!/bin/sh

export PROJECT_DIRECTORY="."

# Add cloudbiolinux to python path and run deployment.
export PYTHONPATH=..:$PYTHONPATH
python $PROJECT_DIRECTORY/../cloudbio/deploy/main.py "$@"
