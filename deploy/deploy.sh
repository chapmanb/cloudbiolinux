#!/bin/bash

# Name of virtualenv to create using virtualenvwrapper
VIRTUALENV_NAME=cbl_deploy

# Ensure working directory is cloudbiolinux/deploy. 
PROJECT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $PROJECT_DIRECTORY

# Ensure virtualenv-burrito has been installed.
if [ ! -e $HOME/.venvburrito/startup.sh ];
then
    wget -qO- https://raw.github.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh | $SHELL
fi

# Configure virtualenv and virtualenvwrapper with virtualenv-burrito
. $HOME/.venvburrito/startup.sh

# If no cbl_deploy virtualenv exists, create it and populate via
# pip-requires
if [ ! `lsvirtualenv | grep $VIRTUALENV_NAME` ];
then
    mkvirtualenv -r requirements.txt $VIRTUALENV_NAME
fi

# Use cbl_deploy virtualenv
workon $VIRTUALENV_NAME

# Add cloudbiolinux to python path and run deployment.
export PYTHONPATH=..:$PYTHONPATH
python $PROJECT_DIRECTORY/../cloudbio/deploy/main.py "$@"
