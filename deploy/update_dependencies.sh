#!/bin/bash

# Name of virtualenv to create using virtualenvwrapper 
VIRTUALENV_NAME=cbl_deploy

# Configure virtualenv and virtualenvwrapper with virtualenv-burrito 
. $HOME/.venvburrito/startup.sh

# Upgrade dependencies
mkvirtualenv -r requirements.txt $VIRTUALENV_NAME

