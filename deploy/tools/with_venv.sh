#!/bin/bash
TOOLS=`dirname $0`
VENV=$TOOLS/../.venv-deploy
source $VENV/bin/activate && $@
