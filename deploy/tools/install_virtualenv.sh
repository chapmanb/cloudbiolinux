#!/bin/bash
export LOCAL_PYTHON=$HOME/.gvl_python
export PYTHON_VERSION=${PYTHON_VERSION:-`python -c "import sys; rev = sys.version_info; str = '%d.%d' % (rev[0], rev[1]); print str"`}
mkdir -p $LOCAL_PYTHON/lib/python$PYTHON_VERSION/site-packages 
export PYTHONPATH=$LOCAL_PYTHON/lib/python$PYTHON_VERSION/site-packages:$PYTHONPATH
export PATH=$LOCAL_PYTHON/bin:$PATH
easy_install --prefix=$LOCAL_PYTHON pip virtualenv
