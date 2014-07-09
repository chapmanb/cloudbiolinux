#!/bin/sh

VIRTUALENV_VERSION=${VIRTUALENV_VERSION:-1.10.1}

cd `dirname $0`
PROJECT_DIRECTORY=${PROJECT_DIRECTORY:-`pwd`}
VENV_DIRECTORY=$PROJECT_DIRECTORY/.venv

if [ ! -e $VENV_DIRECTORY ];
then
    VIRTUALENV_URL="https://pypi.python.org/packages/source/v/virtualenv/virtualenv-${VIRTUALENV_VERSION}.tar.gz"
    DOWNLOAD_TAR_BALL="$PROJECT_DIRECTORY/virtualenv.tar.gz"
    PYTHON_VERSION=${PYTHON_VERSION:-`python -c "import sys; rev = sys.version_info; str = '%d.%d' % (rev[0], rev[1]); print str"`}
    LOCAL_PYTHON=$PROJECT_DIRECTORY/.virtualenv

    VIRTUALENV_PACKAGES_DIR="$LOCAL_PYTHON/lib/python${PYTHON_VERSION}/site-packages"
    VIRTUALENV_SOURCE_DIR=$PROJECT_DIRECTORY/.virtualenv_source
    
    wget -O "$DOWNLOAD_TAR_BALL" ${VIRTUALENV_URL}
    tar xzvf "$DOWNLOAD_TAR_BALL"
    mv virtualenv-${VIRTUALENV_VERSION} $VIRTUALENV_SOURCE_DIR
    mkdir -p $VIRTUALENV_PACKAGES_DIR

    export PYTHONPATH=$VIRTUALENV_PACKAGES_DIR:$PYTHONPATH
    export PATH=$LOCAL_PYTHON/bin:$PATH

    cd $VIRTUALENV_SOURCE_DIR
    python setup.py install --prefix="$LOCAL_PYTHON"
    cd $PROJECT_DIRECTORY
    virtualenv --no-site-packages $VENV_DIRECTORY
    . $VENV_DIRECTORY/bin/activate
    pip install -r $PROJECT_DIRECTORY/requirements.txt
fi

sh $PROJECT_DIRECTORY/deploy_no_deps.sh "$@"
