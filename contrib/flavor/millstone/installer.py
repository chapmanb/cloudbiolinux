from fabric.api import *
from fabric.contrib.files import *
import os.path
import os

INSTALLATION_PATH = "$HOME/millstone"

REPO_URL = "git@github.com:churchlab/millstone.git"
REPO_STABLE_COMMIT = "634db63de2fb275719868839bf44bd9b4b9f016e"

BOOTSTRAP_SCRIPT = """#!/bin/bash
set -x
export RUN_MASTER=%d
export RUN_WORKER=%d
export PROJECT_DIR="%s"
export MILLSTONE_PATH="%s"
export CELERY_MANAGER_PATH="%s"
export JBROWSE_PATH="%s"
export TEMP_FILE_PATH="%s"
export PASSWORD_FILE_PATH="%s"
export BOOTSTRAP_FINISH_PATH="%s"
export RUN_AS_USER="%s"
export EC2_HOSTNAME=$(curl http://169.254.169.254/latest/meta-data/public-hostname)

export MILLSTONE_WEB_PORT=8000
export CELERY_MANAGER_WEB_PORT=8001

genpw() {
    python -c "import string,random;print ''.join(random.choice(string.letters + string.digits) for x in range(10))"
}

cpucount() {
    grep -c ^processor /proc/cpuinfo
}

export NUM_CPU=$(cpucount)

export TIMEOUT=3600

if [ "$RUN_MASTER" == "1" ];
then
    echo "Configuring master..."

    export PASSWORD=$(genpw)
    export RABBITMQ_USER="millstone"
    export POSTGRES_DB="millstone"
    export POSTGRES_USER="millstone"

    echo ${PASSWORD} > ${PASSWORD_FILE_PATH}
    echo "Using '${PASSWORD}' as password for PostgresSQL and RabbitMQ"

    /etc/init.d/rabbitmq-server stop
    cat > /etc/rabbitmq/rabbitmq-env.conf << EOF
NODE_IP_ADDRESS=0.0.0.0
NODE_PORT=5672
NODENAME="rabbit@localhost"
EOF

    /etc/init.d/rabbitmq-server start

    # Setup user in RabbitMQ and make it public.
    rabbitmqctl change_password guest $(genpw)
    RABBITMQ_USERS=$(sudo rabbitmqctl list_users -q)
    if [[ "${RABBITMQ_USERS}" =~ "${RABBITMQ_USER}" ]]
    then
        echo "Deleting existing ${RABBITMQ_USER} user in RabbitMQ."
        rabbitmqctl delete_user ${RABBITMQ_USER}
    fi
    rabbitmqctl add_user ${RABBITMQ_USER} ${PASSWORD}
    rabbitmqctl set_permissions -p / ${RABBITMQ_USER} ".*" ".*" ".*"

    /etc/init.d/rabbitmq-server restart

    # Setup user and database in PostgresSQL
    sudo -u postgres psql -U postgres -d postgres -c "DROP DATABASE IF EXISTS ${POSTGRES_DB};"
    sudo -u postgres psql -U postgres -d postgres -c "DROP USER IF EXISTS ${POSTGRES_USER};"
    sudo -u postgres psql -U postgres -d postgres -c "CREATE USER ${POSTGRES_USER} WITH PASSWORD '${PASSWORD}';"
    sudo -u postgres psql -U postgres -d postgres -c "CREATE DATABASE ${POSTGRES_DB};"
    sudo -u postgres psql -U postgres -d postgres -c "GRANT ALL PRIVILEGES ON DATABASE ${POSTGRES_DB} to ${POSTGRES_USER};"
    sudo -u postgres psql -U postgres -d postgres -c "ALTER USER ${POSTGRES_USER} CREATEDB;"

    POSTGRES_CONF=$(find /etc/postgresql -name "postgresql.conf" | head -1)
    echo "listen_addresses='*'" >> ${POSTGRES_CONF}

    PG_HDA_CONF=$(find /etc/postgresql -name "pg_hba.conf" | head -1)
    echo "host all all  0.0.0.0/0 md5" >> ${PG_HDA_CONF}

    /etc/init.d/postgresql restart

    LOCAL_SETTINGS_PATH="${MILLSTONE_PATH}/conf/local_settings.py"

    if [[ ! -e ${LOCAL_SETTINGS_PATH} ]]
    then
        echo "${LOCAL_SETTINGS_PATH} does not exist!"
    fi

    cat > ${LOCAL_SETTINGS_PATH} << EOF
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': '${POSTGRES_DB}',
        'USER': '${POSTGRES_USER}',
        'PASSWORD': '${PASSWORD}',
        'HOST': '127.0.0.1',
        'PORT': '5432',
        'OS_USER': 'postgres',
    }
}
BROKER_URL = 'amqp://${RABBITMQ_USER}:${PASSWORD}@127.0.0.1:5672//'

EOF
    chown ${RUN_AS_USER} ${LOCAL_SETTINGS_PATH}

    echo "New local settings in ${LOCAL_SETTINGS_PATH}:"
    cat ${LOCAL_SETTINGS_PATH}
fi

chown -R ${RUN_AS_USER} ${PROJECT_DIR}

echo "Reconfiguring nginx and supervisor..."
/etc/init.d/supervisor stop
/etc/init.d/nginx stop

rm -f /etc/nginx/sites-enabled/default
rm -f /etc/nginx/sites-available/millstone
rm -f /etc/nginx/sites-enabled/millstone
rm -f /etc/supervisor/supervisord.conf

cat > /etc/nginx/sites-available/millstone << EOF
server {
  server_name localhost;

  location /jbrowse {
    alias ${JBROWSE_PATH};
  }

  location /tmp {
    alias ${TEMP_FILE_PATH};
  }

  location /static {
    alias ${MILLSTONE_PATH}/main/static;
  }

  location / {
    proxy_pass http://127.0.0.1:${MILLSTONE_WEB_PORT}/;
  }

  location /worker/ {
    proxy_pass http://127.0.0.1:${CELERY_MANAGER_WEB_PORT}/;
  }

  location /worker/static {
    alias ${CELERY_MANAGER_PATH}/celery_manager/static;
  }

  # Override timeouts. Especially relevant to upload requests.
  proxy_connect_timeout ${TIMEOUT};
  proxy_read_timeout ${TIMEOUT};

  # No limit on upload size.
  client_max_body_size 0;
}
EOF

if [ "$RUN_WORKER" == "1" ];
then
    echo "Configuring worker..."
fi


ln -s /etc/nginx/sites-available/millstone /etc/nginx/sites-enabled/millstone
/etc/init.d/nginx start

cat > "/etc/supervisor/supervisord.conf" << EOF
[unix_http_server]
file=/var/run//supervisor.sock   ; (the path to the socket file)
chmod=0700                       ; sockef file mode (default 0700)

[inet_http_server]
port=127.0.0.1:9001

[supervisord]
logfile=/var/log/supervisor/supervisord.log ; (main log file;default $CWD/supervisord.log)
pidfile=/var/run/supervisord.pid ; (supervisord pidfile;default supervisord.pid)
childlogdir=/var/log/supervisor            ; ('AUTO' child log dir, default $TEMP)

[rpcinterface:supervisor]
supervisor.rpcinterface_factory = supervisor.rpcinterface:make_main_rpcinterface

[supervisorctl]
serverurl=http://127.0.0.1:9001

[program:millstone]
command=gunicorn_django -b 127.0.0.1:${MILLSTONE_WEB_PORT} --workers=$(expr ${NUM_CPU} + 1) --timeout=${TIMEOUT}
directory=${MILLSTONE_PATH}
autostart=${RUN_MASTER}
autorestart=True
redirect_stderr=True
user=${RUN_AS_USER}

[program:celery_manager]
command=gunicorn_django -b 127.0.0.1:${CELERY_MANAGER_WEB_PORT} --workers=$(expr ${NUM_CPU} + 1) --timeout=${TIMEOUT}
directory=${CELERY_MANAGER_PATH}
autostart=${RUN_WORKER}
autorestart=True
redirect_stderr=True
user=${RUN_AS_USER}

[program:celery]
command=python manage.py celery worker --loglevel=info
directory=${MILLSTONE_PATH}
autostart=${RUN_WORKER}
autorestart=True
redirect_stderr=True
user=${RUN_AS_USER}

EOF

/etc/init.d/supervisor start

update-rc.d nginx defaults
update-rc.d supervisor defaults

if [ "$RUN_MASTER" == "1" ];
then
    pushd ${MILLSTONE_PATH}
    sudo -u ${RUN_AS_USER} python manage.py syncdb --noinput
    sudo -u ${RUN_AS_USER} python manage.py migrate
    echo "y" | sudo -u ${RUN_AS_USER} python scripts/bootstrap_data.py
    popd
fi

touch ${BOOTSTRAP_FINISH_PATH}

echo "Bootstrap finished!"
"""

BOOTSTRAP_INVOKER_SCRIPT = """#!/bin/bash
### BEGIN INIT INFO
# Provides:          millstone_setup
# Required-Start:    $all
# Required-Stop:     $remote_fs $syslog
# Default-Start:     2 3 4 5
# Default-Stop:      0 1 6
### END INIT INFO

export PROJECT_DIR="%s"
export BOOTSTRAP_FINISH_PATH="%s"
export BOOTSTRAP_SCRIPT_PATH="%s"

# Carry out specific functions when asked to by the system
case "$1" in
  start)
    if [[ ! -e "${BOOTSTRAP_FINISH_PATH}" ]];
      then
        /bin/bash ${BOOTSTRAP_SCRIPT_PATH} 2>&1 > "${PROJECT_DIR}/bootstrap.log"
      fi
    ;;
  stop)
    echo "/etc/init.d/millstone_setup stop"
    ;;
  *)
    echo "Usage: /etc/init.d/millstone_setup {start|stop}"
    exit 1
    ;;
esac

exit 0

"""


def install_millstone(env):
    current_user = env.safe_run("echo $USER").strip()
    home_dir = env.safe_run("echo $HOME").strip()
    installation_dir = env.safe_run("echo %s" % INSTALLATION_PATH).strip()

    VM_MODE = {
        'MASTER': os.getenv("MASTER") is not None,
        'WORKER': os.getenv("WORKER") is not None,
    }
    env.logger.info("VM_MODE: %r" % VM_MODE)

    if env.safe_exists(installation_dir):
        env.logger.warning("%s already exists. Removing the directory. " % installation_dir)
        with cd(installation_dir):
            env.safe_sudo("rm -rf %s" % installation_dir)

    env.safe_run("mkdir -p %s" % installation_dir)
    env.logger.info("Installing Millstone to %s" % installation_dir)

    env.logger.debug("Configure SSH to ignore host checking for Github...")
    env.safe_run("mkdir -p %s" % os.path.join(home_dir, ".ssh"))
    env.safe_run("chmod 700 ~/.ssh")
    env.safe_run("chmod 600 ~/.ssh/*")
    append("~/.ssh/config", "Host github.com\n\tStrictHostKeyChecking no\n")

    with cd(installation_dir):
        env.safe_run("git clone %s %s" % (REPO_URL, installation_dir))
        env.safe_run("git reset --hard %s" % REPO_STABLE_COMMIT)

        project_dir = installation_dir
        jbrowse_dir = os.path.join(project_dir, "jbrowse")
        genome_designer_dir = os.path.join(project_dir, "genome_designer")
        temp_file_dir = os.path.join(genome_designer_dir, "temp_data/tmp")
        celery_manager_dir = os.path.join(project_dir, "celery_manager")
        config_dir = os.path.join(project_dir, "config")
        assert env.safe_exists(project_dir)
        assert env.safe_exists(jbrowse_dir)
        assert env.safe_exists(genome_designer_dir)
        assert env.safe_exists(celery_manager_dir)
        assert env.safe_exists(config_dir)

    with cd(project_dir):
        # clone JBrowse
        env.safe_run("git submodule update --init --recursive")

    with cd(jbrowse_dir):
        # Setup JBrowse.
        env.safe_run("./setup.sh")

    # Install millstone python requirements.
    python_requirements_file = os.path.join(installation_dir,
            'requirements', 'deploy.txt')
    env.safe_sudo("pip install -r %s" % python_requirements_file)

    with cd(genome_designer_dir):
        env.safe_run("ln -s ../jbrowse jbrowse")
        env.safe_run("./millstone_setup.py")

    run_worker = 1 if VM_MODE['WORKER'] else 0
    run_master = 1 if VM_MODE['MASTER'] else 0

    bootstrap_script_path = os.path.join(project_dir, "bootstrap.sh")
    bootstrap_finish_path = os.path.join(project_dir, "BOOTSTRAPPED")
    bootstrap_script = BOOTSTRAP_SCRIPT % (run_master, run_worker,
        project_dir, genome_designer_dir, celery_manager_dir, jbrowse_dir,
        temp_file_dir, os.path.join(project_dir, "password.txt"),
        bootstrap_finish_path, current_user)

    append(bootstrap_script_path, bootstrap_script)
    env.safe_run("chmod +x %s" % bootstrap_script_path)

    """
    /etc/init.d/millstone_setup will check if bootstrap script has
    run before, and execute bootstrap if not.
    """
    append("/etc/init.d/millstone_setup", BOOTSTRAP_INVOKER_SCRIPT % (project_dir,
        bootstrap_finish_path, bootstrap_script_path), use_sudo=True)
    env.safe_sudo("sudo chmod +x /etc/init.d/millstone_setup")
    env.safe_sudo("sudo update-rc.d millstone_setup defaults")
