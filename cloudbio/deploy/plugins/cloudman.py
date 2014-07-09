from datetime import datetime
from os.path import exists, join
from os import listdir
from tempfile import mkdtemp

from cloudbio.deploy.util import eval_template

from boto.exception import S3ResponseError
from boto.s3.key import Key

import yaml

from fabric.api import local, lcd, env

DEFAULT_BUCKET_NAME = 'cloudman'

DEFAULT_CLOUDMAN_PASSWORD = 'adminpass'
DEFAULT_CLOUDMAN_CLUSTER_NAME = 'cloudman'


def bundle_cloudman(vm_launcher, options):
    cloudman_options = options.get('cloudman')
    cloudman_repository_path = cloudman_options['cloudman_repository']
    delete_repository = False
    bucket_source = cloudman_options.get("bucket_source")
    if cloudman_repository_path.startswith("http"):
        # Not a local path, lets clone it out of a remote repostiroy,
        temp_directory = mkdtemp()
        if cloudman_repository_path.endswith(".git"):
            branch_opts = ""
            repository_branch = cloudman_options.get('repository_branch', None)
            if repository_branch:
                branch_opts = "-b '%s'" % repository_branch
            clone_command = "git clone " + branch_opts + " '%s' '%s'"
        else:
            clone_command = "hg clone '%s' '%s'"
        local(clone_command % (cloudman_repository_path, temp_directory))
        cloudman_repository_path = temp_directory
        delete_repository = True
    try:
        with lcd(cloudman_repository_path):
            try:
                local("tar czvf cm.tar.gz *")
                local("mv cm.tar.gz '%s'" % bucket_source)
            finally:
                local("rm -f cm.tar.gz")
    finally:
        if delete_repository:
            local("rm -rf '%s'" % cloudman_repository_path)


def cloudman_launch(vm_launcher, options):
    cloudman_options = options.get('cloudman')
    image_id = cloudman_options.get('image_id', None)
    if str(image_id).lower() == "__use_snaps__":
        # TODO: Make more flexible
        bucket_source = cloudman_options.get("bucket_source")
        snaps_path = join(bucket_source, "snaps.yaml")
        if not exists(snaps_path):
            raise Exception("CloudMan AMI set to __use_snaps__ but now snaps.yaml file could be found with path %s" % snaps_path)
        snaps = {}
        with open(snaps_path, "r") as in_handle:
            snaps = yaml.load(in_handle)
        clouds = snaps["clouds"]
        if len(clouds) != 1:
            raise Exception("Exactly one cloud must be defined snaps.yaml for the deployer's CloudMan launch to work.")
        regions = clouds[0]["regions"]
        if len(regions) != 1:
            raise Exception("Exactly one region must be defined snaps.yaml for the deployer's CloudMan launch to work.")
        deployments = regions[0]["deployments"]
        if len(deployments) != 1:
            raise Exception("Exactly one deployment must be defined snaps.yaml for the deployer's CloudMan launch to work.")
        image_id = deployments[0]["default_mi"]

    size_id = cloudman_options.get('size_id', None)
    user_data = _prepare_user_data(vm_launcher, cloudman_options)
    vm_launcher.create_node('cloudman',
                            image_id=image_id,
                            size_id=size_id,
                            ex_userdata=user_data)


def sync_cloudman_bucket(vm_launcher, options):
    bucket = options.get("target_bucket", None)
    if not bucket:
        bucket = __get_bucket_default(options)
    bucket_source = options.get("cloudman", {}).get("bucket_source", None)
    if not bucket or not bucket_source:
        print "Warning: Failed to sync cloud bucket, bucket or bucket_source is undefined."
        return
    conn = vm_launcher.boto_s3_connection()
    for file_name in listdir(bucket_source):
        _save_file_to_bucket(conn, bucket, file_name, join(bucket_source, file_name))


def _save_file_to_bucket(conn, bucket_name, remote_filename, local_file, **kwargs):
    """ Save the local_file to bucket_name as remote_filename. Also, any additional
    arguments passed as key-value pairs, are stored as file's metadata on S3."""
    # print "Establishing handle with bucket '%s'..." % bucket_name
    b = _get_bucket(conn, bucket_name)
    if b is not None:
        # print "Establishing handle with key object '%s'..." % remote_filename
        k = Key( b, remote_filename )
        print "Attempting to save file '%s' to bucket '%s'..." % (remote_filename, bucket_name)
        try:
            # Store some metadata (key-value pairs) about the contents of the file being uploaded
            # Note that the metadata must be set *before* writing the file
            k.set_metadata('date_uploaded', str(datetime.utcnow()))
            for args_key in kwargs:
                print "Adding metadata to file '%s': %s=%s" % (remote_filename, args_key, kwargs[args_key])
                k.set_metadata(args_key, kwargs[args_key])
            print "Saving file '%s'" % local_file
            k.set_contents_from_filename(local_file)
            print "Successfully added file '%s' to bucket '%s'." % (remote_filename, bucket_name)
            make_public = True
            if make_public:
                k.make_public()
        except S3ResponseError, e:
            print "Failed to save file local file '%s' to bucket '%s' as file '%s': %s" % ( local_file, bucket_name, remote_filename, e )
            return False
        return True
    else:
        return False


def __get_bucket_default(options):
    cloudman_options = options.get("cloudman", {})
    user_data = cloudman_options = cloudman_options.get('user_data', None) or {}
    bucket = user_data.get("bucket_default", None)
    return bucket


def _prepare_user_data(vm_launcher, cloudman_options):
    cloudman_user_data = cloudman_options.get('user_data', None) or {}
    cluster_name = \
        cloudman_options.get('cluster_name', DEFAULT_CLOUDMAN_CLUSTER_NAME)
    password = cloudman_options.get('password', DEFAULT_CLOUDMAN_PASSWORD)
    access_key = vm_launcher.access_id()
    secret_key = vm_launcher.secret_key()

    _set_property_if_needed(cloudman_user_data, 'access_key', access_key)
    _set_property_if_needed(cloudman_user_data, 'secret_key', secret_key)
    cluster_name = eval_template(env, cluster_name)
    _set_property_if_needed(cloudman_user_data, 'cluster_name', cluster_name)
    _set_property_if_needed(cloudman_user_data, 'password', password)

    return yaml.dump(cloudman_user_data)


def _set_property_if_needed(user_data, property, value):
    if property not in user_data:
        user_data[property] = value


def _get_bucket(s3_conn, bucket_name):
    b = None
    for i in range(0, 5):
        try:
            b = s3_conn.get_bucket(bucket_name)
            break
        except S3ResponseError:
            print "Bucket '%s' not found, attempt %s/5" % (bucket_name, i)
            return None
    return b


local_actions = {
    "cloudman_launch": cloudman_launch,
    "sync_cloudman_bucket": sync_cloudman_bucket,
    "bundle_cloudman": bundle_cloudman,
}
