#!/usr/bin/env python
"""
Requires: 
    PyYAML http://pyyaml.org/wiki/PyYAMLDocumentation (easy_install pyyaml)
    boto http://code.google.com/p/boto/ (easy_install boto)

Assumptions:
    DEFAULT_BUCKET_NAME and DEFAULT_BOOT_SCRIPT_NAME are publicly accessible and
    do not require any form of authentication
"""

import os, sys, yaml, urllib2, logging, hashlib, time, subprocess, random
from urlparse import urlparse

from boto.s3.connection import S3Connection
from boto.s3.key import Key
from boto.exception import S3ResponseError
logging.getLogger('boto').setLevel(logging.INFO) # Only log boto messages >=INFO


USER_DATA_URL = 'http://169.254.169.254/latest/user-data'
# USER_DATA_URL = 'http://userwww.service.emory.edu/~eafgan/content/userData.yaml.sample' # used for testing
# USER_DATA_URL = 'http://userwww.service.emory.edu/~eafgan/content/url_ud.txt' # used for testing
LOCAL_PATH = '/tmp/cm' # Local path destination used for storing/reading any files created by this script
USER_DATA_FILE_NAME = 'userData.yaml' # Local file with user data formatted by this script
USER_DATA_FILE = os.path.join(LOCAL_PATH, USER_DATA_FILE_NAME) 
USER_DATA_ORIG = os.path.join(LOCAL_PATH, 'original_%s' % USER_DATA_FILE_NAME) # Local file containing user data in its original format
SERVICE_ROOT = 'http://s3.amazonaws.com/' # Obviously, customized for Amazon's S3
DEFAULT_BUCKET_NAME = 'cloudman' # Ensure this bucket is accessible to anyone!
DEFAULT_BOOT_SCRIPT_NAME = 'cm_boot.py' # Ensure this file is accessible to anyone in the public bucket!
CLOUDMAN_HOME = '/mnt/cm'

# ====================== Utility methods ======================

def _get_user_data():
    for i in range(0, 5):
          try:
              log.info("Getting user data from '%s', attempt %s" % (USER_DATA_URL, i))
              fp = urllib2.urlopen(USER_DATA_URL)
              ud = fp.read()
              fp.close()
              log.debug("Saving user data in its original format to file '%s'" % USER_DATA_ORIG)
              with open(USER_DATA_ORIG, 'w') as ud_orig:
                  ud_orig.write(ud)
              if ud:
                  log.debug("Got user data")
                  return ud
          except IOError:
              log.info("User data not found. Setting it to empty.")
              return ''
    # Used for testing
    # return 'http://s3.amazonaws.com/cloudman/cm_boot'
    # return ''
    # return "gc_dev1|<account_key>|<secret_key>|somePWD"
    # with open('sample.yaml') as ud_yaml:
    #     ud = ud_yaml.read()
    # return ud

def _get_bucket_name(cluster_name, access_key):
    """Compose bucket name based on the user-provided cluster name and user access key""" 
    m = hashlib.md5()
    m.update( cluster_name + access_key )
    return "cm-" + m.hexdigest()

def _isurl(path):
    """Test if path is a net location. Tests the scheme and netloc."""
    # BUG : URLs require a scheme string ('http://') to be used.
    #       www.google.com will fail.
    #       Should we prepend the scheme for those that don't have it and
    #       test that also?
    scheme, netloc, upath, uparams, uquery, ufrag = urlparse(path)
    return bool(scheme and netloc)

def _get_s3_conn(access_key, secret_key):
    log.debug('Establishing boto S3 connection')
    try:
        s3_conn = S3Connection(access_key, secret_key)
        log.debug('Got boto S3 connection.')
    except Exception, e:
        log.error("Exception getting S3 connection: %s" % e)
    return s3_conn

def _bucket_exists(s3_conn, bucket_name):
    bucket = None
    for i in range(1, 6):
		try:
            # log.debug("Looking for bucket '%s'" % bucket_name)
		    bucket = s3_conn.lookup(bucket_name)
		    break
		except S3ResponseError: 
		    log.error ("Bucket '%s' not found, attempt %s/5" % (bucket_name, i+1))
		    time.sleep(2)
		    
    if bucket is not None:
        log.debug("Cluster bucket '%s' found." % bucket_name)
        return True
    else:
        log.debug("Cluster bucket '%s' not found." % bucket_name)
        return False

def _remote_file_exists(s3_conn, bucket_name, remote_filename):
    b = None
    for i in range(0, 5):
		try:
			b = s3_conn.get_bucket(bucket_name)
			break
		except S3ResponseError: 
			log.error ("Problem connecting to bucket '%s', attempt %s/5" % (bucket_name, i))
			time.sleep(2)
	    	
    if b is not None:
		k = Key(b, remote_filename)
		if k.exists():
		    return True
    return False

def _save_file_to_bucket(s3_conn, bucket_name, remote_filename, local_file, force=False):
    local_file = os.path.join(LOCAL_PATH, local_file)
    # log.debug( "Establishing handle with bucket '%s'..." % bucket_name)
    b = None
    for i in range(0, 5):
		try:
			b = s3_conn.get_bucket(bucket_name)
			break
		except S3ResponseError, e:
			log.error ("Problem connecting to bucket '%s', attempt %s/5" % (bucket_name, i))
			time.sleep(2)
    
    if b is not None:
        # log.debug("Establishing handle with key object '%s'..." % remote_filename)
        k = Key(b, remote_filename)
        if k.exists() and not force:
            log.debug("Remote file '%s' already exists. Not overwriting it." % remote_filename)
            return True
        log.debug( "Attempting to save local file '%s' to bucket '%s' as '%s'" % (local_file, bucket_name, remote_filename))
        try:
            k.set_contents_from_filename(local_file)
            log.info( "Successfully saved file '%s' to bucket '%s'." % (remote_filename, bucket_name))
            return True
        except S3ResponseError, e:
             log.error("Failed to save file local file '%s' to bucket '%s' as file '%s': %s" % (local_file, bucket_name, remote_filename, e))
             return False
    else:
        return False

def _get_file_from_bucket(s3_conn, bucket_name, remote_filename, local_filename):
    local_filename = os.path.join(LOCAL_PATH, local_filename)
    try:
        # log.debug("Establishing handle with bucket '%s'" % bucket_name)
        b = s3_conn.get_bucket(bucket_name)
        
        # log.debug("Establishing handle with file object '%s'" % remote_filename)
        k = Key(b, remote_filename)
    
        log.debug("Attempting to retrieve file '%s' from bucket '%s'" % (remote_filename, bucket_name))
        if k.exists():
            k.get_contents_to_filename(local_filename)
            log.info("Successfully retrieved file '%s' from bucket '%s' to '%s'." % (remote_filename, bucket_name, local_filename))
            return True
        else:
            log.error("File '%s' in bucket '%s' not found." % (remote_filename, bucket_name))
            return False
    except S3ResponseError, e:
        log.error("Failed to get file '%s' from bucket '%s': %s" % (remote_filename, bucket_name, e))
        return False

def _get_file_from_url(url):
    local_filename = os.path.join(LOCAL_PATH, os.path.split(url)[1])
    log.info("Getting boot script from '%s' and saving it locally to '%s'" % (url, local_filename))
    try:
        f = urllib2.urlopen(url)
        with open(local_filename, 'w') as local_file:
            local_file.write(f.read())
        os.chmod(local_filename, 0744) 
        if f:
            log.debug("Got boot script from '%s'" % url)
            return True
        return False
    except IOError:
        log.error("Boot script at '%s' not found." % url)
        return False

def _get_boot_script(ud):
    # Test if user's bucket exists; if it does not, resort to the default
    # bucket for downloading the boot script
    use_default_bucket = False
    if ud.has_key('bucket_default'):
        default_bucket_name = ud['bucket_default']
    else:
        default_bucket_name = DEFAULT_BUCKET_NAME
    if ud.has_key('bucket_cluster') and ud['access_key'] is not None and ud['secret_key'] is not None:
        s3_conn = _get_s3_conn(ud['access_key'], ud['secret_key'])
        # Check if user's bucket exists or use the default one
        if not _bucket_exists(s3_conn, ud['bucket_cluster']) or not _remote_file_exists(s3_conn, ud['bucket_cluster'], ud['boot_script_name']):
            log.debug("Using default bucket '%s'" % default_bucket_name)
            use_default_bucket = True
        else:
            log.debug("Using cluster bucket '%s'" % ud['bucket_cluster'])
            use_default_bucket = False
    else:
        log.debug("Defaulting to bucket '%s'" % default_bucket_name)
        use_default_bucket = True
        
    # If using user's bucket, use credentials because the script may not be accessible to everyone
    got_boot_script = False
    if not use_default_bucket:
        log.debug("Trying to get boot script '%s' from bucket '%s'" % (ud['boot_script_name'], ud['bucket_cluster']))
        got_boot_script = _get_file_from_bucket(s3_conn, ud['bucket_cluster'], ud['boot_script_name'], DEFAULT_BOOT_SCRIPT_NAME)
        if got_boot_script:
            os.chmod(os.path.join(LOCAL_PATH, DEFAULT_BOOT_SCRIPT_NAME), 0744)
    # If did not get boot script, fall back on the publicly available one
    if not got_boot_script or use_default_bucket:
        boot_script_url = os.path.join(_get_default_bucket_url(ud), ud['boot_script_name'])
        log.debug("Could not get boot script '%s' from cluster bucket '%s'; retrieving public one from bucket url '%s'" % (ud['boot_script_name'], ud['bucket_cluster'], boot_script_url))
        got_boot_script = _get_file_from_url(boot_script_url)
    if got_boot_script:
        log.debug("Saved boot script to '%s'" % os.path.join(LOCAL_PATH, DEFAULT_BOOT_SCRIPT_NAME))
        # Save downloaded boot script to user's bucket for future invocations 
        if ud.has_key('bucket_cluster') and ud['bucket_cluster']:
            s3_conn = _get_s3_conn(ud['access_key'], ud['secret_key'])
            if _bucket_exists(s3_conn, ud['bucket_cluster']) and not _remote_file_exists(s3_conn, ud['bucket_cluster'], ud['boot_script_name']):
                _save_file_to_bucket(s3_conn, ud['bucket_cluster'], ud['boot_script_name'], DEFAULT_BOOT_SCRIPT_NAME)        
        return True
    log.debug("**Could not get boot script**")
    return False

def _run_boot_script(boot_script_name):
    script = os.path.join(LOCAL_PATH, boot_script_name)
    log.info("Running boot script '%s'" % script)
    process = subprocess.Popen(script, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0:
        log.debug("Successfully ran boot script '%s'" % script)
        return True
    else:
        log.error("Error running boot script '%s'. Process returned code '%s' and following stderr: %s" % (script, process.returncode, stderr))
        return False

def _create_basic_user_data_file():
    # Create a basic YAML file that is expected by CloudMan 
    with open(USER_DATA_FILE, 'w') as ud_file:
        ud_formatted = {'access_key': None,
                        'boot_script_name': DEFAULT_BOOT_SCRIPT_NAME,
                        'boot_script_path': LOCAL_PATH,
                        'bucket_default': DEFAULT_BUCKET_NAME,
                        'bucket_cluster': None,
                        'cloudman_home': CLOUDMAN_HOME,
                        'cluster_name': 'aGalaxyCloudManCluster_%s' % random.randrange(1, 9999999),
                        'role': 'master',
                        'secret_key': None,
                        'no_start': True}
        yaml.dump(ud_formatted, ud_file, default_flow_style=False)
    return ud_formatted

def _get_default_bucket_url(ud=None):
    if ud and ud.has_key('bucket_default'):
        default_bucket_name = ud['bucket_default']
    else:
        default_bucket_name = DEFAULT_BUCKET_NAME
    # TODO: Check if th bucket 'default_bucket_name' is accessible to everyone 
    # because it is being accessed as a URL
    bucket_url = os.path.join(SERVICE_ROOT, default_bucket_name)
    print "Default bucket url: %s" % default_bucket_name
    return bucket_url

# ====================== Actions methods ======================

def _handle_empty():
    log.info("Received empty user data; assuming default contextualization")
    _create_basic_user_data_file() # This file is expected by CloudMan
    # Get & run boot script
    file_url = os.path.join(_get_default_bucket_url(), DEFAULT_BOOT_SCRIPT_NAME)
    log.debug("Resorting to the default bucket to get boot script: %s" % file_url)
    _get_file_from_url(file_url)
    _run_boot_script(DEFAULT_BOOT_SCRIPT_NAME)

def _handle_url(url):
    log.info("Handling user data provided URL: '%s'" % url)
    _get_file_from_url(url)
    boot_script_name = os.path.split(url)[1]
    _run_boot_script(boot_script_name)

def _handle_yaml(user_data):
    """ Process user data in YAML format"""
    log.info("Handling user data in yaml format")
    ud = yaml.load(user_data)

    # Handle bad user data as a string
    if ud == user_data:
        return _handle_empty()
    # Handle freenx passwords and the case with only a password sent
    if ud.has_key("freenxpass"):
        _handle_freenx(ud["freenxpass"])
        if len(ud) == 1:
            return _handle_empty()

    # Create a YAML file from user data and store it as USER_DATA_FILE
    # This code simply ensures fields required by CloudMan are in the 
    # created file. Any other fields that might be included as user data
    # are also included in the created USER_DATA_FILE
    if not ud.has_key('cluster_name'):
        log.warning("The provided user data should contain cluster_name field.")
        ud['cluster_name'] = 'aGalaxyCloudManCluster_%s' % random.randrange(1, 9999999)
    elif ud['cluster_name'] == '':
        log.warning("The cluster_name field of user data should not be empty.")
        ud['cluster_name'] = 'aGalaxyCloudManCluster_%s' % random.randrange(1, 9999999)
    
    if not ud.has_key('access_key'):
        log.info("The provided user data does not contain access_key field; setting it to None..")
        ud['access_key'] = None
    elif ud['access_key'] == '' or ud['access_key'] is None:
        log.warning("The access_key field of user data should not be empty; setting it to None.")
        ud['access_key'] = None
    
    if not ud.has_key('secret_key'):
        log.info("The provided user data does not contain secret_key field; setting it to None.")
        ud['secret_key'] = None
    elif ud['secret_key'] == '' or ud['secret_key'] is None:
        log.warning("The secret_key field of user data should not be empty; setting it to None.")
        ud['secret_key'] = None
    
    if not ud.has_key('password'):
        log.warning("The provided user data should contain password field.")
    elif ud['password'] == '':
        log.warning("The password field of user data should not be empty.")
    
    if not ud.has_key('bucket_default'):
        log.debug("The provided user data does not contain bucket_default field; setting it to '%s'." % DEFAULT_BUCKET_NAME)
        ud['bucket_default'] = DEFAULT_BUCKET_NAME
    elif ud['bucket_default'] == '':
        log.warning("The bucket_default field of user data was empty; setting it to '%s'." % DEFAULT_BUCKET_NAME)
        ud['bucket_default'] = DEFAULT_BUCKET_NAME
    
    if not ud.has_key('bucket_cluster'):
        if ud['access_key'] is not None and ud['secret_key'] is not None:
            ud['bucket_cluster'] = _get_bucket_name(ud['cluster_name'], ud['access_key'])
    
    if not ud.has_key('role'):
        ud['role'] = 'master'
    
    if not ud.has_key('cloudman_home'):
        ud['cloudman_home'] = CLOUDMAN_HOME
    
    if not ud.has_key('boot_script_name'):
        ud['boot_script_name'] = DEFAULT_BOOT_SCRIPT_NAME
    ud['boot_script_path'] = LOCAL_PATH # Marks where boot script was saved
    
    log.debug("Composed user data: %s" % ud)
    with open(USER_DATA_FILE, 'w') as ud_yaml:
        yaml.dump(ud, ud_yaml, default_flow_style=False)

    # Get & run boot script
    if _get_boot_script(ud):
        _run_boot_script(DEFAULT_BOOT_SCRIPT_NAME)

def _handle_freenx(passwd):
    user = "ubuntu"
    p1 = subprocess.Popen(["echo", "%s:%s" % (user, passwd)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["chpasswd"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    out = p2.communicate()[0]
    cl = ["sed", "-i", "s/^PasswordAuthentication .*/PasswordAuthentication yes/",
          "/etc/ssh/sshd_config"]
    subprocess.check_call(cl)
    cl = ["/etc/init.d/ssh", "reload"]
    subprocess.check_call(cl)
    cl = ["dpkg-reconfigure", "-pcritical", "freenx-server"]
    subprocess.check_call(cl)

def _parse_user_data(ud):
    if ud == '':
        _handle_empty()
    elif _isurl(ud):
        _handle_url(ud)
    elif ud.startswith("freenxpass="):
        _handle_freenx(ud.strip().replace("freenxpass=", ""))
    else: # default to yaml
        _handle_yaml(ud)

# ====================== Driver code ======================

if not os.path.exists(LOCAL_PATH):
    os.mkdir(LOCAL_PATH)

# Logging setup
formatter = logging.Formatter("[%(levelname)s] %(module)s:%(lineno)d %(asctime)s: %(message)s")
console = logging.StreamHandler() # log to console - used during testing
# console.setLevel(logging.INFO) # accepts >INFO levels
console.setFormatter(formatter)
log_file = logging.FileHandler(os.path.join(LOCAL_PATH, "%s.log" % sys.argv[0]), 'w') # log to a file
log_file.setLevel(logging.DEBUG) # accepts all levels
log_file.setFormatter(formatter)
log = logging.root
log.addHandler( console )
log.addHandler( log_file )
log.setLevel( logging.DEBUG )

ud = _get_user_data()
_parse_user_data(ud)
log.info("---> %s done <---" % sys.argv[0])
