from fabric.api import run, env
from time import sleep
from boto.exception import EC2ResponseError
from .util import eval_template


def attach_volumes(vm_launcher, options, format=False):
    """
    """
    volumes = options.get("volumes", [])
    if not volumes:
        return
    boto_connection = vm_launcher.boto_connection()
    instance_id = run("curl --silent http://169.254.169.254/latest/meta-data/instance-id")
    for volume in volumes:
        volume_id = volume['id']
        device_id = volume['device']
        if not _get_attached(boto_connection, instance_id, device_id, valid_states=["attached", "attaching"]):
            boto_connection.attach_volume(volume_id, instance_id, device_id)
    for volume in volumes:
        volume_id = volume['id']
        device_id = volume['device']
        path = volume.get("path")

        while True:
            if _get_attached(boto_connection, instance_id, device_id):
                break

            sleep(5)
            print "Waiting for volume corresponding to device %s to attach" % device_id
            break

        # Don't mount if already mounted
        if _find_mounted_device_id(path):
            continue

        format = str(volume.get('format', "False")).lower()
        if format == "true":
            _format_device(device_id)
        env.safe_sudo("mkdir -p '%s'" % path)
        try:
            _mount(device_id, path)
        except:
            if format == "__auto__":
                print "Failed to mount device. format is set to __auto__ so will now format device and retry mount"
                _format_device(device_id)
                _mount(device_id, path)
            else:
                raise


def _mount(device_id, path):
    env.safe_sudo("mount '%s' '%s'" % (device_id, path))


def _format_device(device_id):
    env.safe_sudo("mkfs -t ext3 %s" % device_id)


def detach_volumes(vm_launcher, options):
    volumes = options.get("volumes", [])
    if not volumes:
        return

    boto_connection = vm_launcher.boto_connection()
    instance_id = run("curl --silent http://169.254.169.254/latest/meta-data/instance-id")
    for volume in volumes:
        volume_id = volume['id']
        path = volume.get("path")
        env.safe_sudo("umount '%s'" % path)
        _detach(boto_connection, instance_id, volume_id)


def make_snapshots(vm_launcher, options):
    volumes = options.get("volumes", [])
    for volume in volumes:
        path = volume.get("path")
        desc = volume.get("description", "Snapshot of path %s" % path)
        desc = eval_template(env, desc)
        # Allow volume to specify it should not be snapshotted, e.g. if
        # piggy backing on core teams snapshots for galaxyIndicies for instance.
        snapshot = volume.get("snapshot", True)
        if snapshot:
            _make_snapshot(vm_launcher, path, desc)


def _get_attached(conn, instance_id, device_id, valid_states=['attached']):
    vol_list = conn.get_all_volumes()
    fs_vol = None
    for vol in vol_list:
        if vol.attach_data.instance_id == instance_id and vol.attach_data.device == device_id:
            if vol.attach_data.status in valid_states:
                fs_vol = vol
                break
    return fs_vol


def _make_snapshot(vm_launcher, fs_path, desc):
    """ Create a snapshot of an existing volume that is currently attached to an
    instance, taking care of the unmounting and detaching. If you specify the
    optional argument (:galaxy), the script will pull the latest Galaxy code
    from bitbucket and perform an update before snapshotting. Else, the script
    will prompt for the file system path to be snapshoted.

    In order for this to work, an instance on EC2 needs to be running with a
    volume that wants to be snapshoted attached and mounted. The script will
    unmount the volume, create a snaphost and offer to reattach and mount the
    volume or create a new one from the freshly created snapshot.

    Except for potentially Galaxy, MAKE SURE there are no running processes
    using the volume and that no one is logged into the instance and sitting
    in the given directory.
    """
    instance_id = run("curl --silent http://169.254.169.254/latest/meta-data/instance-id")
    availability_zone = run("curl --silent http://169.254.169.254/latest/meta-data/placement/availability-zone")
    instance_region = availability_zone[:-1]  # Truncate zone letter to get region name
    # Find the device where the file system is mounted to
    # Find the EBS volume where the file system resides
    device_id = _find_mounted_device_id(fs_path)
    ec2_conn = vm_launcher.boto_connection()
    fs_vol = _get_attached(ec2_conn, instance_id, device_id)
    if fs_vol:
        env.safe_sudo("umount %s" % fs_path)
        _detach(ec2_conn, instance_id, fs_vol.id)
        snap_id = _create_snapshot(ec2_conn, fs_vol.id, desc)
        # TODO: Auto Update snaps?
        make_public = True
        if make_public:  # Make option
            ec2_conn.modify_snapshot_attribute(snap_id, attribute='createVolumePermission', operation='add', groups=['all'])
        reattach = True
        if reattach:
            _attach(ec2_conn, instance_id, fs_vol.id, device_id)
            env.safe_sudo("mount %s %s" % (device_id, fs_path))
        delete_old_volume = False
        if delete_old_volume:
            _delete_volume(ec2_conn, fs_vol.id)
        print "----- Done snapshoting volume '%s' for file system '%s' -----" % (fs_vol.id, fs_path)
    else:
        print "ERROR: Failed to find require file system, is boto installed? Is it not actually mounted?"


def _find_mounted_device_id(path):
    # Adding dollar sign to grep to distinguish between /mnt/galaxy and /mnt/galaxyIndices
    device_id = env.safe_sudo("df | grep '%s$' | awk '{print $1}'" % path)
    return device_id


def _attach(ec2_conn, instance_id, volume_id, device):
    """
    Attach EBS volume to the given device (using boto).
    Try it for some time.
    """
    try:
        print "Attaching volume '%s' to instance '%s' as device '%s'" % (volume_id, instance_id, device)
        volumestatus = ec2_conn.attach_volume(volume_id, instance_id, device)
    except EC2ResponseError, e:
        print "Attaching volume '%s' to instance '%s' as device '%s' failed. Exception: %s" % (volume_id, instance_id, device, e)
        return False

    for counter in range(30):
        print "Attach attempt %s, volume status: %s" % (counter, volumestatus)
        if volumestatus == 'attached':
            print "Volume '%s' attached to instance '%s' as device '%s'" % (volume_id, instance_id, device)
            break
        if counter == 29:
            print "Volume '%s' FAILED to attach to instance '%s' as device '%s'. Aborting." % (volume_id, instance_id, device)
            return False
        volumes = ec2_conn.get_all_volumes([volume_id])
        volumestatus = volumes[0].attachment_state()
        sleep(3)
    return True


def _detach(ec2_conn, instance_id, volume_id):
    """
    Detach EBS volume from the given instance (using boto).
    Try it for some time.
    """
    try:
        volumestatus = ec2_conn.detach_volume( volume_id, instance_id, force=True )
    except EC2ResponseError, ( e ):
        print "Detaching volume '%s' from instance '%s' failed. Exception: %s" % ( volume_id, instance_id, e )
        return False

    for counter in range( 30 ):
        print "Volume '%s' status '%s'" % ( volume_id, volumestatus )
        if volumestatus == 'available':
            print "Volume '%s' successfully detached from instance '%s'." % ( volume_id, instance_id )
            break
        if counter == 29:
            print "Volume '%s' FAILED to detach to instance '%s'." % ( volume_id, instance_id )
        sleep(3)
        volumes = ec2_conn.get_all_volumes( [volume_id] )
        volumestatus = volumes[0].status


def _delete_volume(ec2_conn, vol_id):
    try:
        ec2_conn.delete_volume(vol_id)
        print "Deleted volume '%s'" % vol_id
    except EC2ResponseError, e:
        print "ERROR deleting volume '%s': %s" % (vol_id, e)


def _create_snapshot(ec2_conn, volume_id, description=None):
    """
    Create a snapshot of the EBS volume with the provided volume_id.
    Wait until the snapshot process is complete (note that this may take quite a while)
    """
    snapshot = ec2_conn.create_snapshot(volume_id, description=description)
    if snapshot:
        while snapshot.status != 'completed':
            sleep(6)
            snapshot.update()
        print "Creation of snapshot for volume '%s' completed: '%s'" % (volume_id, snapshot)
        return snapshot.id
    else:
        print "Could not create snapshot from volume with ID '%s'" % volume_id
        return False
