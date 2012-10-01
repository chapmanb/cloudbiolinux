import yaml

DEFAULT_CLOUDMAN_PASSWORD = 'adminpass'
DEFAULT_CLOUDMAN_CLUSTER_NAME = 'cloudman'


def cloudman_launch(vm_launcher, options):
    cloudman_options = options.get('cloudman')
    image_id = cloudman_options.get('image_id', None)
    size_id = cloudman_options.get('size_id', None)
    user_data = _prepare_user_data(vm_launcher, cloudman_options)
    vm_launcher.create_node('cloudman',
                            image_id=image_id,
                            size_id=size_id,
                            ex_userdata=user_data)


def _prepare_user_data(vm_launcher, cloudman_options):
    cloudman_user_data = cloudman_options.get('user_data', {})
    cluster_name = \
        cloudman_options.get('cluster_name', DEFAULT_CLOUDMAN_CLUSTER_NAME)
    password = cloudman_options.get('password', DEFAULT_CLOUDMAN_PASSWORD)
    access_key = vm_launcher.access_id()
    secret_key = vm_launcher.secret_key()

    _set_property_if_needed(cloudman_user_data, 'access_key', access_key)
    _set_property_if_needed(cloudman_user_data, 'secret_key', secret_key)
    _set_property_if_needed(cloudman_user_data, 'cluster_name', cluster_name)
    _set_property_if_needed(cloudman_user_data, 'password', password)

    return yaml.dump(cloudman_user_data)


def _set_property_if_needed(user_data, property, value):
    if property not in user_data:
        user_data[property] = value
