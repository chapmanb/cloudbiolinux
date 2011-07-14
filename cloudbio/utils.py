import logging

def _setup_logging(env):
    env.logger = logging.getLogger("cloudbiolinux")
    env.logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    env.logger.addHandler(ch)

