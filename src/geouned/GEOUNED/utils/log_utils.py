import logging

formatter = logging.Formatter("%(asctime)s :: %(levelname)s :: %(funcName)s :: %(lineno)d :: %(message)s")


def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger
