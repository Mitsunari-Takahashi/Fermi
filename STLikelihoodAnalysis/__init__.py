#!/usr/bin/env python
import logging


def get_module_logger(modname):
    logger = logging.getLogger(modname)
    sh = logging.StreamHandler()
    #formatter = logging.Formatter('my-format')
    #handler.setFormatter(formatter)
    logger.addHandler(sh)
    logger.setLevel(logging.INFO)
    return logger
