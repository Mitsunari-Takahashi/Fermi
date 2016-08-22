#!/usr/bin/env python

import sys
import command
import subprocess

def ls_list(path_files):
    """Returns the results of 'ls <path_files>' command in list format.
    """
    ret_ls = commands.getoutput("{0}".format(path_files))
    li_path_files = ret_ls.split("\n")
    return li_path_files
