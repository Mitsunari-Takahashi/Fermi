#!/usr/bin/env python

import sys
import commands
import subprocess

def ls_list(path_files):
    """Returns the results of 'ls <path_files>' command in list format.
    """
    print path_files
    ret_ls = commands.getoutput("ls {0}".format(path_files))
    li_path_files = ret_ls.split("\n")
    print li_path_files
    return li_path_files
