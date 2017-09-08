#!/usr/bin/env python

import sys
import os
import pickle

keys = sys.argv[2:]


def load(path_file):
    with open(path_file, mode='rb') as f:
        a = pickle.load(f)
    return a


def dump(path_file, obj):
    """Serialize your object to path_file.
"""
    with open(path_file, 'wb') as f:
        pickle.dump(obj, f)


def show(path_file, keys):
    """Print deserialized data.
"""
    nkey = len(keys)
    data = load(path_file)
    if nkey<1:
        for k, v in data.items():
            print '='*len(k)
            print k
            print '-'*len(k)
            print v
    else:
        parent = data
        for i in range(nkey):
            child = parent[keys[i]]
            parent = child
            #if i<nkey:
            print ' '*i, keys[i]
            #else:
        print ' '*(i+1), child


if __name__ == '__main__':
    show(sys.argv[1], keys)
