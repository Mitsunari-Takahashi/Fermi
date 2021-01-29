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


def extract_from_container(container, lst_keys, bprint=False):
    extracted = container
    for i, key in enumerate(lst_keys):
        if (isinstance(extracted, list) or isinstance(extracted, tuple)) and isinstance(key, basestring):
            if key in [str(m) for m in range(-len(extracted)-1, len(extracted)+1)]:
                extracted = extracted[int(key)]
            else:
                print 'Wrong key {0} for list or tuple!!'.format(key)
        else:
            extracted = extracted[key]
        if bprint==True:
            print ' '*i, key
    if bprint==True:
        print ' '*(i+1), extracted
    return extracted
    

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
            key = keys[i]
            if key in [str(m) for m in range(-len(parent)-1, len(parent)+1)]:
                key = int(key)
            child = parent[key]
            parent = child
            #if i<nkey:
            print ' '*i, key
            #else:
        print ' '*(i+1), child


if __name__ == '__main__':
    show(sys.argv[1], keys)
