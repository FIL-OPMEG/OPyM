#!/usr/bin/python

import inspect
import os.path as op

def root():
    import opym as test
    path = op.dirname(op.abspath(inspect.getfile(test)))
    return path