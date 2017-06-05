#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

def print_diff( frame_num ):
    '''Write some help documentation here
'''

    folder1 = (os.getcwd() + '/output/')
    folder2 = (os.getcwd() + '/output_2d/')

    filename = '/q00%02i.dat' % frame_num
    q1  = np.loadtxt(folder1 + filename )[1:]
    q2  = np.loadtxt(folder2 + filename )[1:]

    print('len q1 = %d, len q2 = %d' % (len(q1), len(q2) ) )
    diff = q1-q2

    diff = np.sqrt( np.dot(diff,diff) )

    print 'files %s differ by %2.15e' % (filename, diff )

def main():
    
    print_diff(0)
    print_diff(1)

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
