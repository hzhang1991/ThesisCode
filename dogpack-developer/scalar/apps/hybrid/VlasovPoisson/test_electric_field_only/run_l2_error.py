#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

def main( ):
    '''Simple python script for computing the l2-error of the coefficients
    printed to output.
'''

    print "# leading comments can be given a '#' character"
    my_dictionary = {}
    old_err = i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output%(dir_num)03i/') % my_dictionary

        if( not os.path.exists(folder) ):
            print 'did not find folder: %s' % folder
            break

        my_dictionary['curr_folder'] = folder
        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = i

        print(folder)
        try:
            qex  = np.loadtxt(folder + "/c0000.dat")[1:]
            qapp = np.loadtxt(folder + "/c0001.dat")[1:]
        except IOError:
            print('Trial # %d has not finished running' % i )
            break

        # difference and new error:
        diff = qex - qapp
        new_err = np.sqrt( np.dot(diff,diff) ) / np.sqrt( np.dot(qex,qex) )

        r1 = '    Error = %(new).3e  ' % {'old': old_err, 'new' : new_err}
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 2) } 
        else:
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 

        print result

        old_err = new_err

        i = i + 1

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
