#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

def main( ):
    '''Write some help documentation here
'''

    print "# leading comments can be given a '#' character"
    my_dictionary = {}
    old_err = i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        #folder = (os.getcwd() + '/outputSL4_4_%(dir_num)03i/') % my_dictionary
        folder = (os.getcwd() + '/output_%(dir_num)03i/') % my_dictionary

        print folder
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

        qex  = np.loadtxt(folder + "/q0000.dat")[1:]
        qapp = np.loadtxt(folder + "/q0010.dat")[1:]
      
        diff = qex - qapp

        new_err = np.sqrt( np.dot(diff,diff) ) / np.sqrt( np.dot(qex,qex) )

        r1 = '%(new).3e  ' % {'old': old_err, 'new' : new_err}
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 2) } 
        else:
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 

 # This is exactly the format I want:
#{\normalsize $25$}   &   {\normalsize $1.747\times 10^{-4}$}   & {\normalsize --} &  {\normalsize $8.292\times 10^{-5}$} & {\normalsize --}  \\
       
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
