#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
import math

run_matlab_template = """
addpath('%(dogpack_dir)s/viz/matlab/');
run_convergence2( %(num_pts)i, '%(curr_folder)s' );
"""

def main( ):
    '''Write some help documentation here
'''
    input_file = 'tmpfile.m'
    dogpack_home_dir = os.environ['DOGPACK']
    t_order    = 2
    s_order    = 4
    number_pts = 1
    iterations = 8

    my_dictionary = {'dogpack_dir' : dogpack_home_dir, 'num_pts' : number_pts ,
        's_order' : s_order, 't_order' : t_order }
    print "# leading comments can be given a '#' character"
    print "# running t_order = ", t_order, " space_order = ", s_order
    old_err = i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i
        folder = 'outputSL%(s_order)i_%(t_order)i_00%(dir_num)i' % my_dictionary
        if( not os.path.exists(folder) ):
            print 'Did Not find folder: ' + folder
            break

        my_dictionary['curr_folder'] = folder
        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = i
        with closing(open(input_file,'w')) as data:
            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything
            print >> data, run_matlab_template % my_dictionary
        
        s = Popen('matlab -nodesktop -nosplash < ' + input_file, shell=True, stdout=PIPE).communicate()[0]
        new_err = float( s[s.find('>>')-1:len(s)].replace(">>","") )

        r1 = '%(old)e    %(new)e    ' % {'old': old_err, 'new' : new_err}

        if( new_err > 0 and old_err > 0 ):
            rat = math.log( old_err / new_err, 2 )
        else:
            rat = 1.0

        result = r1 + folder + '   log2( ratio ) = %(rat)8e' % \
            {'old' : old_err, 'new' : new_err, 'rat' : rat }
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
