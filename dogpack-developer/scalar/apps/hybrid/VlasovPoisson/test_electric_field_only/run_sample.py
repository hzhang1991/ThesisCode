#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

from params_template import dogpack_data_template, input2d_template

def main(space_order, iterations, mx_start, dt_start, h0_start ):
    '''Runs time stepping method given by ts_method for 'iterations' number of
    times.  The resolution gets increased by a factor of 'mx_ratio' for each
    iteration.  The starting resolution is given by 'mx_start'.   Output is
    written to folders 'output0000n' for n=0...iterations-1.
'''

    data_file = 'parameters.ini'
    ratio = 2

    my_start = mx_start
    time_order = space_order

    for i in range(iterations):

        mx_now = mx_start * ratio**i
        my_now = my_start * ratio**i
        dt_now = dt_start / ratio**i
        h0_now = h0_start / ratio**i

        my_dictionary = {'s_order' : space_order, 't_order' : time_order, \
                    'mx' : mx_now, 'my' : my_now, 'dt' : dt_now, "i_now": i }
        my_dictionary['h0'] = h0_now

#       # Create new mesh:
#       with closing(open('Unstructured_Mesh/input2D.data','w')) as data:
#           # print >> data, dogpack_data_template % locals() 
#           ## if we had used same names for everything
#           print >> data, input2d_template % my_dictionary
#       cmd = 'cd Unstructured_Mesh; ./mesh.exe'
#       call(cmd, shell=True)

        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        with closing(open(data_file,'w')) as data:
            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything
            print >> data, dogpack_data_template % my_dictionary
        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        cmd = './dog.exe -o output%(i_now)03i' % my_dictionary
        print cmd
        call(cmd, shell=True)
        print ''' 
//////////////////////// finished running output directory output%03i //////////
''' % (i)

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h | [-i ITERATIONS] [-x MX_START] [-w MY_START] 
                        [-z H0_START] [-t DT_START] [-s SPACE_ORDER] )
    
%s''' % main.__doc__)

    parser.add_option('-i', '--iterations', type='int', default=6, 
                       help='''number of times we wish to run computation.
                       Default = 6''')
    parser.add_option('-x', '--mx-start', type='int', default=4, 
                       help='''MX_START = number of grid points to start
                       computation with. Default = 4.''')
    parser.add_option('-t', '--dt-start', type='float', default=0.44,
                       help='''DT_START = starting dt value''')
    parser.add_option('-s', '--space-order', type='int', default=3,
                       help='''SPACE_ORDER = spatial order of accuracy. default = 3.''')
    parser.add_option('-z', '--h0-start', type='float', default=0.32,
                       help='''H0_START = starting value of h0.  Default = 0.32.''')
    opts, args = parser.parse_args()

#    check that the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')
    main(opts.space_order, opts.iterations, opts.mx_start,
        opts.dt_start, opts.h0_start)
