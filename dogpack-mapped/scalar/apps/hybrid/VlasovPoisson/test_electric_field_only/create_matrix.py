#!/usr/bin/env python

def main(sorder, h0):

    from contextlib import closing
    from subprocess import call, Popen, PIPE

    my_dictionary = {'h0':h0, 'sorder':sorder}

    # input2D.data template:
    from params_template import input2d_template, matlab_script_template

    # Create new mesh:
    with closing(open('Unstructured_Mesh/input2D.data','w')) as data:
        print >> data, input2d_template % my_dictionary
    cmd = 'cd Unstructured_Mesh; make; ./mesh.exe'
    call(cmd, shell=True)

    ## Create script file to call matlab with:
    # Create new mesh:
    with closing(open('matlab/script.m','w')) as data:
        print >> data, matlab_script_template % my_dictionary
    cmd = 'cd matlab; matlab -nodesktop -nosplash -nojvm -r script'
    call(cmd, shell=True)

if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser( description='''Create a new mesh, call the
    matalab code responsible for setting up matrices for Poisson solver.''')

    parser.add_argument('-s', '--space_order', type=int, default=3, \
        help="Spatial order of accuracy.  (default=3)" )
    parser.add_argument('-h0', '--h0', type=float, default=0.20, \
        help="Initial grid spacing for unstructured mesh.  (default=0.20)" )

    args = parser.parse_args()

    main( args.space_order, args.h0 )
