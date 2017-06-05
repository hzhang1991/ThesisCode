## ------------------------------------------------------------------------- ##
def plotdog1(points_per_dir=1, output_directory="output", point_type=1, q_name="a", aux_name="a", plotq1_name="plotq1" ):
    """Generic code for 1D plotting in DoGPack with matplotlib.

Execute via

    $ python $DOGPACK/viz/python/plotdog1.py

from an application directory.   For help type

    $ python $DOGPACK/viz/python/plotdog1.py -h
 
to see a list of options.

Parameters:
----------

    points_per_dir   = points per direction (spatial dimension)

    output_directory = location of output directory

    point_type       = 1:   uniform points on each element  ( Default )
                     = 2:   Gauss-Legendre points on each element

    q_name           = name of variable filename (i.e., [q_name]0000.dat, [q_name]0001.dat, ...)

    aux_name         = name of aux variable filename (i.e., [aux_name]0000.dat, [aux_name]0001.dat, ...)
    
    plotq1_name      = name of file for additional ploting options (i.e.,
                       filename is [plotq1_name].py).  Given that each
                       application needs to define limits for axes, color
                       options, etc., this is the single hook an application has
                       to set additional options. Typically, this is done
                       through a single call to plotq1.py.

See also: plotdog2.py and plotdog3.py for 2- and 3D plotting routines.
"""

    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from helper1 import SampleBasis1
    from helper1 import read_qfile

    from parse_parameters import parse_ini_parameters       #.ini parser
    
    TF = os.path.exists( output_directory )
    if TF==False:
        print("\n    Directory not found, output_directory = %s\n" % output_directory )
        exit()

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
    plotq1_file  = os.path.abspath("".join((plotq1_name,".py")))
    local_plotq1 = os.path.exists(plotq1_file)
    if local_plotq1==False:
        from plotq1_default import plotq1
    else:
        #from plotq1 import plotq1
        plotq1 = getattr(__import__(plotq1_name, fromlist=["plotq1"]), "plotq1")

    ini_params   = parse_ini_parameters( output_directory+'/parameters.ini' )
    #print('ini_params = ', ini_params)
    GridType = ini_params['mesh_type']
    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    nplot    = ini_params['nout']
    space_order    = ini_params['space_order']
    datafmt  = ini_params['datafmt']
    mx       = ini_params['mx']
    xlow     = ini_params['xlow']
    xhigh    = ini_params['xhigh']


#   params = np.zeros(8,  dtype=np.float64) 
#   GridType = read_params(output_directory,params,qhelpname)    
#   meqn    = int(params[0])
#   maux    = int(params[1])
#   nplot   = int(params[2])
#   space_order   = int(params[3])
#   datafmt = int(params[4])
#   mx      = int(params[5])
#   xlow    = params[6]
#   xhigh   = params[7]

    print("")
    print("         GridType = %s"  % GridType          )
    print("   points_per_dir = %s"  % points_per_dir    )
    print("       point_type = %s"  % point_type        )
    print(" output_directory = %s " % output_directory  )
    print("           q_name = %s"  % q_name            )
    print("         aux_name = %s"  % aux_name          )
    print("      plotq1_name = %s"  % plotq1_name       )
    print("")

    # Grid information
    mx_old  = mx
    mx      = mx*points_per_dir
    dx_old  = (xhigh-xlow)/np.float64(mx_old)
    dx      = (xhigh-xlow)/np.float64(mx)

    # Construct point values on entire grid
    if( point_type == 1 ):

        # Location of points on Canonical element
        ds  = 2.0 / (1 + points_per_dir )
        
        s1d = np.linspace( -1.+ds, 1.-ds, num=points_per_dir )
        
        # Cell centers
        xc = np.linspace( xlow+0.5*dx, xhigh-0.5*dx, mx )

        xc_old = np.zeros(mx_old, dtype=np.float64 )
        xc_old[0] = xlow + 0.5*dx_old
        for i in range(1, mx_old):
            xc_old[i] = xc_old[i-1] + dx_old

    else:

        # Accessors
        sq3, sq5, sq7 = ( np.sqrt(3.), np.sqrt(5.), np.sqrt(7.) )
   
        # quadrautre points (TODO - add in 6th order story ... )
        if(points_per_dir==1):
            s1d = np.array( [0.0] )
        elif(points_per_dir==2):
            s1d = np.array( [-1.0/sq3, 1.0/sq3] )
        elif(points_per_dir==3):
            s1d = np.array( [-sq3/sq5, 0.0e0, sq3/sq5] )
        elif(points_per_dir==4):
            s1d = np.array( [-sqrt(3.0+sqrt(4.8))/sq7, -sqrt(3.0-sqrt(4.8))/sq7, sqrt(3.0-sqrt(4.8))/sq7,  sqrt(3.0+sqrt(4.8))/sq7] )
        elif(points_per_dir==5):
            s1d = np.array( [-np.sqrt(5.0 + np.sqrt(40.0/7.0))/3.0,
                   -np.sqrt(5.0 - np.sqrt(40.0/7.0))/3.0,
                   0.0,
                   np.sqrt(5.0 - np.sqrt(40.0/7.0))/3.0,
                   np.sqrt(5.0 + np.sqrt(40.0/7.0))/3.0] )

        xc    = np.zeros( mx_old*points_per_dir, dtype=np.float64 )
        xline = np.linspace( xlow+0.5*dx, xhigh-0.5*dx, mx_old )
        xline = np.linspace( xlow, xhigh, mx_old+1 )

        kk = 0
        for i in range( mx_old ):
            for k3 in range( kk, kk+points_per_dir ):
                xc[ k3 ] = xline[i]+(dx_old/2.)*(s1d[k3-kk]+1.0)
            kk = kk + points_per_dir

    ## --------------------------------------------------------------------- ##
    # Sample basis functions on mesh
    # size of phi = (points_per_dir, space_order)
    ## --------------------------------------------------------------------- ##
    phi = SampleBasis1(s1d, space_order )

    quit = -1
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print("")
    if (not m):
        m = 1
    else:
        m = int(m)

    if(m<1):
        print("")
        print("  Error, need m > 1,  m = %d" % m )
        print("")
        exit(1)
    elif(m>meqn):   
        print("")
        print("  Error, need m <= %d, m = %d " % (meqn, m ) )
        print("")
        exit(1)

    nf = 0     # Frame number
    n1 = -1    # Frame number

    plt.ion()
    while (nf!=-1):
        tmp1 = "".join((" Plot which frame ( 0 - ",str(nplot)))
        tmp2 = "".join((tmp1," ) [type -1 or type 'quit' to quit] ? "))
        nf   = raw_input(tmp2)
        if(not nf):
            n1 = n1 + 1
            nf = 0
        elif(nf=="quit" or nf=="q"):
            nf = -1
        else:
            nf = int(nf)
            n1 = nf

        if( n1 > nplot ):
            print("")
            print(" End of plots " )
            print("")
            n1 = nplot

        if( nf != -1 ):

            ## ------------------------------------------------------------------------- ##
            # Solution -- q
            # solution should be found in file
            #     output_directory/q[n1].dat
            ## ------------------------------------------------------------------------- ##


            # Q-file name (with directory thrown in the mix)
            qfile   = output_directory + "/" + q_name + "%04d" % n1 + ".dat"

            # Read in the coefficients from file
            mtmp       = mx_old*meqn*space_order
            time, qtmp = read_qfile(mtmp, qfile )
            qcoeffs    = np.reshape(qtmp,(space_order,meqn,mx_old))

            # Evaluate the basis functions
            qsoln = np.zeros((mx,meqn), dtype=np.float64)
            v1 = np.zeros(space_order, dtype=np.float64)
            v2 = np.zeros(space_order, dtype=np.float64)
            for i in range(1,mx_old+1):
                for me in range(1,meqn+1):
                    for ii in range(1,points_per_dir+1):
                        v1[:] = phi[ii-1,:]
                        v2[:] = qcoeffs[:,me-1,i-1]
                        tmp = 0.0
                        for k in range(0,space_order):
                            tmp = tmp + v1[k]*v2[k]
                        qsoln[(i-1)*points_per_dir+ii-1,me-1] = tmp

            ## ------------------------------------------------------------------------- ##
            # Aux arrays -- aux
            # solution should be found in file
            #     output_directory/aux_name[n1].dat
            ## ------------------------------------------------------------------------- ##
            if( maux > 0):

                # Auxiliary file name (with directory thrown in the mix)
                afile   = output_directory + "/" + aux_name + "%04d" % n1 + ".dat"

                # Read in the coefficients from file
                mtmp       = mx_old*maux*space_order
                time, atmp = read_qfile(mtmp, afile)
                acoeffs    = np.reshape(atmp, (space_order,maux,mx_old) )

                # Evaluate the coefficients
                auxsoln = np.zeros((mx,maux), dtype=np.float64)
                v1 = np.zeros(space_order, dtype=np.float64)
                v2 = np.zeros(space_order, dtype=np.float64)
                for i in range(1,mx_old+1):
                    for me in range(1,maux+1):
                        for ii in range(1,points_per_dir+1):
                            v1[:] = phi[ii-1,:]
                            v2[:] = acoeffs[:,me-1,i-1]
                            tmp = 0.0
                            for k in range(0,space_order):
                                tmp = tmp + v1[k]*v2[k]
                            auxsoln[(i-1)*points_per_dir+ii-1,me-1] = tmp
            else:
                auxsoln = False

            ###################################################################
            #
            # USER SUPPLIED FUNCTION (or default function )
            #
            # This is your one single hook into the code for adding extra
            # information you would like to plot.
            #
            ###################################################################
            plotq1(m-1, space_order, meqn, mx, time, xc, qsoln, auxsoln)

    plt.ioff()
    print("")
## ------------------------------------------------------------------------- ##
    

## ------------------------------------------------------------------------- ##
def parse_input( help_message ):
    """Parse command line arguments for 1D plotting routines."""

    import argparse, sys

    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    parser.add_argument('-p', '--points-per-dir', 
        type    = int, 
        default = 1, 
        help    =
'''Number of points per cell to be plotted.
(default: 1)''')

    parser.add_argument('-o', '--output-directory', 
        type        = str, 
        default     = 'output', 
        help        =
'''Location of the output directory where coefficients can be located.
(default: output)''')

    parser.add_argument('-t', '--point-type', 
        type    = int, 
        default = 1, 
        help    =
'''Type of points to be plotted.  
    point-type==1: Linearly spaced.
    point-type==2: Gauss-Legendre.
(default: 1)''')

    parser.add_argument( '-q', '--q-name', type = str, default='q', 
        help        =
'''Name of variable used in the filename.  In most routines this is 'q' but in some cases
one may wish to save additional data.  For example, a 2D code may want to save
1D data objects, and then plot them.
(default: q).''')

    parser.add_argument( '-x', '--aux-name', type = str, default='a', 
        help        =
'''Name of aux variable used in the filename.  In most routines this is 'a' but in some cases
one may wish to save additional data.  For example, a 2D code may want to save
1D data objects, and then plot them.
(default: a).''')

    parser.add_argument( '-a', '--plotq1-name', 
        type = str, 
        default     = 'plotq1', 
        help        =
'''filename for file with additional plotting options.
(default: plotq1)''')

    return parser.parse_args()
#----------------------------------------------------------

if __name__== '__main__':
    """Main python plotting routine for 1D simulations in DoGPack.

    When run from the command line, this script parses user supplied command
    line arguments, if any, and then executes plotdog1.  To see a list of
    options, type 

        python $DOGPACK/viz/python/plotq1.py -h

    """

    # Parse input arguments
    args = parse_input( plotdog1.__doc__ )
    #print(args)
    #print('')

    # Call the main 1D plotting routine
    plotdog1(args.points_per_dir, args.output_directory, args.point_type, args.q_name, args.aux_name, args.plotq1_name )
