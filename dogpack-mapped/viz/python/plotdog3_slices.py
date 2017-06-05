def plotdog3_slices(points_per_dir_in="1",outputdir="output"):
    """
Generic code for plotting DoGPack output in matplotlib.
    """

    """
Execute via

    $ python $DOGPACK/viz/python/plotdog3_slices.py

from an application directory.   For help type

    $ python $DOGPACK/viz/python/plotdog3_slices.py -h
 
to see a list of options.

    Usage: plotdog3_slices( points_per_dir_in="1", outputdir="output")
 
   points_per_dir = points per direction (spatial dimension)
 
   outputdir = location of output directory
   
"""

    import os
    import sys
    from math import sqrt
    import numpy as np
    import matplotlib.pyplot as plt
    from helper3 import get_grid_type
    from helper3 import read_params
    from helper3 import get_kmax
    from helper3 import GetCart2Legendre
    from helper3 import read_qfile
    from helper3 import sample_state2_cart_mod
    from helper3 import get_numslices
    from helper3 import read_slice_info
    
    points_per_dir = int(points_per_dir_in)
    point_type = 1

    TF = os.path.exists(outputdir)
    if TF==False:
        print ""
        print "    Directory not found, outputdir =",outputdir
        print ""
        return -1
    else:
        print "found outputdir = ", outputdir
    
    GridType = get_grid_type(outputdir)
    
    if (GridType=="Cartesian"):
        params = np.zeros(14, float)
        ndims = read_params(outputdir,params)
        meqn    = int(params[0])
        maux    = int(params[1])
        nplot   = int(params[2])
        meth1   = int(params[3])
        datafmt = int(params[4])
        mx      = int(params[5])
        my      = int(params[6])
        mz      = int(params[7])
        xlow    = params[8]
        xhigh   = params[9]
        ylow    = params[10]
        yhigh   = params[11]
        zlow    = params[12]
        zhigh   = params[13]

        numslices = get_numslices(outputdir)
        slicekind = np.zeros(numslices,int)
        sliceidx  = np.zeros(numslices,int)
        read_slice_info(outputdir,slicekind,sliceidx)
        
    elif (GridType=="Unstructured"):
        print " In 3D currently only GridType=Cartesian is supported. GridType = ",GridType
        print " "
        return -1

    if numslices<1:
        print " "
        print " NOTE: this code won't do anything because numslices = ",numslices
        print " "
        return -1
    
    kmax = get_kmax(meth1,2)
        
    print ""
    print "        GridType = ",GridType
    print "  points_per_dir = ",points_per_dir
    print "      point_type = ",point_type
    print "       outputdir = ",outputdir
    for ns in range(0,numslices):
        print "  slicekind(",ns+1,") = ",slicekind[ns],", sliceidx(",ns+1,") = ",sliceidx[ns]
    print ""

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
        
    if (GridType=="Cartesian"):
        plotq3_file  = os.path.abspath("plotq3_cart_slices.py")
        local_plotq3 = os.path.exists(plotq3_file)
        if local_plotq3==False:
            from plotq3_cart_slices_default import plotq3_cart_slices
        else:
            from plotq3_cart_slices import plotq3_cart_slices

        # Grid information
        mx_old = mx
        my_old = my
        mz_old = mz
        dx_old = (xhigh-xlow)/mx_old
        dy_old = (yhigh-ylow)/my_old
        dz_old = (zhigh-zlow)/mz_old
        mx = mx*points_per_dir
        my = my*points_per_dir
        mz = mz*points_per_dir
        dx = (xhigh-xlow)/mx
        dy = (yhigh-ylow)/my
        dz = (zhigh-zlow)/mz

        xc = np.zeros((mx,my,mz),float)
        yc = np.zeros((mx,my,mz),float)
        zc = np.zeros((mx,my,mz),float)
        for i in range(0,mx):
            for j in range(0,my):
                for k in range(0,mz):
                    xc[i,j,k] = xlow + (i+0.5)*dx
                    yc[i,j,k] = ylow + (j+0.5)*dy
                    zc[i,j,k] = zlow + (k+0.5)*dz

        xl = np.zeros((mx+1,my+1,mz+1),float)
        yl = np.zeros((mx+1,my+1,mz+1),float)
        zl = np.zeros((mx+1,my+1,mz+1),float)
        for i in range(0,mx+1):
            for j in range(0,my+1):
                for k in range(0,mz+1):
                    xl[i,j,k] = xlow + i*dx
                    yl[i,j,k] = ylow + j*dy
                    zl[i,j,k] = zlow + k*dz

        xc_z = np.zeros((mx,my),float)
        yc_z = np.zeros((mx,my),float)
        for i in range(0,mx):
            for j in range(0,my):
                xc_z[i,j] = xc[i,j,0]
                yc_z[i,j] = yc[i,j,0]

        xc_y = np.zeros((mx,mz),float)
        zc_y = np.zeros((mx,mz),float)
        for i in range(0,mx):
            for k in range(0,mz):
                xc_y[i,k] = xc[i,0,k]
                zc_y[i,k] = zc[i,0,k]

        yc_x = np.zeros((my,mz),float)
        zc_x = np.zeros((my,mz),float)
        for j in range(0,my):
            for k in range(0,mz):
                yc_x[j,k] = yc[0,j,k]
                zc_x[j,k] = zc[0,j,k]
                
        xl_z = np.zeros((mx+1,my+1),float)
        yl_z = np.zeros((mx+1,my+1),float)
        for i in range(0,mx+1):
            for j in range(0,my+1):
                xl_z[i,j] = xl[i,j,0]
                yl_z[i,j] = yl[i,j,0]

        xl_y = np.zeros((mx+1,mz+1),float)
        zl_y = np.zeros((mx+1,mz+1),float)
        for i in range(0,mx+1):
            for k in range(0,mz+1):
                xl_y[i,k] = xl[i,0,k]
                zl_y[i,k] = zl[i,0,k]

        yl_x = np.zeros((my+1,mz+1),float)
        zl_x = np.zeros((my+1,mz+1),float)
        for j in range(0,my+1):
            for k in range(0,mz+1):
                yl_x[j,k] = yl[0,j,k]
                zl_x[j,k] = zl[0,j,k]
                
        #####################################################

        # 1D points 
        dxi = 1.0/float(points_per_dir)
        s1d = -1 + dxi + 2*dxi * np.arange(points_per_dir)

        kk=-1;
        s2d = np.zeros((points_per_dir*points_per_dir,2),float)
        for jj in range(0,points_per_dir):
            for ii in range(0,points_per_dir):
                kk = kk+1
                s2d[kk,0] = s1d[ii]
                s2d[kk,1] = s1d[jj]        

        # Sample Legendre polynomial on the midpoint of each sub-element
        p2 = points_per_dir*points_per_dir
        LegVals = np.zeros((kmax,p2),float)
        GetCart2Legendre(meth1,points_per_dir,s2d,LegVals)
  
        q=-1;
        tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
        tmp2 = "".join((tmp1," ) ? "))
        m = raw_input(tmp2)
        print ""
        if (not m):
            m = 1
        else:
            m = int(m)

        if m<1:
            print ""
            print "  Error, need m > 1,  m = ",m
            print ""
            return -1
        elif m>meqn:
            print ""
            print "  Error, need m <=",meqn,",  m = ",m
            print ""
            return -1

        kn = 0;
    
        n = 0;
        nf = 0;
        n1 = -1;

        plt.ion()
        while (nf!=-1):        
            tmp1 = "".join((" Plot which frame ( 0 - ",str(nplot)))
            tmp2 = "".join((tmp1," ) [type -1 or q to quit] ? "))
            nf = raw_input(tmp2)
            if (not nf):
                n1 = n1 + 1
                nf = 0
            elif nf=="q":
                nf = -1
            else:
                nf = int(nf)
                n1 = nf

            if n1>nplot:
                print ""
                print " End of plots "
                print ""
                n1 = nplot

            if (nf!=-1):
                # dimension each to be able to store each slice
                mmax = np.max([mx,my,mz])
                qaug = np.zeros((mmax+1,mmax+1,meqn,numslices),float)
                #if maux>0:
                #    aux_aug = np.zeros((mmax+1,mmax+1,maux,numslices),float)
                
                ms1 = np.zeros(numslices,int);
                ms2 = np.zeros(numslices,int);

                # Loop over each slice
                for ns in range(0,numslices):
                    # Get dimensions
                    if slicekind[ns]==1:                    
                        m1 = mx
                        m2 = my
                        m1_old = mx_old
                        m2_old = my_old
                    elif slicekind[ns]==2:
                        m1 = mx
                        m2 = mz
                        m1_old = mx_old
                        m2_old = mz_old
                    else:
                        m1 = my
                        m2 = mz
                        m1_old = my_old
                        m2_old = mz_old
          
                    ms1[ns] = m1
                    ms2[ns] = m2

                    qfile_tmp1 = "".join((str(n1+10000),"_slic"))
                    qfile_tmp2 = "q" + qfile_tmp1[1:]
                    qfile_tmp3 = "".join((str(ns+10001),".dat"))
                    qfile_tmp4 = "e" + qfile_tmp3[1:]
                    qfile_tmp5 = "".join((qfile_tmp2,qfile_tmp4))
                    qfile = "".join(("".join((outputdir,"/")),qfile_tmp5))

                    mtmp = m1_old*m2_old*meqn*kmax
                    qtmp = np.zeros(mtmp,float)   
                    time = read_qfile(mtmp,qfile,qtmp)
                    qcoeffs = np.reshape(qtmp,(kmax,meqn,m2_old,m1_old))
                
                    qsoln = np.zeros((m1*m2,meqn),float)
                    sample_state2_cart_mod(m1_old,m2_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln)
                    qaug[0:m1,0:m2,0:meqn,ns] = np.reshape(qsoln,(m1,m2,meqn),'F')            

                # USER SUPPLIED FUNCTION
                plotq3_cart_slices(outputdir,n1,
                                   numslices,slicekind,sliceidx,
                                   m-1,meth1,meqn,ms1,ms2,time,s2d,
                                   xlow,xhigh,ylow,yhigh,zlow,zhigh,
                                   mx,my,mz,dx,dy,dz,
                                   mx_old,my_old,mz_old,dx_old,dy_old,dz_old,
                                   points_per_dir,LegVals,
                                   xc_z,yc_z,xc_y,zc_y,yc_x,zc_x,
                                   xl_z,yl_z,xl_y,zl_y,yl_x,zl_x,qaug)
            
        plt.ioff()
        print ""
        
    else:
        
        print ""
        print " Error in plotdog3_slices.py: GridType = ",GridType," is not supported."
        print ""
        return -1
            
#----------------------------------------------------------

    

#----------------------------------------------------------
if __name__=='__main__':
    """
    If executed at command line prompt, simply run plotdog3_slices
    """

    import optparse
    parser = optparse.OptionParser(
        usage=''' %%prog (-h | [-p POINTS_PER_DIR] [-o OUTPUT_DIRECTORY] [-t POINT_TYPE] )
    
%s''' % plotdog3_slices.__doc__)

    parser.add_option('-p', '--points-per-dir', type='int', default=1, 
                       help='''number of points per cell to be plotted.
                       Default = 1''')
    parser.add_option('-o', '--output-directory', type='string', default='output', 
                       help='''OUTPUT_DIR = output directory where coefficients
                       can be located''')

    opts, args = parser.parse_args()
    plotdog3_slices(opts.points_per_dir, opts.output_directory )

#    This is how you could check if the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')

### The old way of calling plotdog3_slices: ###
#   import sys
#   args = sys.argv[1:]   # any command line arguments
#   plotdog3_slices(*args)

