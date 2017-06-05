def plotdog2(points_per_dir_in="1",outputdir="output",point_type_in="1"):
    """
Generic code for plotting DoGPack output in matplotlib.
    """

    """
Execute via

    $ python $DOGPACK/viz/python/plotdog2.py

from an application directory.   For help type

    $ python $DOGPACK/viz/python/plotdog2.py -h
 
to see a list of options.

    Usage: plotdog2( points_per_dir_in="1", outputdir="output", point_type_in="1")
 
   points_per_dir = points per direction (spatial dimension)
 
   outputdir = location of output directory

   point_type = 1:   uniform points on each element
              = 2:   Gauss-Legendre points on each element
"""
    def phi(xin,yin,x1,y1,dx1,dy1,a):
	 c0=0.0
	 c1=0.0
	 c2=0.0
	 c3=0.0
	 c4=0.0
	 c5=0.0

	 x=2.0*(xin-x1)/dx1-1.0;
	 y=2.0*(yin-y1)/dy1-1.0;

	 if(a==1):
	    c0=-1.0/6.0;
	    c1=-1.0/6.0;
	    c2=-1.0/6.0;
	    c3=1.0/3.0;
	    c4=1.0/6.0;
	    c5=1.0/3.0;
	   
	 if(a==2):
	    c0=1.0;
	    c1=0.0;
	    c2=-1.0/3.0;
	    c3=0.0;
	    c4=-2.0/3.0;
	    c5=-2.0/3.0;
	    
	 if(a==3):
	    c0=-1.0/6.0;
	    c1=1.0/6.0;
	    c2=-1.0/6.0;
	    c3=-1.0/3.0;
	    c4=1.0/6.0;
	    c5=1.0/3.0;
	    
	 if(a==4):
	    c0=1.0/6.0;
	    c1=1.0/3.0;
	    c2=1.0/6.0;
	    c3=1.0/3.0;
	    c4=1.0/3.0;
	    c5=-1.0/3.0;
	    
	 if(a==5):
	    c0=0.0;
	    c1=0.0;
	    c2=1.0/3.0;
	    c3=0.0;
	    c4=-1.0/3.0;
	    c5=2.0/3.0;
	    
	 if(a==6):
	    c0=1.0/6.0;
	    c1=-1.0/3.0;
	    c2=1.0/6.0;
	    c3=-1.0/3.0;
	    c4=1.0/3.0;
	    c5=-1.0/3.0;
	    
	  
	 phi1=c0+c1*x+c2*y+c3*x*y+c4*x*x+c5*y*y;

	 return phi1;

    import os
    import sys
    from math import sqrt
    import matplotlib.pyplot as plt
    from helpermap2 import get_grid_type
    from helpermap2 import read_params
    from helpermap2 import get_kmax
    from helpermap2 import GetCart2Legendre
    from helpermap2 import read_qfile
    from helpermap2 import sample_state2_cart_mod_map,sample_state2_cart_mod
    from helpermap2 import read_mesh_params
    from helpermap2 import read_node
    from helpermap2 import read_tnode
    from helpermap2 import GetMonomialToLegendre
    from helpermap2 import GetUnstLegendre
    from helpermap2 import set_tcounter
    from helpermap2 import set_soln_at_node_values
    from helpermap2 import sample_state2_unst
    from helpermap2 import DivideUnstMesh

    import numpy as np

    try:
     from mapc2p import mapc2p # Import DNSPython
    except ImportError:
     def mapc2p(xc,yc):
         """
      Specifies the mapping to curvilinear coordinates -- should be consistent
      with mapc2p.f
         """
         return xc,yc



    
    points_per_dir = int(points_per_dir_in)
    point_type = int(point_type_in)

    TF = os.path.exists(outputdir)
    if TF==False:
        print ""
        print "    Directory not found, outputdir =",outputdir
        print ""
        return -1
    else:
        print "found outputdir = ", outputdir

    if (point_type!=1 and point_type!=2):
        point_type = 1
    
    if (point_type==2 and points_per_dir>5):
        print 'setting points_per_dir = 5 instead of requested %d' % points_per_dir
        points_per_dir = 5

    
    GridType = get_grid_type(outputdir)
    if (GridType=="Cartesian"):
        params = np.zeros(11, float)
        read_params(outputdir,params)
        meqn    = int(params[0])
        maux    = int(params[1])
        nplot   = int(params[2])
        meth1   = int(params[3])
        mx      = int(params[5])
        my      = int(params[6])
        xlow    = params[7]
        xhigh   = params[8]
        ylow    = params[9]
        yhigh   = params[10]
        datafmt = int(params[4])
    elif (GridType=="Unstructured"):
        params = np.zeros(4, int)
        read_params(outputdir,params)
        meqn    = params[0]
        maux    = params[1]
        nplot   = params[2]
        meth1   = params[3]
    kmax = get_kmax(meth1,2)
        
    print ""
    print "        GridType = ",GridType
    print "  points_per_dir = ",points_per_dir
    print "      point_type = ",point_type
    print "       outputdir = ",outputdir
    print ""

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
        
    if (GridType=="Cartesian"):
        plotq2_file  = os.path.abspath("plotq2_cart.py")
        local_plotq2 = os.path.exists(plotq2_file)
        if local_plotq2==False:
            from plotq2_cart_default import plotq2_cart
        else:
            from plotq2_cart import plotq2_cart

        # Grid information
        mx_old = mx
        my_old = my
        dx_old = (xhigh-xlow)/mx_old
        dy_old = (yhigh-ylow)/my_old
        mx = mx*points_per_dir
        my = my*points_per_dir
        dx = (xhigh-xlow)/mx
        dy = (yhigh-ylow)/my
        print "mx=",mx,"my=",my
        if (point_type==1):
            xl = np.zeros((mx+1,my+1),float)
            yl = np.zeros((mx+1,my+1),float)
            phil=np.zeros((mx,my,points_per_dir,kmax),float)
            xc = np.zeros((mx,my),float)
            yc = np.zeros((mx,my),float)


            for j in range(0,my+1):
                xl[:,j] = xlow + dx*np.arange(mx+1)[:]
            for i in range(0,mx+1):
                yl[i,:] = ylow + dy*np.arange(my+1)[:]
            for j in range(0,my):
                xc[:,j] = (xlow+0.5*dx) + dx*np.arange(mx)[:]
            for i in range(0,mx):
                yc[i,:] = (ylow+0.5*dy) + dy*np.arange(my)[:]

            xpl,ypl=mapc2p(xl,yl)

            for i in range(0,mx):
             for j in range(0,my):

              xp1 = xpl[i,j];
              yp1=  ypl[i,j];
              xp2 = xpl[i+1,j];
              yp2=  ypl[i+1,j];
              xp3 = xpl[i+1,j+1];
              yp3=  ypl[i+1,j+1];
              xp4 = xpl[i,j+1];
              yp4=  ypl[i,j+1];

              xmax=max(max(xp1,xp2),max(xp3,xp4));xmin=min(min(xp1,xp2),min(xp3,xp4));
              ymax=max(max(yp1,yp2),max(yp3,yp4));ymin=min(min(yp1,yp2),min(yp3,yp4));
              dx1=xmax-xmin;dy1=ymax-ymin;
              xin=0.5*(xmax+xmin);
              yin=0.5*(ymax+ymin);
              for l in range(0,kmax):
                phil[i,j,0,l]=phi(xin,yin,xmin,ymin,dx1,dy1,l+1)

              

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
                    
        else:

            sq3 = sqrt(3.0)
            sq5 = sqrt(5.0)
            sq7 = sqrt(7.0)

            # 1D quadrature points
            s1d = np.zeros(points_per_dir,float)
            if (points_per_dir==1):
                s1d[0] = 0.0
            elif (points_per_dir==2):
                s1d[0] = -1.0/sq3
                s1d[1] =  1.0/sq3
            elif (points_per_dir==3):
                s1d[0] = -sq3/sq5
                s1d[1] =  0.0
                s1d[2] =  sq3/sq5
            elif (points_per_dir==4):
                s1d[0] = -sqrt(3.0+sqrt(4.8))/sq7
                s1d[1] = -sqrt(3.0-sqrt(4.8))/sq7
                s1d[2] =  sqrt(3.0-sqrt(4.8))/sq7
                s1d[3] =  sqrt(3.0+sqrt(4.8))/sq7
            elif (points_per_dir==5):
                s1d[0] = -sqrt(5.0 + sqrt(40.0/7.0))/3.0
                s1d[1] = -sqrt(5.0 - sqrt(40.0/7.0))/3.0
                s1d[2] =  0.0
                s1d[3] =  sqrt(5.0 - sqrt(40.0/7.0))/3.0
                s1d[4] =  sqrt(5.0 + sqrt(40.0/7.0))/3.0

            kk=-1;
            s2d = np.zeros((points_per_dir*points_per_dir,2),float)
            for jj in range(0,points_per_dir):
                for ii in range(0,points_per_dir):
                    kk = kk+1
                    s2d[kk,0] = s1d[ii]
                    s2d[kk,1] = s1d[jj]

            xx = np.zeros(mx,float)
            yy = np.zeros(my,float)

            kk=0
            xtmp = xlow-0.5*dx_old
            for i in range(0,mx_old):
                xtmp = xtmp + dx_old
                for m in range(0,points_per_dir):
                    xx[kk+m] = xtmp + 0.5*dx_old*s1d[m]
                kk = kk + points_per_dir

            kk=0
            ytmp = ylow-0.5*dy_old
            for j in range(0,my_old):
                ytmp = ytmp + dy_old
                for m in range(0,points_per_dir):
                    yy[kk+m] = ytmp + 0.5*dy_old*s1d[m]
                kk = kk + points_per_dir

            xc = np.zeros((mx,my),float)
            yc = np.zeros((mx,my),float)

            for i in range(0,mx):
                for j in range(0,my):
                    xc[i,j] = xx[i]
                    yc[i,j] = yy[j]

            xxx = np.zeros(mx+1,float)
            yyy = np.zeros(my+1,float)

            xxx[0] = xlow
            for i in range(1,mx):
                xxx[i] = 0.5*(xx[i]+xx[i-1])
            xxx[mx] = xhigh

            yyy[0] = ylow
            for j in range(1,my):
                yyy[j] = 0.5*(yy[j]+yy[j-1])
            yyy[my] = yhigh

            xl = np.zeros((mx+1,my+1),float)
            yl = np.zeros((mx+1,my+1),float)

            for i in range(0,mx+1):
                for j in range(0,my+1):
                    xl[i,j] = xxx[i]
                    yl[i,j] = yyy[j]

        # Sample Legendre polynomial on the midpoint of each sub-element
        p2 = points_per_dir*points_per_dir
        LegVals = np.zeros((kmax,p2),float)
        GetCart2Legendre(meth1,points_per_dir,s2d,LegVals)
  
        q=-1;
        tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
        tmp2 = "".join((tmp1," ) ? "))
        m = raw_input(tmp2)
        print "herro?",m
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
        #while (nf!=-1): 
        for n1 in range(nplot+1): 
            """      
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
            """
            if n1>nplot:
                print ""
                print " End of plots "
                print ""
                n1 = nplot

            if (nf!=-1):
                # Solution -- q
                # solution should be found in file
                #     outputdir/q[n1].dat

                qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
                qfile_tmp = "q" + qfile_tmp_tmp[1:]
                qfile = "".join(("".join((outputdir,"/")),qfile_tmp))

                mtmp = mx_old*my_old*meqn*kmax
                qtmp = np.zeros(mtmp,float)   
                time = read_qfile(mtmp,qfile,qtmp)
                qcoeffs = np.reshape(qtmp,(kmax,meqn,my_old,mx_old))
                
                qsoln = np.zeros((mx*my,meqn),float)
                sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,phil,qsoln)
                qsoln = np.reshape(qsoln,(mx,my,meqn),'F')

                # USER SUPPLIED FUNCTION
                print xc.shape,xl.shape
                #plotq2_cart(m-1,meth1,meqn,mx,my,time,xc,yc,xl,yl,qsoln)
                plotq2_cart(outputdir,n1,
                            m-1,meth1,meqn,time,
                            points_per_dir,LegVals,
                            xlow,xhigh,ylow,yhigh,
                            mx,my,dx,dy,
                            mx_old,my_old,dx_old,dy_old,
                            xc,yc,xl,yl,qsoln);
            
        plt.ioff()
        print ""

    elif (GridType=="Unstructured"):
        plotq2_file  = os.path.abspath("plotq2_unst.py")
        local_plotq2 = os.path.exists(plotq2_file)
        if local_plotq2==False:
            from plotq2_unst_default import plotq2_unst
        else:
            from plotq2_unst import plotq2_unst

        print " Creating mesh ... "

        # READ-IN MESH INFO
        meshdir = "".join((outputdir,"/mesh_output"))
        TF = os.path.exists(meshdir)
        if TF==False:
            print ""
            print "    Directory not found, meshdir =",meshdir
            print ""
            return -1

        mesh_params = np.zeros(7,int)
        read_mesh_params(meshdir,mesh_params)
        NumElems      = mesh_params[0]
        NumPhysElems  = mesh_params[1]
        NumGhostElems = mesh_params[2]
        NumNodes      = mesh_params[3]
        NumPhysNodes  = mesh_params[4]
        NumBndNodes   = mesh_params[5]
        NumEdges      = mesh_params[6]

        tnode = np.zeros((NumPhysElems,3),int)
        x  = np.zeros(NumPhysNodes,float)
        y  = np.zeros(NumPhysNodes,float)
        read_node(meshdir,NumPhysNodes,x,y)
        read_tnode(meshdir,NumPhysElems,tnode)

        xlow  = min(x)
        ylow  = min(y)
        xhigh = max(x)
        yhigh = max(y)

        tcounter = np.zeros(NumPhysNodes,int)
        set_tcounter(NumPhysElems,tnode,tcounter)

        # Add extra points and elements if points_per_dir>1
        p2 = points_per_dir*points_per_dir
        points_per_elem = ((points_per_dir+1)*(points_per_dir+2))/2
        zx = np.zeros(p2,float)
        zy = np.zeros(p2,float)
        
        if (points_per_dir>1):
            x_tmp = np.zeros(NumPhysElems*points_per_elem,float)
            y_tmp = np.zeros(NumPhysElems*points_per_elem,float)
            tnode_new = np.zeros((NumPhysElems*p2,3),int)
            NumPhysElems_new = 0
            NumPhysNodes_new = 0
            newsizes = np.zeros(2,int)            
            DivideUnstMesh(points_per_dir,NumPhysElems,x,y,tnode,
                           x_tmp,y_tmp,tnode_new,zx,zy,newsizes)
            NumPhysElems_new = newsizes[0]
            NumPhysNodes_new = newsizes[1]
            x_new = np.zeros(NumPhysNodes_new,float)
            y_new = np.zeros(NumPhysNodes_new,float)
            for ijk in range(0,NumPhysNodes_new):
                x_new[ijk] = x_tmp[ijk]
                y_new[ijk] = y_tmp[ijk]
            tcounter_new = np.zeros(NumPhysNodes_new,int)
            set_tcounter(NumPhysElems_new,tnode_new,tcounter_new)
        else:
            NumPhysElems_new = NumPhysElems
            NumPhysNodes_new = NumPhysNodes
            x_new = x
            y_new = y
            tnode_new = tnode
            tcounter_new = tcounter

        # Get physical midpoints of each element
        xmid = np.zeros(NumPhysElems_new,float)
        ymid = np.zeros(NumPhysElems_new,float)
        onethird = 1.0/3.0
        for i in range(0,NumPhysElems_new):
            xmid[i] = onethird*(x_new[tnode_new[i,0]]+x_new[tnode_new[i,1]]+x_new[tnode_new[i,2]])
            ymid[i] = onethird*(y_new[tnode_new[i,0]]+y_new[tnode_new[i,1]]+y_new[tnode_new[i,2]])

        # Sample Legendre polynomial on the midpoint of each element
        Mon2Leg = np.zeros((kmax,kmax),float)
        GetMonomialToLegendre(kmax,Mon2Leg)
        MonVals = np.zeros((kmax,p2),float)
        LegVals = np.zeros((kmax,p2),float)
        GetUnstLegendre(meth1,kmax,points_per_dir,zx,zy,Mon2Leg,MonVals,LegVals)
  
        print " Finished creating mesh. "
        print ""

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
                # Solution -- q
                # solution should be found in file
                #     outputdir/q[n1].dat

                qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
                qfile_tmp = "q" + qfile_tmp_tmp[1:]
                qfile = "".join(("".join((outputdir,"/")),qfile_tmp))

                mtmp = NumElems*meqn*kmax
                qtmp = np.zeros(mtmp,float)   
                time = read_qfile(mtmp,qfile,qtmp)                
                qcoeffs_tmp = np.reshape(qtmp,(NumElems,meqn,kmax),'fortran')
                qcoeffs = qcoeffs_tmp[0:NumPhysElems,:,:]

                # Solution -- q
                qsoln_elem = np.zeros((NumPhysElems_new,meqn),float)
                qsoln      = np.zeros((NumPhysNodes_new,meqn),float)
                sample_state2_unst(NumPhysElems,p2,meqn,kmax,LegVals,qcoeffs,qsoln_elem)
                set_soln_at_node_values(NumPhysElems_new,NumPhysNodes_new,
                                        meqn,tcounter_new,tnode_new,qsoln_elem,qsoln)

                # USER SUPPLIED FUNCTION: Plotting function
                plotq2_unst(m-1,meqn,NumPhysElems_new,NumPhysNodes_new,
                            xlow,xhigh,ylow,yhigh,time,x_new,y_new,tnode_new,
                            qsoln,xmid,ymid,qsoln_elem)
            
        plt.ioff()
        print ""
        
    else:
        
        print ""
        print " Error in plotdog2.py: GridType = ",GridType," is not supported."
        print ""
        return -1
            
#----------------------------------------------------------

    

#----------------------------------------------------------
if __name__=='__main__':
    """
    If executed at command line prompt, simply run plotdog2
    """

    import optparse
    parser = optparse.OptionParser(
        usage=''' %%prog (-h | [-p POINTS_PER_DIR] [-o OUTPUT_DIRECTORY] [-t POINT_TYPE] )
    
%s''' % plotdog2.__doc__)

    parser.add_option('-p', '--points-per-dir', type='int', default=1, 
                       help='''number of points per cell to be plotted.
                       Default = 1''')
    parser.add_option('-o', '--output-directory', type='string', default='output', 
                       help='''OUTPUT_DIR = output directory where coefficients
                       can be located''')
    parser.add_option('-t', '--point-type', type='int', default=1,
                       help='''POINT_TYPE = type of points to be plotted''')

    opts, args = parser.parse_args()
    plotdog2(opts.points_per_dir, opts.output_directory, opts.point_type )

#    This is how you could check if the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')

### The old way of calling plotdog2: ###
#   import sys
#   args = sys.argv[1:]   # any command line arguments
#   plotdog2(*args)

