def plotdog2np(points_per_dir_in="1",outputdir="output",point_type_in="1",component_of_solution="1"):
    """
Generic code for plotting DoGPack output in matplotlib.

Execute via

    $ python $DOGPACK/viz/python/plotdog2np.py

from an application directory.   For help type

    $ python $DOGPACK/viz/python/plotdog2np.py -h
 
to see a list of options.

    Usage: plotdog2np( points_per_dir_in="1", outputdir="output", point_type_in="1", component_of_solution"1")
 
   points_per_dir = points per direction (spatial dimension)
 
   outputdir = location of output directory

   point_type = 1:   uniform points on each element
              = 2:   Gauss-Legendre points on each element

   component_of_solution = which component of the solution q to plot
"""



    import os
    import sys
    from math import sqrt
    import numpy as np
    import matplotlib.pyplot as plt
    from helper2 import get_grid_type
    from helper2 import read_params
    from helper2 import get_kmax
    from helper2 import GetCart2Legendre
    from helper2 import read_qfile
    from helper2 import sample_state2_cart_mod
    from helper2 import read_mesh_params
    from helper2 import read_node
    from helper2 import read_tnode
    from helper2 import GetMonomialToLegendre
    from helper2 import GetUnstLegendre
    from helper2 import set_tcounter
    from helper2 import set_soln_at_node_values
    from helper2 import sample_state2_unst
    from helper2 import DivideUnstMesh
    
    points_per_dir = int(points_per_dir_in)
    point_type = int(point_type_in)
    component_of_solution = int(component_of_solution)

    TF = os.path.exists(outputdir)
    if TF==False:
        print ""
        print "    Directory not found, outputdir =",outputdir
        print ""
        return -1

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
        mx      = int(params[4])
        my      = int(params[5])
        xlow    = params[6]
        xhigh   = params[7]
        ylow    = params[8]
        yhigh   = params[9]
        datafmt = int(params[10])
    elif (GridType=="Unstructured"):
        params = np.zeros(4, int)
        read_params(outputdir,params)
        meqn    = params[0]
        maux    = params[1]
        nplot   = params[2]
        meth1   = params[3]

    kmax = get_kmax(meth1,2)
    if (component_of_solution<1):
        print ""
        print "  Error, need component_of_solution > 1,  component_of_solution = ",component_of_solution
        print ""
        return -1
    elif (component_of_solution>meqn):
        print ""
        print "  Error, need component_of_solution <=",meqn,",  component_of_solution = ",component_of_solution
        print ""
        return -1
        
    print ""
    print "              GridType = ",GridType
    print "        points_per_dir = ",points_per_dir
    print "            point_type = ",point_type
    print "             outputdir = ",outputdir
    print " component_of_solution = ",component_of_solution
    print ""

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)

    if (GridType=="Cartesian"):
        plotq2_file  = os.path.abspath("plotq2np_cart.py")
        local_plotq2 = os.path.exists(plotq2_file)
        if local_plotq2==False:
            from plotq2np_cart_default import plotq2np_cart
        else:
            from plotq2np_cart import plotq2np_cart

        # Grid information
        mx_old = mx
        my_old = my
        dx_old = (xhigh-xlow)/mx_old
        dy_old = (yhigh-ylow)/my_old
        mx = mx*points_per_dir
        my = my*points_per_dir
        dx = (xhigh-xlow)/mx
        dy = (yhigh-ylow)/my

        if (point_type==1):
            xl = np.zeros((mx+1,my+1),float)
            yl = np.zeros((mx+1,my+1),float)
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
  
        m = component_of_solution;
        n1 = 0;

        while (n1<=nplot):
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
            sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln)
            qsoln = np.reshape(qsoln,(mx,my,meqn),'F')

            # USER SUPPLIED FUNCTION
            plotq2np_cart(m-1,meth1,meqn,mx,my,time,xc,yc,xl,yl,qsoln,n1)
            n1=n1+1
          
        print ""

    elif (GridType=="Unstructured"):
        plotq2_file  = os.path.abspath("plotq2np_unst.py")
        local_plotq2 = os.path.exists(plotq2_file)
        if local_plotq2==False:
            from plotq2np_unst_default import plotq2np_unst
        else:
            from plotq2np_unst import plotq2np_unst

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

        m = component_of_solution
        n1 = 0;

        while (n1<=nplot):        
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
            plotq2np_unst(m-1,meqn,NumPhysElems_new,NumPhysNodes_new,
                          xlow,xhigh,ylow,yhigh,time,x_new,y_new,tnode_new,
                          qsoln,xmid,ymid,qsoln_elem,n1)
            n1 = n1+1
            
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
    If executed at command line prompt, simply run plotdog2np
    """

    # TODO - swap out this parser for argparse.  See
    # http://docs.python.org/dev/library/argparse.html
    # for full details.  (optparse is now deprecated).
    import optparse
    parser = optparse.OptionParser(
        usage=''' %%prog (-h | [-p POINTS_PER_DIR] [-o OUTPUT_DIRECTORY] [-t POINT_TYPE] [-m COMPONENT_OF_SOLUTION])
    
%s''' % plotdog2np.__doc__)

    parser.add_option('-p', '--points-per-dir', type='int', default=1, 
                       help='''number of points per cell to be plotted.
                       Default = 1''')
    parser.add_option('-o', '--output-directory', type='string', default='output', 
                       help='''OUTPUT_DIR = output directory where coefficients
                       can be located''')
    parser.add_option('-t', '--point-type', type='int', default=1,
                       help='''POINT_TYPE = type of points to be plotted''')
    parser.add_option('-n', '--component', type='int', default=1,
                       help='''COMPONENT_OF_SOLUTION = component of solution to be plotted''')

    opts, args = parser.parse_args()

    # TODO - add an option for which component to plot
    plotdog2np(opts.points_per_dir, opts.output_directory, opts.point_type, 1)



