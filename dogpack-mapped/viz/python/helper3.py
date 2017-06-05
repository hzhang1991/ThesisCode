#----------------------------------------------------------
def get_grid_type(outputdir):

    import string

    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    ndims = int(linelist[0])

    if ndims!=3:
        print ""
        print " Incorrect dimension, ndims must be 3. ndims = ",ndims
        print ""
        return -1

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    GridType = str(linelist[0])

    Rparams.close()

    return GridType
#----------------------------------------------------------
    

#----------------------------------------------------------
def read_params(outputdir,params):

    import string

    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    ndims = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    GridType = str(linelist[0])

    if (GridType=="Cartesian"):
        for k in range (0,14):
            linestring = Rparams.readline()
            linelist = string.split(linestring)
            params[k] = float(linelist[0])
    elif (GridType=="Unstructured"):
        print " In 3D currently only GridType=Cartesian is supported. GridType = ",GridType
        print " "
        return -1

    Rparams.close()

    return ndims
#----------------------------------------------------------


#----------------------------------------------------------
def get_kmax(meth1,ndims):

    if (ndims==1):
        return meth1
    elif (ndims==2):
        return int((meth1*(meth1+1))/2)
    elif (ndims==3):
        return int((meth1*(meth1+1)*(meth1+2))/6)
    else:
        print ""
        print " Incorrect dimension in get_kmax, ndims must be 1, 2, or 3. ndims = ",ndims
        print ""
        return -1
    
#----------------------------------------------------------


#----------------------------------------------------------
def read_qfile(mtmp,qfile,qtmp):

    import string

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    time = float(linelist[0])
    
    # store all Legendre coefficients in qtmp
    for k in range (0,mtmp):
        linestring = Rqfile.readline()
        linelist = string.split(linestring)
        qtmp[k] = float(linelist[0])        

    # close file
    Rqfile.close()

    # return time
    return time
#----------------------------------------------------------


#----------------------------------------------------------
#  Sample Legendre polynomial on the midpoint of each element
def GetCart2Legendre(meth1,points_per_dir,s2d,LegVals):

    from math import sqrt
  
    sq3 = sqrt(3.0)
    sq5 = sqrt(5.0)
    sq7 = sqrt(7.0)

    for m in range(0,points_per_dir*points_per_dir):
      xi  = s2d[m,0]
      eta = s2d[m,1]

      xi2 = xi*xi
      xi3 = xi2*xi
      xi4 = xi3*xi

      eta2 = eta*eta
      eta3 = eta2*eta
      eta4 = eta3*eta

      if (meth1==5):
        LegVals[0,m]  = 1.0      
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta      
        LegVals[3,m]  = 3.0*xi*eta 
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[9,m]  = sq7*(2.5*eta3 - 1.5*eta)
        LegVals[10,m] = sq3*sq7*(2.5*xi3 - 1.5*xi)*eta
        LegVals[11,m] = sq3*sq7*(2.5*eta3 - 1.5*eta)*xi
        LegVals[12,m] = 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0)      
        LegVals[13,m] = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0      
        LegVals[14,m] = 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0

      elif (meth1==4):
        LegVals[0,m]  = 1.0      
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta      
        LegVals[3,m]  = 3.0*xi*eta 
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[9,m]  = sq7*(2.5*eta3 - 1.5*eta)

      elif (meth1==3):
        LegVals[0,m]  = 1.0
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta
        LegVals[3,m]  = 3.0*xi*eta
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)

      elif (meth1==2):
        LegVals[0,m]  = 1.0  
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta

      elif (meth1==1):
        LegVals[0,m]  = 1.0
#----------------------------------------------------------


#----------------------------------------------------------
# Turn coefficients into point values on 2D Cartesian grid
#
def sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln):

    index = 0

    for j in range(1,my_old+1):
        for m1 in range(1,points_per_dir+1):
            for i in range(1,mx_old+1):
                for m2 in range(1,points_per_dir+1):
    
                    m = m2 + points_per_dir*(m1-1)
                    index = index + 1
    
                    for n in range(1,meqn+1):
                        qsoln[index-1,n-1] = 0.0

                        for k in range(1,kmax+1):
                            qsoln[index-1,n-1] = qsoln[index-1,n-1] + qcoeffs[k-1,n-1,j-1,i-1]*LegVals[k-1,m-1]
   
#----------------------------------------------------------

#----------------------------------------------------------
# Find out how many slices need to be plotted
#
def get_numslices(outputdir):

    import string

    Fparams = "".join((outputdir,"/qhelp_slice.dat"))
    Sparams = open(Fparams,'r')

    linestring = Sparams.readline()
    linelist = string.split(linestring)
    numslices = int(linelist[0])

    return numslices

#----------------------------------------------------------

#----------------------------------------------------------
# Read-in slice information
#
def read_slice_info(outputdir,slicekind,sliceidx):

    import string

    Fparams = "".join((outputdir,"/qhelp_slice.dat"))
    Sparams = open(Fparams,'r')

    linestring = Sparams.readline()
    linelist = string.split(linestring)
    numslices = int(linelist[0])

    for ns in range(0,numslices):
        linestring = Sparams.readline()
        linelist = string.split(linestring)
        slicekind[ns] = int(linelist[0])

        linestring = Sparams.readline()
        linelist = string.split(linestring)
        sliceidx[ns]  = int(linelist[0])

#----------------------------------------------------------
