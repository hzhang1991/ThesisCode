

#----------------------------------------------------------
def read_params(outputdir,params,qhelpname):

    import string
    import numpy as np

    Fparams = "".join(("".join((outputdir,"/")),qhelpname));

    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    ndims = int(linelist[0])

    if ndims!=1:
        print ""
        print " Incorrect dimension, ndims must be 1. ndims = ",ndims
        print ""
        return -1

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    GridType = str(linelist[0])

    for k in range (0,8):
        linestring = Rparams.readline()
        linelist = string.split(linestring)
        params[k] = np.float64(linelist[0])

    Rparams.close()

    return GridType
#----------------------------------------------------------


#----------------------------------------------------------
# Sample basis functions on mesh
# size of phi = (numpts,meth1)
def SampleBasis1(numpts,meth1,phi):

    from math import sqrt
    import numpy as np
    
    sq3 = sqrt(3.0)
    sq5 = sqrt(5.0)
    sq7 = sqrt(7.0)
    sq11 = sqrt(11.0)

    for n1 in range(0,numpts):
        xi = -1.0 + (2.0*n1+1.0)/np.float64(numpts)
        xi2 = xi*xi
        xi3 = xi2*xi
        xi4 = xi3*xi
        xi5 = xi4*xi

        phi[n1,0] = 1
        
        if meth1>1:
            phi[n1,1] = sq3*xi

        if meth1>2:
            phi[n1,2] = sq5*(1.5*xi2 - 0.5)

        if meth1>3:
            phi[n1,3] = sq7*(2.5*xi3 - 1.5*xi)

        if meth1>4:
            phi[n1,4] = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0

        if meth1>5:
            phi[n1,5] = (63.0/8.0)*sq11 * ( xi5 - (10.0/9.0)*xi3 + (5.0/21.0)*xi )
            
#----------------------------------------------------------


#----------------------------------------------------------
def read_qfile(mtmp,qfile,qtmp):

    import string
    import numpy as np

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    time = np.float64(linelist[0])
    
    # store all Legendre coefficients in qtmp
    for k in range (0,mtmp):
        linestring = Rqfile.readline()
        linelist = string.split(linestring)
        qtmp[k] = np.float64(linelist[0])

    # close file
    Rqfile.close()

    # return time
    return time
#----------------------------------------------------------
