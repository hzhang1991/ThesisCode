#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np
import matplotlib.pyplot as plt
import mapc2p as mp

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
# Turn coefficients into point values on 2D Cartesian grid
#

"""
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
                            qsoln[index-1,n-1] = qsoln[index-1,n-1] + qcoeffs[k-1,n-1,j-1,i-1]*LegVals[i-1,j-1,0,k-1]
"""
#----------------------------------------------------------


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
        LegVals[0,m]  = 1.0/4.0*(1.0-xi)*(1.0-eta) 
        LegVals[1,m]  = 1.0/4.0*(1.0+xi)*(1.0-eta)
        LegVals[2,m]  = 1.0/4.0*(1.0+xi)*(1.0+eta)
        LegVals[3,m]  = 1.0/4.0*(1.0-xi)*(1.0+eta)

      elif (meth1==1):
        LegVals[0,m]  = 1.0
#----------------------------------------------------------

def read_params(outputdir,params):

    import string

    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    linestring = Rparams.readline()
    GridType = 'Cartesian'
    if (GridType=="Cartesian"):
        for k in range (0,11):
            linestring = Rparams.readline()
            linelist = string.split(linestring)
            params[k] = float(linelist[0])
    elif (GridType=="Unstructured"):
        for k in range (0,4):
            linestring = Rparams.readline()
            linelist = string.split(linestring)
            params[k] = int(linelist[0])

    Rparams.close()

    return GridType

def finc(xx,yy,theta,c,t,Iv,x0,y0):
 w=5.;
 xi=Iv[0]*(xx-x0)+Iv[1]*(yy-y0)-c*t;
 M=((xi<0.0)*(xi>-1.0/w)*(xx<0.5));
 P=np.zeros(xx.shape);
 P[M]=1.0/2.*(1.-np.cos(2.*np.pi*w*xi[M]))**4;
 return P

def fref(xx,yy,theta,c,t,Iv,A,x0,y0):
 w=5.;
 xi=Iv[0]*(xx-x0)+Iv[1]*(yy-y0)-c*t;
 M=((xi<0.0)*(xi>-1.0/w)*(xx<0.5));
 P=np.zeros(xx.shape);
 P[M]=A*1.0/2.*(1.-np.cos(2.*np.pi*w*xi[M]))**4;
 return P

def ftran(xx,yy,theta,c,t,Iv,A,x0,y0):
 w=5.;
 xi=Iv[0]*(xx-x0)+Iv[1]*(yy-y0)-c*t;
 M=((xi<0.0)*(xi>-0.5/w)*(xx>=0.5));
 P=np.zeros(xx.shape);
 P[M]=A*1.0/2.*(1.-np.cos(4.*np.pi*w*xi[M]))**4;
 return P


def interface(xx,yy): 
 P0=np.zeros(xx.shape);
 H=(xx>0.0)*(xx<.2)

 P0[H]=(1.0/2.0*(1.0+np.cos(np.pi*(xx[H]-0.1)*10.0)))**5;
 z=xx[P0==max(P0)]
 #P0=(1.0/2.0*(1.0+np.cos(np.pi*(xx-0.35)*4.0)))**3
 #print xx[P0==max(P0)],P0[P0==max(P0)],(1.0/2.0*(1.0+np.cos(np.pi*(z-0.35)*4.0)))**3;
 return P0

def main( ):
    '''Write some help documentation here
'''
    err1=[];nnn1=[];
    err2=[];nnn2=[];
    err3=[];nnn3=[];

    my_dictionary = {}
    old_err = l = 0
    old_err=1.0
    new_err1=1.0
    new_err2=1.0
    new_err3=1.0
    old_n =1.0
    new_n=1.0
      ###FRAME TO TEST
    n1=10
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  l

        folder = 'output_00%s' % l
        if( not os.path.exists(folder) ):
            break
        my_dictionary['curr_folder'] = folder
        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = l

        try:
            qpre = np.loadtxt(folder + "/q0010.dat")[1:]
        except IOError:
            print "Could not find data in folder %s" % folder
            return
        print directory_num,folder
###############################################
        params = np.zeros(11, float)
        read_params(folder,params)
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
        # Grid information
        mx_old = mx
        my_old = my
        dx_old = (xhigh-xlow)/mx_old
        dy_old = (yhigh-ylow)/my_old
        point_type=2
        points_per_dir=5;
        mx = mx*points_per_dir
        my = my*points_per_dir
        dx = (xhigh-xlow)/mx
        dy = (yhigh-ylow)/my

        kmax=int((meth1*(meth1+1)/2))
        phil=np.zeros((mx,my,points_per_dir,kmax),float)
        #print 'mx=',mx,'my=',my,'xlow=',xlow,'xhigh=',xhigh,'dx=',dx,params

        xc = np.zeros((mx,my),float)
        yc = np.zeros((mx,my),float)

        xol = np.zeros((mx,my),float)
        yol = np.zeros((mx,my),float)

        xor = np.zeros((mx,my),float)
        yor = np.zeros((mx,my),float)

        sxc = np.zeros((mx,my),float)
        syc = np.zeros((mx,my),float)
        
        if (point_type==1):
            xl = np.zeros((mx+1,my+1),float)
            yl = np.zeros((mx+1,my+1),float)


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

            sq3 = np.sqrt(3.0)
            sq5 = np.sqrt(5.0)
            sq7 = np.sqrt(7.0)

            # 1D quadrature points
            s1d = np.zeros(points_per_dir,float)
            if (points_per_dir==1):
                s1d[0] = -1.0
            elif (points_per_dir==2):
                s1d[0] = -1.0/sq3
                s1d[1] =  1.0/sq3
            elif (points_per_dir==3):
                s1d[0] = -sq3/sq5
                s1d[1] =  0.0
                s1d[2] =  sq3/sq5
            elif (points_per_dir==4):
                s1d[0] = -np.sqrt(3.0+np.sqrt(4.8))/sq7
                s1d[1] = -np.sqrt(3.0-np.sqrt(4.8))/sq7
                s1d[2] =  np.sqrt(3.0-np.sqrt(4.8))/sq7
                s1d[3] =  np.sqrt(3.0+np.sqrt(4.8))/sq7
            elif (points_per_dir==5):
                s1d[0] = -np.sqrt(5.0 + np.sqrt(40.0/7.0))/3.0
                s1d[1] = -np.sqrt(5.0 - np.sqrt(40.0/7.0))/3.0
                s1d[2] =  0.0
                s1d[3] =  np.sqrt(5.0 - np.sqrt(40.0/7.0))/3.0
                s1d[4] =  np.sqrt(5.0 + np.sqrt(40.0/7.0))/3.0

            kk=-1;
            s2d = np.zeros((points_per_dir*points_per_dir,2),float)
            for jj in range(0,points_per_dir):
                for ii in range(0,points_per_dir):
                    kk = kk+1
                    s2d[kk,0] = s1d[ii]
                    s2d[kk,1] = s1d[jj]

            xx = np.zeros(mx,float)
            yy = np.zeros(my,float)

            sx = np.zeros(mx,float)
            sy = np.zeros(my,float)

            xor1 = np.zeros(mx+points_per_dir,float)
            yor1 = np.zeros(my+points_per_dir,float)

            kk=0
            xtmp = xlow-0.5*dx_old
            xtmp1= xlow-dx_old
            for i in range(0,mx_old):
                xtmp = xtmp + dx_old
                xtmp1 = xtmp1 + dx_old
                for m in range(0,points_per_dir):
                    xx[kk+m] = xtmp + 0.5*dx_old*s1d[m]
                    sx[kk+m]=s1d[m]
                    xor1[kk+m]=xtmp1
                kk = kk + points_per_dir
            xor1[mx:mx+points_per_dir]=xtmp1+dx_old

            kk=0
            ytmp = ylow-0.5*dy_old
            ytmp1= ylow-dy_old
            for j in range(0,my_old):
                ytmp = ytmp + dy_old
                ytmp1 = ytmp1 + dy_old
                for m in range(0,points_per_dir):
                    yy[kk+m] = ytmp + 0.5*dy_old*s1d[m]
                    sy[kk+m] = s1d[m]
                    yor1[kk+m]=ytmp1
                kk = kk + points_per_dir
            yor1[my:my+points_per_dir]=ytmp1+dy_old



            for i in range(0,mx):
                for j in range(0,my):
                    xc[i,j] = xx[i]
                    yc[i,j] = yy[j]
                    sxc[i,j] = sx[i]
                    syc[i,j] = sy[j]

            ix=-1
            for i in range(0,mx_old):
              for l1 in range(points_per_dir):
                ix=ix+1
                iy=-1
                for j in range(0,my_old):
                  for l2 in range(points_per_dir):
                    iy=iy+1 
                    xol[ix,iy] = xor1[points_per_dir*(i)]
                    yol[ix,iy] = yor1[points_per_dir*(j)]


            ix=-1
            for i in range(0,mx_old):
              for l1 in range(points_per_dir):
                ix=ix+1
                iy=-1
                for j in range(0,my_old):
                  for l2 in range(points_per_dir):
                    iy=iy+1
                    xor[ix,iy] = xor1[points_per_dir*(i)+points_per_dir]
                    yor[ix,iy] = yor1[points_per_dir*(j)+points_per_dir]

   

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


        qfile_tmp_tmp = "".join((str(0+10000),".dat"))
        qfile_tmp = "q" + qfile_tmp_tmp[1:]
        qfile = "".join(("".join((folder,"/")),qfile_tmp))
        mtmp = mx_old*my_old*meqn*kmax
        qtmp = np.zeros(mtmp,float)   
        time = read_qfile(mtmp,qfile,qtmp)
        qcoeffsold = np.reshape(qtmp,(kmax,meqn,my_old,mx_old))


        qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
        qfile_tmp = "q" + qfile_tmp_tmp[1:]
        qfile = "".join(("".join((folder,"/")),qfile_tmp))

        mtmp = mx_old*my_old*meqn*kmax
        qtmp = np.zeros(mtmp,float)   
        time = read_qfile(mtmp,qfile,qtmp)
        qcoeffs = np.reshape(qtmp,(kmax,meqn,my_old,mx_old))
        #qdiff=abs(qcoeffs-qcoeffsold)        
        qsoln = np.zeros((mx*my,meqn),float)
        xc1=np.reshape(xc,(mx*my))
        yc1=np.reshape(yc,(mx*my))

      
        xleg=np.array(xc1)
        yleg=np.array(yc1)
        
          
        xp1,yp1=mp.mapc2p(xc1,yc1)
        xpl,ypl=mp.mapc2p(xl,yl)

        xpllo,ypllo=mp.mapc2p(xol,yol)
        xplro,yplro=mp.mapc2p(xol,yor)
        xprlo,yprlo=mp.mapc2p(xor,yol)
        xprro,yprro=mp.mapc2p(xor,yor)

        xpy=np.array(xc)
        ypy=np.array(yc)
        xpy1=np.array(xc)
        ypy1=np.array(yc)
        coeff=np.array(xc)
        mcoeff=0.0
        for i in range(0,mx):
             for j in range(0,my):

              xpp1= xpllo[i,j];
              ypp1= ypllo[i,j];
              xpp2= xprlo[i,j];
              ypp2= yprlo[i,j];
              xpp3= xprro[i,j];
              ypp3= yprro[i,j];
              xpp4= xplro[i,j];
              ypp4= yplro[i,j];
              
              xmax=max(max(xpp1,xpp2),max(xpp3,xpp4));xmin=min(min(xpp1,xpp2),min(xpp3,xpp4));
              ymax=max(max(ypp1,ypp2),max(ypp3,ypp4));ymin=min(min(ypp1,ypp2),min(ypp3,ypp4));
              dx1=xmax-xmin;dy1=ymax-ymin;
              xin=xp1[i*my+j]#xmin+dx1/2*(1+1)#xp1[i*my+j]#xmin+dx1/2*(1-0)#xp1[i*my+j];
              yin=yp1[i*my+j]#ymin+dy1/2*(1-0.5)#yp1[i*my+j]#ymin+dy1/2*(1-0)#yp1[i*my+j];

 
              xo,yo=mp.squaremap(sxc[i,j],syc[i,j],xpp1,xpp2,xpp3,xpp4,ypp1,ypp2,ypp3,ypp4)
              #print sxc[i,j],syc[i,j],xo,yo,xpp1,xpp2,xpp3,xpp4,ypp1,ypp2,ypp3,ypp4
              #adsf
              xpy[i,j]=xo
              ypy[i,j]=yo
              xpy1[i,j]=xmin
              ypy1[i,j]=ymin
              #mn=qdiff[:,0,j,i]
              #coeff[i,j]=max(mn)
              #mcoeff=max(mcoeff,coeff[i,j])
              #for l1 in range(0,kmax):
              #  phil[i,j,0,l1]=phi(xin,yin,xmin,ymin,dx1,dy1,l1+1)
 
   
        sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln)
        qsoln = np.reshape(qsoln,(mx,my,meqn),'F')
        qex1=interface(np.reshape(xpy,(mx*my)),np.reshape(ypy,(mx*my)))#(xp1,yp1)
        qex=np.reshape(qex1,(mx,my)) 
        
        plt.figure(1)
        plt.clf()
        xc2=xpy#xpl[0:mx,0:my]#np.reshape(xp1,(mx,my))
        yc2=ypy#ypl[0:mx,0:my]#np.reshape(yp1,(mx,my))
        #plt.gca().set_aspect('equal')
        #plt.gca().set_xlim(xl[0,0],xl[mx,0])
        #plt.gca().set_ylim(yl[0,0],yl[0,my])

       # qex[yc<0.2]=0.0
       # qsoln[:,:,0][yc<0.2]=0.0
        #plt.pcolor(xl,yl,qex)
        #plt.show()
        #plt.pcolor(xpl,ypl,qex)
        #print qcoeffs[0,0,:,6],xpy1[6,:]
        #plt.colorbar();

        #plt.figure(2)
        #plt.clf()
        #plt.pcolor(xpl,ypl,qsoln[:,:,0])
        #plt.colorbar();

        #plt.figure(3)
        #plt.clf()
        #plt.pcolor(xpl,ypl,qex-qsoln[:,:,0])
        #for  i1 in range(mx):
        # for  j1 in range(my):
        #   #print [xpl[i1,j1],xpl[i1+1,j1]]
        #   plt.plot([xpl[i1,j1], xpl[i1+1,j1]],[ypl[i1,j1], ypl[i1+1,j1]],'k')
        #   plt.plot([xpl[i1,j1], xpl[i1,j1+1]],[ypl[i1,j1], ypl[i1,j1+1]],'k')
        #plt.colorbar();
        #plt.show()
        
    
####################################

        diff = np.reshape((qsoln[:,:,0]-qex),mx*my)

        h=dx_old
        old_n=new_n
        new_n=mx_old
        old_err1=new_err1
        old_err2=new_err2
        old_err3=new_err3
        #print "h=",h,mx,my
        new_err1 = max(abs(diff))#h**2*np.linalg.norm(diff,1) #np.linalg.norm(diff,np.inf)#
        new_err2 = h**2*np.linalg.norm(diff,1)
        new_err3 = h*np.linalg.norm(diff,2)


        """
        if( old_err > 0 and new_err > 0 ):
            result = r1 + folder + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), float(new_n)/float(old_n)) } 
            errors['rat'] = log( ( old_err/new_err), float(new_n)/float(old_n) )

        else:
            result = r1 + folder + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 
            errors['rat'] = old_err/new_err
        """
        #print 'ratio=',old_err/new_err,log(old_err/new_err)/log(new_n/old_n)
        rat1=np.log(old_err1/new_err1)/np.log(float(new_n)/float(old_n))
        rat2=np.log(old_err2/new_err2)/np.log(float(new_n)/float(old_n))
        rat3=np.log(old_err3/new_err3)/np.log(float(new_n)/float(old_n))
        #print 'old error',old_err,'old n',old_n,old_err/new_err,float(new_n)/float(old_n)
        print "    Error = %s; Order = %s " %(new_err1,rat1)
        print "    Error = %s; Order = %s " %(new_err2,rat2)
        print "    Error = %s; Order = %s " %(new_err3,rat3)
        err1=err1+[new_err1]
        nnn1=nnn1+[new_n]
        err2=err2+[new_err2]
        nnn2=nnn2+[new_n]
        err3=err3+[new_err3]
        nnn3=nnn3+[new_n]

        l = l + 1
        print "here! ",l
    print nnn1,err1
    plt.loglog(nnn1,err1)
    plt.show()

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
