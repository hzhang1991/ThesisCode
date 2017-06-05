

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
def plotq2_cart(outputdir,nframe,
                m,meth1,meqn,time,
                points_per_dir,LegVals,
                xlow,xhigh,ylow,yhigh,
                mx,my,
                dx,dy,
                mx_old,my_old,
                dx_old,dy_old,
                xc,yc,
                xl,yl,qsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import mapc2p as mp

    plt.figure(1)
    plt.clf()
    #xp=xl
    #yp=yl
    #xpl,ypl=mp.mapc2p(xl,yl)
    n4=xl.shape[0]
    m4=xl.shape[1]
    data=np.loadtxt('mapped.txt')
    
    xp=data[:,0]
    yp=data[:,1]
    #print xp.shape,yp.shape,n4,m4
    #xp=np.reshape(xp,(n4,m4))
    #yp=np.reshape(yp,(n4,m4))
    #plt.gca().set_aspect('equal')
    #plt.gca().set_xlim(xl[0,0],xl[mx,0])
    #plt.gca().set_ylim(yl[0,0],yl[0,my])
    #plt.contourf(xl,yl,qsoln[:,:,m])
    plt.pcolor(xl,yl,qsoln[:,:,m])
    #plt.pcolor(xp,yp,qsoln[:,:,m])

    print xl


    """
    plt.plot(xp,yp,'ko')
    for  i1 in range(mx):
     for  j1 in range(my):
       #print [xpl[i1,j1],xpl[i1+1,j1]]
       if yp[i1,j1]!=yp[i1+1,j1]:
        plt.plot([xp[i1,j1], xp[i1+1,j1]],[yp[i1,j1], yp[i1+1,j1]],'k')
       if yp[i1,j1]!=yp[i1,j1+1]:
        plt.plot([xp[i1,j1], xp[i1,j1+1]],[yp[i1,j1], yp[i1,j1+1]],'k')
    """
    plt.colorbar();
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
