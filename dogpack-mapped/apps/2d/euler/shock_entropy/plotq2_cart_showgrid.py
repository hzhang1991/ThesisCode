
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
                xl,yl,
                qaug):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import mapc2p as mp

    data=np.loadtxt('mapped.txt')
    
    xp=data[:,0]
    yp=data[:,1]
    xpl,ypl=mp.mapc2p(xl,yl)
    n4=xpl.shape[0]
    m4=xpl.shape[1]
    print n4,m4,xp.shape

    #xp=np.reshape(xp,(n4,m4))
    #yp=np.reshape(yp,(n4,m4))

    xp=np.array(xpl)
    yp=np.array(ypl)

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    #plt.gca().set_xlim(xl[0,0],xl[mx,0])
    #plt.gca().set_ylim(yl[0,0],yl[0,my])
    #im=plt.pcolor(xp,yp,qaug[:,:,m], cmap='OrRd')
    for  i1 in range(mx):
     for  j1 in range(my):
       #print [xpl[i1,j1],xpl[i1+1,j1]]
       print i1,j1
       if xp[i1,j1]!=xp[i1+1,j1]:
        plt.plot([xp[i1,j1], xp[i1+1,j1]],[yp[i1,j1], yp[i1+1,j1]],'k')
       if xp[i1,j1]!=xp[i1,j1+1]:
        plt.plot([xp[i1,j1], xp[i1,j1+1]],[yp[i1,j1], yp[i1,j1+1]],'k')
       if yp[i1,j1]!=yp[i1+1,j1]:
        plt.plot([xp[i1,j1], xp[i1+1,j1]],[yp[i1,j1], yp[i1+1,j1]],'k')
       if yp[i1,j1]!=yp[i1,j1+1]:
        plt.plot([xp[i1,j1], xp[i1,j1+1]],[yp[i1,j1], yp[i1,j1+1]],'k')
    for  i1 in range(mx):
       if xp[i1,my]!=xp[i1+1,my]:
        plt.plot([xp[i1,my], xp[i1+1,my]],[yp[i1,my], yp[i1+1,my]],'k')

    for  j1 in range(my):
       if yp[mx,j1]!=yp[mx,j1+1]:
        plt.plot([xp[mx,j1], xp[mx,j1+1]],[yp[mx,j1], yp[mx,j1+1]],'k')

    #plt.colorbar(im, use_gridspec=True)
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,""))
    plt.draw()
    plt.tight_layout()
    plt.savefig('plots/figure_%s.pdf'%nframe,transparent=True) 
    
#----------------------------------------------------------
