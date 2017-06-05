
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
    from matplotlib.colors import ListedColormap


    def weird_colormap():
       from numpy import linspace, fliplr, flipud, empty
       n = 500
       colmap = empty((n+n+1, 3))
       for i in range(n+n+1):
          colmap[i, :] = (1.0 - (1.0 / (2*n) * i)) ** 10
#    colmap1 = empty((n+1, 3))
#    colmap1[:, 0] = linspace(0.0, 1.0, num = n + 1)
#    colmap1[:, 1] = linspace(0.0, 1.0, num = n + 1)
#    colmap1[:, 2] = 1.0
#    colmap3 = fliplr(flipud(colmap1))
#    colmap = empty((n + n + 1, 3))
#    colmap[0:n, :] = colmap1[0:n, :]
#    colmap[n:,  :] = colmap3
       return colmap


    qsoln=qaug;
    # HACK:
    gamma = 1.4
    gp1 = gamma+1.
    gm1 = gamma-1.
    print ('hi')
    # DENSITY
    if 1:
     #cmap = ListedColormap(weird_colormap())
     rho=qaug[:,:,m]
     grad = np.gradient(rho)
     schlieren = np.sqrt(grad[0]**2 + grad[1]**2)
     m=0
     fig=plt.figure(1)#,figsize=(15,5))
     plt.clf()
     plt.gca().set_aspect('equal')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim([0.0, 7.0])
     #colors=np.linspace(0.0662,7.07,20)
     #plt.contour(xl,yl,qsoln[:,:,m],40)#colors)
     plt.pcolor(xl,yl,qsoln[:,:,m])

     #plt.contourf(xl, yl, schlieren)#,100, cmap = cmap)
     xo=np.linspace(0.0,1.0,40)
     yo1=0.0+0*xo
     yo2=6.0+0*xo
     #plt.plot(xo,yo1,xo,yo2)
     #plt.fill_between(xo,yo1,yo2, color='grey',zorder=3)
     tmp1 = "".join(("Density at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     plt.colorbar()
     plt.plot([0.6],[0.3],'*')
     #fig.tight_layout()
     plt.show()
     #plt.savefig('rarefaction_density.pdf')
    """
    print "HERE!"
    # PRESSURE
    if 1:
     m=1
     press = (gamma-1)*(qsoln[:,:,4]-0.5/qsoln[:,:,0]*(
        qsoln[:,:,1]**2 + qsoln[:,:,2]**2 + qsoln[:,:,3]**2))
     print press[press != 0].min()
     fig=plt.figure(2)
     plt.clf()
     plt.gca().set_aspect('equal')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim( [0.0, 0.2] )
     #colors=np.linspace(0.091,37,38)
     #plt.contour(xl,yl,press,20)#colors)
     plt.pcolor(xl,yl,press)
     xo=np.linspace(0.0,1.0,40)
     yo1=0.0+0*xo
     yo2=6.0+0*xo
     #plt.plot(xo,yo1,xo,yo2)
     #plt.fill_between(xo,yo1,yo2, color='grey',zorder=3)
     tmp1 = "".join(("Pressure at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     plt.colorbar()
     #fig.tight_layout()
     plt.show()
     #plt.savefig('rarefaction_pressure.pdf')

 
    print "HERE2!"
    # Velocity
    if 1:
     m=2
     u = np.zeros(mx,float)
     u = qsoln[:,:,1] #/ qsoln[:,:,0]
        
     fig=plt.figure(3)
     plt.clf()
     plt.gca().set_aspect('equal')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim([-1.0, 1.0])
     levels=np.linspace(0.1,4.2,41)
     plt.contour(xl,yl,rho,levels,colors='k')
     #plt.pcolor(xl,yl,u)
     tmp1 = "".join(("Velocity at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     plt.colorbar()
     #fig.tight_layout()
     plt.show()
     #plt.savefig('rarefaction_velocity.pdf')
    #plt.draw()
    """
#----------------------------------------------------------
