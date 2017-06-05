
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
    qsoln=qaug;
    # HACK:
    gamma = 1.4
    gp1 = gamma+1.
    gm1 = gamma-1.
    print ('hi')
    # DENSITY
    if 1:
     m=0
     fig=plt.figure(1)
     plt.clf()
     plt.gca().set_aspect('auto')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim([0.0, 7.0])
     colors=np.linspace(0.0662,7.07,20)
     plt.contour(xl,yl,qsoln[:,:,m],colors)
     xo=np.linspace(0.0,1.0,40)
     yo1=0.0+0*xo
     yo2=6.0+0*xo
     plt.plot(xo,yo1,xo,yo2)
     plt.fill_between(xo,yo1,yo2, color='grey',zorder=3)
     tmp1 = "".join(("Density at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     fig.tight_layout()
     plt.savefig('rarefaction_density.pdf')


    # PRESSURE
    if 1:
     m=1
     press = (gamma-1)*(qsoln[:,:,4]-0.5/qsoln[:,:,0]*(
        qsoln[:,:,1]**2 + qsoln[:,:,2]**2 + qsoln[:,:,3]**2))
     print press[press != 0].min()
     fig=plt.figure(2)
     plt.clf()
     plt.gca().set_aspect('auto')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim( [0.0, 0.2] )
     colors=np.linspace(0.091,37,38)
     plt.contour(xl,yl,press,colors)
     xo=np.linspace(0.0,1.0,40)
     yo1=0.0+0*xo
     yo2=6.0+0*xo
     plt.plot(xo,yo1,xo,yo2)
     plt.fill_between(xo,yo1,yo2, color='grey',zorder=3)
     tmp1 = "".join(("Pressure at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     fig.tight_layout()
     plt.savefig('rarefaction_pressure.pdf')

 
    # Velocity
    if 1:
     m=2
     u = np.zeros(mx,float)
     u = qsoln[:,:,1] #/ qsoln[:,:,0]
        
     fig=plt.figure(3)
     plt.clf()
     plt.gca().set_aspect('auto')
     #plt.gca().set_xlim([xc[0],xc[mx-1]])
     #plt.gca().set_ylim([-1.0, 1.0])
     colors=np.linspace(0.091,37,40)
     plt.contourf(xl,yl,u, 20)
     tmp1 = "".join(("Velocity at t = ",str(time)))
     title = "".join((tmp1,"     [DoGPack]"))
     plt.title(title)
     fig.tight_layout()
     #plt.savefig('rarefaction_velocity.pdf')
    
#----------------------------------------------------------
