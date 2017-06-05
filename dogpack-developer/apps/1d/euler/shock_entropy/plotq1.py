from __future__ import print_function  # used for printing to file

#-----------------------------------------------------------------------------#
def plotq1(m,meth1,meqn,mx,time,xc,qsoln,auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # HACK:
    gamma = 1.4
    gp1 = gamma+1.
    gm1 = gamma-1.

    # DENSITY
    m=0
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.5,4.5])
    plt.plot(xc,qsoln[:,m],'bo')
    tmp1 = "".join(("Density at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    # PRESSURE
    m=1

    press = (gamma-1)*(qsoln[:,4]-0.5*(
        qsoln[:,1]**2 + qsoln[:,2]**2 + qsoln[:,3]**2)/qsoln[:,0])
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0, 12])
    plt.plot(xc,press,'bo')
    tmp1 = "".join(("Pressure at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
 
    # Velocity
    m=2
    u = np.zeros(mx,float)
    u = qsoln[:,1] / qsoln[:,0]
        
    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.5,3])
    plt.plot(xc,u,'bo')
    tmp1 = "".join(("Velocity at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    data = np.array( [xc, qsoln[:,0], qsoln[:,1], qsoln[:,4] ] )
#   fmt = '%.15e'
#   with open( 'dg_picture.dat', 'wb' ) as f:
#     print( fmt % time, file=f )         # time instant on first row
#     np.savetxt( 'dg_picture.dat', np.transpose( data ), fmt=fmt )
    np.savetxt( 'dg_picture.dat', np.transpose( data ) )
 
#-----------------------------------------------------------------------------#
