
#----------------------------------------------------------
def plotq2_unst(m,meqn,NumPhysElems,NumPhysNodes,
                xlow,xhigh,ylow,yhigh,time,x,y,
                tnode,qsoln,xmid,ymid,qsoln_elem):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import sqrt
    from math import pow
    from math import cos
    from math import sin
    from math import pi

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim([xlow,xhigh])
    plt.gca().set_ylim([ylow,yhigh])
    p1=plt.tripcolor(x, y, tnode, qsoln[:,m], shading='faceted', vmin=0.0, vmax=1.0)
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.colorbar()
    plt.draw()


    x0 = -0.25*cos(2.0*pi*time) + 0.50
    y0 =  0.25*sin(2.0*pi*time) + 0.50
    r = np.zeros(NumPhysElems,float)
    for i in range(0,NumPhysElems):
        r[i] = sqrt(pow(xmid[i]-x0,2)+pow(ymid[i]-y0,2))
    ind = r.argsort()
    qscat_ex = np.zeros((NumPhysElems,meqn),float)
    qex(NumPhysElems,meqn,r,qscat_ex)

    err = np.linalg.norm(qscat_ex[:,m]-qsoln_elem[:,m])/np.linalg.norm(qscat_ex[:,m])
    print ""
    print "  Error = ",'{:e}'.format(err)
    print ""
    
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([0.0,0.5])
    plt.gca().set_ylim([0.0,1.0])
    plt.plot(r[ind],qscat_ex[ind,m],'k-')
    plt.plot(r[ind],qsoln_elem[ind,m],'bo')
    tmp1 = "".join(("Scattor plot of q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------


#----------------------------------------------------------
def qex(NumPhysElems,meqn,r,qscat_ex):

    from math import pow
    from math import cos
    from math import pi
    
    for i in range(NumPhysElems):
        if (r[i]<0.2):
            qscat_ex[i,0] = pow( cos(5.0/2.0*pi*r[i]) ,6)

#----------------------------------------------------------
