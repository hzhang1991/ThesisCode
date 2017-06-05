
#----------------------------------------------------------
def plotq2_unst(outputdir,nframe,
                m,meqn,NumPhysElems,NumPhysNodes,
                xlow,xhigh,ylow,yhigh,time,
                x,y,tnode,qsoln,
                xmid,ymid,qsoln_elem):
    
    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim([xlow,xhigh])
    plt.gca().set_ylim([ylow,yhigh])
    print min(qsoln[:,0])
    p1=plt.tripcolor(x, y, tnode, qsoln[:,m], shading='flat')
    plt.colorbar()
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
