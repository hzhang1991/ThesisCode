
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

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim(xl[0,0],xl[mx,0])
    plt.gca().set_ylim(yl[0,0],yl[0,my])
    plt.pcolor(xl,yl,qaug[:,:,m])
    plt.colorbar()
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
