
#----------------------------------------------------------
def plotq2np_cart(m,meth1,meqn,mx,my,time,xc,yc,xl,yl,qsoln,FRAME):
    
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(FRAME)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim(xl[0,0],xl[mx,0])
    plt.gca().set_ylim(yl[0,0],yl[0,my])
    plt.pcolor(xl,yl,qsoln[:,:,m])
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)

    stmp = list(str(1000+FRAME));
    stmp[0]="0";
    fname = "".join(("figure","".join(stmp),".jpg"));
    plt.savefig(fname);

    print " Finished creating file for FRAME = ",FRAME,"  filename = ",fname
    
#----------------------------------------------------------
