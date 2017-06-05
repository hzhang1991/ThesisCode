import numpy as np


def mapx(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4):
    x=0.25*((1.0-xi)*(1.0-eta)*xp1+(1.0+xi)*(1.0-eta)*xp2+(1.0-xi)*(1.0+eta)*xp3+(1.0+xi)*(1.0+eta)*xp4);
    return x;

def mapy(xi,eta,xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4):
    y=0.25*((1.0-xi)*(1.0-eta)*yp1+(1.0+xi)*(1.0-eta)*yp2+(1.0-xi)*(1.0+eta)*yp3+(1.0+xi)*(1.0+eta)*yp4);
    return y;


def mapc2p(xc,yc):
        """
    Specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
        """

        import string
        import os 
        from numpy import maximum,minimum,sqrt,ones,sign
        #p = os.system('ls')


	xp = np.array(xc);
	yp = np.array(yc);

	xp2 = np.array(xc);
	yp2 = np.array(yc);

	for i in range(1):
	  r1 = 1.5
	  d = maximum(abs(xp),abs(yp));            # value on diagonal of computational grid
	  d = maximum(d, 1e-14);           # to avoid divide by zero at center
	  d1 = d*r1/sqrt(2.);
	  R = r1*d;
	  #R = r1*ones(d1.shape);
	  #R = r1^2 ./ (r2*d);
	  
          xc1=abs(xc)
          yc1=abs(yc)
          
	   
	  xp2 = d1/d * xc1;
	  yp2 = d1/d * yc1;
	  center = d1 - sqrt(R**2 - d1**2);
	  
	  ij = (yc1 > xc1-1.0e-15);
	  yp2[ij] = center[ij] + sqrt(R[ij]**2 - xp2[ij]**2);

	  ij = (xc1 > yc1-1.0e-15);
	  xp2[ij] = center[ij] + sqrt(R[ij]**2 - yp2[ij]**2);
	  
	  
	  xp = sign(xp) * xp2;
	  yp = sign(yp) * yp2;
	  
        xpm=np.array(xc);ypm=np.array(yc);
  
        return xp,yp

def squaremap(xi,eta,x1,x2,x3,x4,y1,y2,y3,y4):

    xout=(1.0/4.0)*x1+(1.0/4)*x2+(1.0/4.)*x3+(1.0/4.0)*x4+(-.2500000000*x1+.2500000000*x2+.2500000000*x3-.2500000000*x4)*xi+(-.2500000000*x1-.2500000000*x2+.2500000000*x3+.2500000000*x4)*eta+(.2500000000*x1-.2500000000*x2+.2500000000*x3-.2500000000*x4)*xi*eta

    yout=(1.0/4.0)*y1+(1.0/4)*y2+(1.0/4.)*y3+(1.0/4.0)*y4+(-.2500000000*y1+.2500000000*y2+.2500000000*y3-.2500000000*y4)*xi+(-.2500000000*y1-.2500000000*y2+.2500000000*y3+.2500000000*y4)*eta+(.2500000000*y1-.2500000000*y2+.2500000000*y3-.2500000000*y4)*xi*eta

    return xout,yout

