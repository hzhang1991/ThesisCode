import numpy as np

def mapc2p(xc,yc):
        """
    Specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
        """

        import string
        import os 
        from numpy import maximum,minimum,sqrt,ones,sign
        #p = os.system('ls')
        datapre=open('setprob.data','r')
        data1=datapre.readlines()
#print data1.split()

        data=[float(da) for da1 in data1 for da in sum([da1.split()],[])]

        #print data

	#print "here",data

	zout = data[0];
	cout = data[1];
	ncirc = int(data[2]);
	x0circ=list (range(ncirc))
	y0circ=list (range(ncirc))
	r1circ=list (range(ncirc))
	r2circ=list (range(ncirc))
	zin=list (range(ncirc))
	cin=list (range(ncirc))
	for i in range(ncirc):
	   x0circ[i]= data[3+6*(i)];
	   y0circ[i]= data[4+6*(i)];
	   r1circ[i] = data[5+6*(i)];
	   r2circ[i] = data[6+6*(i)];
	   zin[i] = data[7+6*(i)];
	   cin[i] = data[8+6*(i)];



	xp = np.array(xc);
	yp = np.array(yc);


	for i in range(ncirc):
	  x0 = x0circ[i];
	  y0 = y0circ[i];
	  r1 = r1circ[i];
	  r2 = r2circ[i];
	  xc0 = abs(xc-x0);
	  yc0 = abs(yc-y0);
          ij1 = (xc0<=r2)*(yc0<=r2)
#	  ij1 = (max(xc0,yc0)<=r2);  # portion of grid deformed for this circle
	  xc0 = (xc[ij1] - x0)/r2;       # subset of gridpoints on this portion
	  yc0 = (yc[ij1] - y0)/r2;
	  
	  xc1 = abs(xc0);
	  yc1 = abs(yc0);
	  d = maximum(xc1,yc1);            # value on diagonal of computational grid
	  d = maximum(d, 1e-14);           # to avoid divide by zero at center
	  d = minimum(d, 0.99999);         # to avoid divide by zero at d=1
	  d1 = d*r2/sqrt(2.);
	  R = sqrt(2.) * d1;
	  R = r1*ones(d1.shape);
	  #R = r1^2 ./ (r2*d);
	  
	  
	  # modify d1 and R outside circle to morph back to square:
	  ij = (d>r1/r2);
	  d1[ij] = r1/sqrt(2.) + (d[ij]-r1/r2)*(r2-r1/sqrt(2.0))/(1.-r1/r2);
          #print 'problem=','d=',d,'x1=',xc1,'yc1=',yc1,d[ij]
	  R[ij] = r1 * ((1.-r1/r2) / (1.-d[ij]))**(r2/r1 + .5);
	  
	  xp2 = d1/d * xc1;
	  yp2 = d1/d * yc1;
	  center = d1 - sqrt(R**2 - d1**2);
	  
	  ij = (xc1 > yc1-1.0e-15);
	  xp2[ij] = center[ij] + sqrt(R[ij]**2 - yp2[ij]**2);
	  
	  ij = (yc1 > xc1-1.0e-15);
	  yp2[ij] = center[ij] + sqrt(R[ij]**2 - xp2[ij]**2);
	  
	  xp2 = sign(xc0) * xp2;
	  yp2 = sign(yc0) * yp2;
	  
	  xp[ij1] = x0 + xp2;      # grid is modified only on portion for this circle
	  yp[ij1] = y0 + yp2;
        xpm=np.array(xc);ypm=np.array(yc);
  
        return xp,yp
