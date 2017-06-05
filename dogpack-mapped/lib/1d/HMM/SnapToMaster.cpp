#include "../defs.h"
#include "GridPatch.h"

// Set the list "map" and the values "leftBC" and "rightBC"
// so that the slave grid and master grid can communicate
void SnapToMaster(GridPatch& Master, GridPatch& Slave,
		  iTensor1& map, iTensor2& leftBC, iTensor2& rightBC)
{
    int i;
    int mskip1 = 0;
    int mskip2 = 0;
    double dx,xlow,xhigh,dx_Master,xlow_Master,xhigh_Master;
    double x_edge,xlow_new,xhigh_new;
    int mstop;
    int istart = -100;
    int iend   = -100;
    int mbc = Master.getmbc();
    int mbc_sub = Slave.getmbc();

    // Determine the starting and ending
    // points of the Slave grid
    xlow  = Slave.getxlow();
    xhigh = Slave.getxhigh();
    xlow_Master  = Master.getxlow();
    xhigh_Master = Master.getxhigh();
    dx_Master    = Master.getdx();

    // basic error check
    if (xlow>=(xhigh_Master-1.0e-12) || xhigh<=(xlow_Master+1.0e-12))
    {
        cout << endl;
	cout << " Error in SnapToMaster.cpp ... " << endl;
	cout << "    xlow and xhigh are incorrectly set ... " << endl;
	cout << "             xlow = " << xlow << endl;
	cout << "            xhigh = " << xhigh << endl;
	cout << "      xlow_Master = " << xlow_Master << endl;
	cout << "     xhigh_Master = " << xhigh_Master << endl;	
	cout << endl;
	exit(1);
    }

    // if xlow is out-of-bounds, set to xlow_Master
    if (xlow<=(xlow_Master+1.0e-12))
    {  
	xlow_new = xlow_Master;  
	istart = 1;
	mskip1 = 1;
    }
    
    // if xhigh is out-of-bounds, set to xhigh_Master
    if (xhigh>=(xhigh_Master-1.0e-12))
    {
	xhigh_new = xhigh_Master;
	iend = Master.getmx();
	mskip2 = 1;
    }

    // find the new xlow so that it exists as
    // an edge on the Master grid    
    if (mskip1==0)
    {
        istart = 1;
        mstop = 0;
	
	x_edge = xlow_Master;

	while (mstop==0)
	{
	    if ( fabs(xlow-x_edge) <= (0.5*dx_Master+1.0e-12)  )
	    { 
	        xlow_new = x_edge;
		mstop = 1;
	    }
	    else
	    {
  	        x_edge = x_edge + dx_Master;
		istart = istart + 1;
	    }
	}
    }
    else
    {
	x_edge = xlow_Master;
    }

    // make sure there is enough space to allow
    // for ghost cells
    if (istart>1 && (istart-mbc_sub)<1)
    {
	istart = mbc_sub+1;
	xlow_new = xlow_Master + double(mbc_sub)*dx_Master;
	x_edge = xlow_new;
    }

    // find the new xhigh so that it exists as
    // an edge on the Master grid  
    if (mskip2==0)
    {
        iend = istart-1;
        mstop = 0;
	
	while (mstop==0)
	{
	  if ( fabs(xhigh - x_edge) < 0.5*dx_Master  )
	    { 
		mstop = 1;

		if (fabs(fabs(xhigh-x_edge)-0.5*dx_Master)<=1.0e-12)
		  {  
		    xhigh_new = x_edge + dx_Master;
		    iend = iend+1;
		  }
		else
		  {
		    xhigh_new = x_edge;
		  }
	    }
	    else
	    {
  	        x_edge = x_edge + dx_Master;
		iend = iend + 1;
	    }
	}
    }

    // make sure there is enough space to allow
    // for ghost cells
    if (iend<Master.getmx() && (Master.getmx()-iend)<mbc_sub)
    {
        iend = Master.getmx()-mbc_sub;
	xhigh_new = Master.getxhigh() - double(mbc_sub)*dx_Master;
    }

    // Modify "Slave" grid patch
    Slave.resize(iend+1-istart);
    Slave.setxlow(xlow_new);
    Slave.setxhigh(xhigh_new);

    // set the "map" array, which will be useful in 
    // communicating information from the slave grid
    // to the master grid
    map.resize(iend+1-istart);

    for (i=1; i<=map.getsize(); i++)
    {  map.set(i, i+istart-1 );  }

    // set the "leftBC" and "rightBC" arrays, which will
    // be useful in communicating boundary info from
    // master grid to the slave grid
    for (i=1; i<=mbc_sub; i++)
    {  
        leftBC.set(i,1,  map.get(1) - i );  
	rightBC.set(i,1, map.get(map.getsize()) + i );

	if (map.get(1)==1)
	  {  leftBC.set(i,2, -1 );  }
	else
	  {  leftBC.set(i,2,  1 );  }

	if (map.get(map.getsize())==Master.getmx())
	  {  rightBC.set(i,2, -1 );  }
	else
	  {  rightBC.set(i,2,  1 );  }
    }
    
}
