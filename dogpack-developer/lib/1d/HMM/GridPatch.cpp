// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (GridPatch.cpp)
//    A class that describes a 1d grid patch
// --------------------------------------------------------------------------

#include "../defs.h"
#include "GridPatch.h"

// Private members of class:
//    int mx,meqn,maux,kmax,mbc;
//    double xlow,xhigh;
//    dTensorBC3   q(mx,meqn,kmax,mbc);
//    dTensorBC3 aux(mx,maux,kmax,mbc);

GridPatch::GridPatch(int mx_in, double xlow_in, double xhigh_in,
		     int meqn_in, int maux_in, 
		     int kmax_in, int mbc_in) :   
  q(mx_in,meqn_in,kmax_in,mbc_in), aux(mx_in,maux_in,kmax_in,mbc_in)
// Constructor
{
    mx    = mx_in;
    xlow  = xlow_in;
    xhigh = xhigh_in;
    meqn  = meqn_in;
    maux  = maux_in;
    kmax  = kmax_in;
    mbc   = mbc_in;

    // Some simple error checking
    if (mx<1 || meqn<1 || maux<1 || kmax<1 || mbc<1)
    {
        cout << endl;
	cout << " Error in GridPatch Constructor ... " << endl;
	cout << "      mx = " << mx << endl;
	cout << "    meqn = " << meqn << endl;
	cout << "    maux = " << maux << endl;
	cout << "    kmax = " << kmax << endl;
	cout << "     mbc = " << mbc << endl;
	cout << endl;
	exit(1);
    }

    if (xlow>=xhigh)
    {
        cout << endl;
	cout << " Error in GridPatch Constructor ... " << endl;
	cout << "    xlow = " << xlow  << endl;
	cout << "   xhigh = " << xhigh << endl;
	cout << endl;
	exit(1);
    }
}

GridPatch::GridPatch(const GridPatch& another) :   q(mx,meqn,kmax,mbc), 
						   aux(mx,maux,kmax,mbc)
// Copy constructor
{
    mx    = another.mx;
    xlow  = another.xlow;
    xhigh = another.xhigh;
    meqn  = another.meqn;
    maux  = another.maux;
    kmax  = another.kmax;
    mbc   = another.mbc;
}

GridPatch::~GridPatch()
// Destructor
{
}

void GridPatch::getq(dTensorBC3& qout)
// Put the values of "q" into "qout"
{
    int i,m,k;
    int mx1   = qout.getsize(1);
    int meqn1 = qout.getsize(2);
    int kmax1 = qout.getsize(3);
    int mbc1  = qout.getmbc();
 
    if (mbc!=mbc1 || mx!=mx1 || meqn!=meqn1 || kmax!=kmax1)
    {
        cout << endl;
	cout << " Error in GridPatch::GetQ " << endl;
	cout << "  mbc = " << mbc  << "   mbc1 = " << mbc1  << endl;
	cout << "   mx = " << mx   << "    mx1 = " << mx1   << endl;
	cout << " meqn = " << meqn << "  meqn1 = " << meqn1 << endl;
	cout << " kmax = " << kmax << "  kmax1 = " << kmax1 << endl;
	cout << endl;
	exit(1);
    }

    for (i=(1-mbc); i<=(mx+mbc); i++)
      for (m=1; m<=meqn; m++)
	for (k=1; k<=kmax; k++)
	  {   qout.set(i,m,k, q.get(i,m,k) ); }
}

void GridPatch::getaux(dTensorBC3& auxout)
// Put the values of "aux" into "auxout"
{
    int i,m,k;
    int mx1   = auxout.getsize(1);
    int maux1 = auxout.getsize(2);
    int kmax1 = auxout.getsize(3);
    int mbc1  = auxout.getmbc();
 
    if (mbc!=mbc1 || mx!=mx1 || maux!=maux1 || kmax!=kmax1)
    { 
        cout << endl;
	cout << " Error in GridPatch::GetAux " << endl;
	cout << "  mbc = " << mbc  << "   mbc1 = " << mbc1  << endl;
	cout << "   mx = " << mx   << "    mx1 = " << mx1   << endl;
	cout << " maux = " << maux << "  maux1 = " << maux1 << endl;
	cout << " kmax = " << kmax << "  kmax1 = " << kmax1 << endl;
	cout << endl;
	exit(1);
    }

    for (i=(1-mbc); i<=(mx+mbc); i++)
      for (m=1; m<=maux; m++)
	for (k=1; k<=kmax; k++)
	  {   auxout.set(i,m,k, aux.get(i,m,k) ); }
}

void GridPatch::setq(dTensorBC3 qin)
// Put the values of "qin" into "q"
{
    int i,m,k;
    int mx1   = qin.getsize(1);
    int meqn1 = qin.getsize(2);
    int kmax1 = qin.getsize(3);
    int mbc1  = qin.getmbc();
 
    if (mbc!=mbc1 || mx!=mx1 || meqn!=meqn1 || kmax!=kmax1)
    {
        cout << endl;
	cout << " Error in GridPatch::SetQ " << endl;
	cout << "  mbc = " << mbc  << "   mbc1 = " << mbc1  << endl;
	cout << "   mx = " << mx   << "    mx1 = " << mx1   << endl;
	cout << " meqn = " << meqn << "  meqn1 = " << meqn1 << endl;
	cout << " kmax = " << kmax << "  kmax1 = " << kmax1 << endl;
	cout << endl;
	exit(1);
    }

    for (i=(1-mbc); i<=(mx+mbc); i++)
      for (m=1; m<=meqn; m++)
	for (k=1; k<=kmax; k++)
	  {   q.set(i,m,k, qin.get(i,m,k) ); }
}

void GridPatch::setaux(dTensorBC3 auxin)
// Put the values of "auxin" into "aux"
{
    int i,m,k;
    int mx1   = auxin.getsize(1);
    int maux1 = auxin.getsize(2);
    int kmax1 = auxin.getsize(3);
    int mbc1  = auxin.getmbc();
 
    if (mbc!=mbc1 || mx!=mx1 || maux!=maux1 || kmax!=kmax1)
    {
        cout << endl;
	cout << " Error in GridPatch::SetAux " << endl;
	cout << "  mbc = " << mbc  << "   mbc1 = " << mbc1  << endl;
	cout << "   mx = " << mx   << "    mx1 = " << mx1   << endl;
	cout << " maux = " << maux << "  maux1 = " << maux1 << endl;
	cout << " kmax = " << kmax << "  kmax1 = " << kmax1 << endl;
	cout << endl;
	exit(1);
    }

    for (i=(1-mbc); i<=(mx+mbc); i++)
      for (m=1; m<=maux; m++)
	for (k=1; k<=kmax; k++)
	  {   aux.set(i,m,k, auxin.get(i,m,k) ); }
}

double GridPatch::getxlow()
// get parameter xlow
{
  return xlow;
}

double GridPatch::getxhigh()
// get parameter xhigh
{
  return xhigh;
}

double GridPatch::getdx()
// get parameter dx
{
  return (xhigh-xlow)/double(mx);
}

int GridPatch::getmx()
// get parameter mx
{
  return mx;
}

void GridPatch::setxlow(double xlow_in)
// set parameter xlow
{
  xlow = xlow_in;
}

void GridPatch::setxhigh(double xhigh_in)
// set parameter xhigh
{
  xhigh = xhigh_in;
}

void GridPatch::resize(int mx_in)
// reset mx and modify q and aux
// according to this new mx
{
    mx = mx_in;
    
    q.resize(mx);
    aux.resize(mx);
}

int GridPatch::getmeqn()
// get parameter meqn
{
  return meqn;
}

int GridPatch::getmaux()
// get parameter maux
{
  return maux;
}

int GridPatch::getkmax()
// get parameter kmax
{
  return kmax;
}

int GridPatch::getmbc()
// get parameter mbc
{
  return mbc;
}
