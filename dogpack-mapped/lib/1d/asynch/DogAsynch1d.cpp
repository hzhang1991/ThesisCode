// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (DogAsynch1d.cpp)
//    1d unstructured mesh
// --------------------------------------------------------------------------

#include "DogAsynch1d.h"
#include "dog_math.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
#include <assert.h>
using namespace std;

DogAsynch1d::DogAsynch1d(int inmx, int inmeqn, int inmaux, int inkmax, int inmbc)
// Constructor
// POST: Creates a mesh
{
  mx   = inmx;
  meqn = inmeqn;
  maux = inmaux;
  kmax = inkmax;
  mbc  = inmbc;

  if (mx<1 || meqn<1 || maux<0 || kmax<1 || mbc<0)
    {
      cout << endl;
      cout << " Error in mesh constructor ... "       << endl;
      cout << "   mx = " << mx   << endl;
      cout << " meqn = " << meqn << endl;
      cout << " maux = " << maux << endl;
      cout << " kmax = " << kmax << endl;
      cout << "  mbc = " << mbc  << endl;
      cout << endl;
      exit(1);
    }

  // tstart and tend (these must be reset later)
  tstart = 0.0;
  tend   = 0.0;
  dt_max = 0.0;

  // node: list of x coordinates of nodes
  node = new dTensorBC1(mx+1,mbc);
  
  // element lengths
  dx = new dTensorBC1(mx,mbc);

  // element times
  time = new dTensorBC1(mx,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    { time->set(i, 0.0 ); }

  time_old = new dTensorBC1(mx,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    { time_old->set(i, 0.0 ); }
    
  // element time-steps
  dt = new dTensorBC1(mx,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    { dt->set(i, 0.0 ); }

  // solution array
  q = new dTensorBC3(mx,meqn,kmax,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  q->set(i,m,k, 0.0 );
	}

  // old solution array
  qold = new dTensorBC3(mx,meqn,kmax,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  qold->set(i,m,k, 0.0 );
	}

  // auxiliary variable array  
  aux = new dTensorBC3(mx,maux,kmax,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  aux->set(i,m,k, 0.0 );
	}

  // old auxiliary variable array
  auxold = new dTensorBC3(mx,maux,kmax,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  auxold->set(i,m,k, 0.0 );
	}

  // vector of (2*u/dx)^(k-1) and (-2*u/dx)^(k-1) for each element i
  facp = new dTensorBC2(mx,kmax,mbc);
  facm = new dTensorBC2(mx,kmax,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int k=1; k<=kmax; k++)
      {
	facp->set(i,k, 0.0 );
	facm->set(i,k, 0.0 );
      }

  // list of update elements
  num_to_be_updated = mx;
  curr_to_be_updated = 1;
  elements_to_be_updated = new iTensor1(mx);

  // Correction matrices
  Mmat = new dTensor2(meqn,kmax);     // (meqn,1..kmax)
  Rmat = new dTensor2(meqn,kmax);     // (meqn,1..kmax)
  Lmat = new dTensor2(meqn,kmax);     // (meqn,1..kmax)
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	Mmat->set(m,k, 0.0 );
	Rmat->set(m,k, 0.0 );
	Lmat->set(m,k, 0.0 );
      }

  // Helpers for correction matrices
  Mmat_s = new dTensor3(meqn,kmax,kmax-1);   // (meqn,kmax,kmax-1);
  Ltmp = new dTensor2(meqn,kmax);            // (meqn,kmax)
  Rtmp = new dTensor2(meqn,kmax);            // (meqn,kmax)
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	for (int s=1; s<=(kmax-1); s++)
	  {  Mmat_s->set(m,k,s, 0.0 );  }
	Ltmp->set(m,k, 0.0 );
	Rtmp->set(m,k, 0.0 );
      }

  // Time derivative matrix
  Tmat = new dTensorBC4(mx,meqn,kmax,kmax-1,mbc);
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	for (int s=1; s<=(kmax-1); s++)
	  {
	    Tmat->set(i,m,k,s, 0.0 );	    
	  }

  // matrix of phi derivatives at xi=-1 (left) and xi=1 (right)
  dphiR = new dTensor2(5,5);    // (l=1..5,s+1=1..5)
  dphiL = new dTensor2(5,5);    // (l=1..5,s+1=1..5)
  
  dphiR->set(1,1,  1.0  );
  dphiR->set(2,1,  sq3  );
  dphiR->set(3,1,  sq5  );
  dphiR->set(4,1,  sq7  );
  dphiR->set(5,1,  3.0  );

  dphiR->set(1,2,  0.0  );
  dphiR->set(2,2,  sq3  );
  dphiR->set(3,2,  3.0*sq5  );
  dphiR->set(4,2,  6.0*sq7  );
  dphiR->set(5,2,  30.0  );

  dphiR->set(1,3,  0.0  );
  dphiR->set(2,3,  0.0  );
  dphiR->set(3,3,  3.0*sq5  );
  dphiR->set(4,3,  15.0*sq7  );
  dphiR->set(5,3,  135.0  );

  dphiR->set(1,4,  0.0  );
  dphiR->set(2,4,  0.0  );
  dphiR->set(3,4,  0.0  );
  dphiR->set(4,4,  15.0*sq7  );
  dphiR->set(5,4,  315.0  );

  dphiR->set(1,5,  0.0  );
  dphiR->set(2,5,  0.0  );
  dphiR->set(3,5,  0.0  );
  dphiR->set(4,5,  0.0  );
  dphiR->set(5,5,  315.0  );

  dphiL->set(1,1,  1.0  );
  dphiL->set(2,1, -sq3  );
  dphiL->set(3,1,  sq5  );
  dphiL->set(4,1, -sq7  );
  dphiL->set(5,1,  3.0  );

  dphiL->set(1,2,  0.0  );
  dphiL->set(2,2,  sq3  );
  dphiL->set(3,2, -3.0*sq5  );
  dphiL->set(4,2,  6.0*sq7  );
  dphiL->set(5,2, -30.0  );

  dphiL->set(1,3,  0.0  );
  dphiL->set(2,3,  0.0  );
  dphiL->set(3,3,  3.0*sq5  );
  dphiL->set(4,3, -15.0*sq7  );
  dphiL->set(5,3,  135.0  );

  dphiL->set(1,4,  0.0  );
  dphiL->set(2,4,  0.0  );
  dphiL->set(3,4,  0.0  );
  dphiL->set(4,4,  15.0*sq7  );
  dphiL->set(5,4, -315.0  );

  dphiL->set(1,5,  0.0  );
  dphiL->set(2,5,  0.0  );
  dphiL->set(3,5,  0.0  );
  dphiL->set(4,5,  0.0  );
  dphiL->set(5,5,  315.0  );

  // Stored fluxes at interfaces --> this is need for conservation fix-ups
  Flft = new dTensor3(mx+1,meqn,kmax);  // (mx+1,meqn,kmax)
  Frgt = new dTensor3(mx+1,meqn,kmax);  // (mx+1,meqn,kmax)
  ftally_lft = new iTensor1(mx+1);      // (mx+1)
  ftally_rgt = new iTensor1(mx+1);      // (mx+1)

  for (int i=1; i<=(mx+1); i++)
    {
      ftally_lft->set(i, 0 );
      ftally_rgt->set(i, 0 );

      for (int m=1; m<=meqn; m++)    
	for (int k=1; k<=kmax; k++)
	  {
	    Flft->set(i,m,k, 0.0 );
	    Frgt->set(i,m,k, 0.0 );
	  }
    }

  // vector of u^s
  us = new double[2*kmax];
  us[0] = 1.0;
  for (int s=1; s<=(2*kmax-1); s++)
    {  us[s] = 0.0;  } 

  // vector of dt^s
  dts = new double[2*kmax];
  dts[0] = 1.0;
  for (int s=1; s<=(2*kmax-1); s++)
    {  dts[s] = 0.0;  }
}

DogAsynch1d::DogAsynch1d(const DogAsynch1d& amesh)
// Copy constructor
// POST: New mesh created with size and contents same as amesh
{
  mx   = amesh.mx;
  meqn = amesh.meqn;
  maux = amesh.maux;
  kmax = amesh.kmax;
  mbc  = amesh.mbc;

  CFL     = amesh.CFL;
  CFL_MAX = amesh.CFL_MAX;
  tstart  = amesh.tstart;
  tend    = amesh.tend;
  dt_max  = amesh.dt_max;
  
  for (int i=(1-mbc); i<=(mx+1+mbc); i++)
    {
      node->set(i, amesh.node->get(i) );
    }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      dx->set(i, amesh.dx->get(i) );
    }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      time->set(i, amesh.time->get(i) );
    }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      time_old->set(i, amesh.time_old->get(i) );
    }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      dt->set(i, amesh.dt->get(i) );
    }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  q->set(i,m,k, amesh.q->get(i,m,k) );
	}

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  qold->set(i,m,k, amesh.qold->get(i,m,k) );
	}

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  aux->set(i,m,k, amesh.aux->get(i,m,k) );
	}

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  auxold->set(i,m,k, amesh.auxold->get(i,m,k) );
	}

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int k=1; k<=kmax; k++)
      {
	facp->set(i,k, amesh.facp->get(i,k) );
	facm->set(i,k, amesh.facm->get(i,k) );
      }

  num_to_be_updated = amesh.num_to_be_updated;
  curr_to_be_updated = amesh.curr_to_be_updated;
  for (int i=1; i<=mx; i++)
    {
      elements_to_be_updated->set(i, amesh.elements_to_be_updated->get(i) );
    }
  
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	Mmat->set(m,k, amesh.Mmat->get(m,k) );
	Rmat->set(m,k, amesh.Rmat->get(m,k) );
	Lmat->set(m,k, amesh.Lmat->get(m,k) );
      }

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	for (int s=1; s<=(kmax-1); s++)
	  {
	    Tmat->set(i,m,k,s, amesh.Tmat->get(i,m,k,s) );	    
	  }

  for (int m=1; m<=5; m++)
    for (int k=1; k<=5; k++)
      {
	dphiR->set(m,k, amesh.dphiR->get(m,k) );
	dphiL->set(m,k, amesh.dphiL->get(m,k) );
      }

  for (int i=1; i<=(mx+1); i++)
    {
      ftally_lft->set(i, amesh.ftally_lft->get(i) );
      ftally_rgt->set(i, amesh.ftally_rgt->get(i) );
      
      for (int m=1; m<=meqn; m++)    
	for (int k=1; k<=kmax; k++)
	  {
	    Flft->set(i,m,k, amesh.Flft->get(i,m,k) );
	    Frgt->set(i,m,k, amesh.Frgt->get(i,m,k) );
	  }
    }

  for (int s=0; s<=(2*kmax-1); s++)
    {  us[s] = amesh.us[s];  } 
}

DogAsynch1d::~DogAsynch1d()
// Destructor
// POST: mesh no longer exists
{
}

void DogAsynch1d::Output(string outputdir,double t,int nframe)
// Ouput all mesh, q, and aux information
{
  OutputMesh(outputdir);
  OutputSoln(outputdir,t,nframe);
}

void DogAsynch1d::OutputMesh(string outputdir)
// Ouput mesh information
{
  string fname;

  // output: constants
  fname = outputdir + "/mesh_params.dat";
  ofstream write1(fname.c_str(), ios::out);
  write1 << setw(8) << mx  << "  : mx  "  << endl;
  write1 << setw(8) << mbc << "  : mbc "  << endl;
  write1.close();
  write1.clear();
  
  // output: NODE
  fname = outputdir + "/mesh_node.dat";
  ofstream write2(fname.c_str(), ios::out);
  write2 << setprecision(16); 
  for (int i=1; i<=(mx+1); i++)
    {
      write2 << setw(24) << scientific << node->get(i) << endl;
    }
  write2.close();
  write2.clear();
  
  // output: DX
  fname = outputdir + "/mesh_dx.dat";
  ofstream write3(fname.c_str(), ios::out);
  write3 << setprecision(16); 
  for (int i=1; i<=mx; i++)
    {
      write3 << setw(24) << scientific << dx->get(i) << endl;
    }
  write3.close();
  write3.clear();

  // output: DT
  fname = outputdir + "/mesh_dt.dat";
  ofstream write4(fname.c_str(), ios::out);
  write4 << setprecision(16); 
  for (int i=1; i<=mx; i++)
    {
      write4 << setw(24) << scientific << dt->get(i) << endl;
    }
  write4.close();
  write4.clear();

  // output: TIME
  fname = outputdir + "/mesh_time.dat";
  ofstream write5(fname.c_str(), ios::out);
  write5 << setprecision(16); 
  for (int i=1; i<=mx; i++)
    {
      write5 << setw(24) << scientific << time->get(i) << endl;
    }
  write5.close();
  write5.clear();

}

void DogAsynch1d::OutputSoln(string outputdir,double t,int nframe)
// Ouput q and aux arrays
{
  // Open file -- q
  ostringstream fname1;
  fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
  ofstream q_file(fname1.str().c_str(), ios::out );
   
  q_file << setprecision(16);
  q_file << setw(24) << scientific << t << endl;
  
  // Output each coefficient
  for (int k=1; k<=kmax; k++)
    for (int m=1; m<=meqn; m++)
      for (int i=1; i<=mx; i++)      
        {
          double tmp = q->get(i,m,k);
          
          if (fabs(tmp) < 1.0e-99) {tmp = 0.0;}
          q_file << setw(24) << scientific << tmp << endl;
        }
  q_file.close();

  // Open file -- aux
  ostringstream fname2;
  fname2 << outputdir << "/" << "a" << setfill('0') 
        << setw(4) << nframe << ".dat";
  ofstream aux_file(fname2.str().c_str(), ios::out );
   
  aux_file << setprecision(16);
  aux_file << setw(24) << scientific << t << endl;
  
  // Output each coefficient
  for (int k=1; k<=kmax; k++)
    for (int m=1; m<=maux; m++)
      for (int i=1; i<=mx; i++)      
        {
          double tmp = aux->get(i,m,k);
          
          if (fabs(tmp) < 1.0e-99) {tmp = 0.0;}
          aux_file << setw(24) << scientific << tmp << endl;
        }
  aux_file.close();
}

void DogAsynch1d::InputMesh(double xlow, double xhigh)
// Create uniform mesh
{
  double dxtmp = (xhigh-xlow)/double(mx);

  for (int i=(1-mbc); i<=(mx+1+mbc); i++)
    {
      node->set(i, xlow + (double(i)-1.0)*dxtmp );
    } 

  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      dx->set(i, dxtmp );
    }
}

void DogAsynch1d::InputMesh(string inputdir)
// Input all mesh information
{
  string fname;
  char buffer[256];
  int mx_in, mbc_in;
  double tmp_double;

  // input: constants
  fname = inputdir + "/mesh_params.dat";
  ifstream read1(fname.c_str(), ios::in);
  read1 >> mx_in;   read1.getline(buffer,256);
  read1 >> mbc_in;  read1.getline(buffer,256);
  read1.close();   
  read1.clear();
  
  // input: NODE
  fname = inputdir + "/mesh_node.dat";
  ifstream read2(fname.c_str(), ios::in);
  for (int i=1; i<=(mx+1); i++)
    {
      read2 >> tmp_double;
      node->set(i, tmp_double );
    }
  read2.close();
  read2.clear();
  
  // input: DX
  fname = inputdir + "/mesh_dx.dat";
  ifstream read3(fname.c_str(), ios::in);
  for (int i=1; i<=mx; i++)
    {
      read3 >> tmp_double;
      dx->set(i, tmp_double);
    }
  read3.close();
  read3.clear();

  // input: DT
  fname = inputdir + "/mesh_dt.dat";
  ifstream read4(fname.c_str(), ios::in);
  for (int i=1; i<=mx; i++)
    {
      read4 >> tmp_double;
      dt->set(i, tmp_double );
    }
  read4.close();
  read4.clear();

  // input: TIME
  fname = inputdir + "/mesh_time.dat";
  ifstream read5(fname.c_str(), ios::in);
  for (int i=1; i<=mx; i++)
    {
      read5 >> tmp_double;
      time->set(i, tmp_double );
    }
  read5.close();
  read5.clear();
}

// Returns "CFL"
const double& DogAsynch1d::get_CFL() const
{
  return CFL;
}

// Sets "CFL"
void DogAsynch1d::set_CFL(double in_CFL)
{
  CFL = in_CFL;
}

// Returns "CFL_MAX"
const double& DogAsynch1d::get_CFL_MAX() const
{
  return CFL_MAX;
}

// Sets "CFL_MAX"
void DogAsynch1d::set_CFL_MAX(double in_CFL_MAX)
{
  CFL_MAX = in_CFL_MAX;
}

const int& DogAsynch1d::get_mx() const
// Returns "mx"
{
  return mx;
}

const int& DogAsynch1d::get_meqn() const
// Returns "meqn"
{
  return meqn;
}

const int& DogAsynch1d::get_maux() const
// Returns "maux"
{
  return maux;
}

const int& DogAsynch1d::get_kmax() const
// Returns "kmax"
{
  return kmax;
}

const int& DogAsynch1d::get_mbc() const
// Returns "mbc"
{
  return mbc;
}

const double& DogAsynch1d::get_tstart() const
// Returns "tstart"
{
  return tstart;
}

const double& DogAsynch1d::get_tend() const
// Returns "tend"
{
  return tend;
}

const int& DogAsynch1d::get_num_to_be_updated() const
// Returns "num_to_be_updated"
{
  return num_to_be_updated;
}

void DogAsynch1d::dx_compute()
// Set the correct values of dx once node values have been set
{
  for (int i=1; i<=mx; i++)
    {
      double tmp = node->get(i+1)-node->get(i);
      if (tmp<=0.0)
	{
	  cout << endl;
	  cout << " Error in DogAsynch1d::dx_compute() " << endl;
	  cout << "   i = " << i << "    dx->get(i) = " << tmp << endl;
	  cout << "       node->get(i+1) = " << node->get(i+1) << endl;
	  cout << "       node->get(i)   = " << node->get(i) << endl;
	  cout << endl;
	  exit(1);
	}
      else
	{
	  dx->set(i, tmp );
	}
    }
}

// ---------------------------------------------
//  KEY FUNCTIONS FOR ASYNCHRONOUS TIME-STEPPING
// ---------------------------------------------

// Compute the time step, dt, that will
//       give the slowest moving element a courant
//       number of CFL.
double DogAsynch1d::get_max_dt_step()
{
  double u = fabs(aux->get(1,1,1));
  double dt_out = CFL*dx->get(1)/u;

  for (int i=2; i<=mx; i++)
    {
      dt_out = Max(dt_out, CFL*dx->get(i)/u );
    }

  return dt_out;
}

// Compute the time step, dt, that will
//       give the fastest moving element a courant
//       number of CFL.
double DogAsynch1d::get_min_dt_step()
{
  double u = fabs(aux->get(1,1,1));
  double dt_out = CFL*dx->get(1)/u;

  for (int i=2; i<=mx; i++)
    {
      dt_out = Min(dt_out, CFL*dx->get(i)/u );
    }

  return dt_out;
}


void DogAsynch1d::initialize_for_time_stepping(double in_tstart, double in_tend)
// Initialize mesh for time-stepping
{
  // make sure tend>tstart
  assert(in_tend>in_tstart);
  tstart = in_tstart;
  tend   = in_tend;
  dt_max = tend - tstart;

  // set all element times to tstart
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      time->set(i, tstart );  
      time_old->set(i, tstart );
    }

  // create list of elements that will eventually be updated
  num_to_be_updated = mx;
  curr_to_be_updated = 1;
  for (int i=1; i<=mx; i++)
    {  elements_to_be_updated->set(i, i );  }

  // set qold values to q
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=meqn; m++)
      for (int k=1; k<=kmax; k++)
	{
	  qold->set(i,m,k, q->get(i,m,k) );
	}

  // set auxold values to aux
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  auxold->set(i,m,k, aux->get(i,m,k) );
	}

  // set auxold values to aux
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	for (int s=1; s<=(kmax-1); s++)
	  {
	    Tmat->set(i,m,k,s, 0.0 );
	  }	

  // set boundary conditions
  SetBndValues(1);
  SetBndValues(mx);

  // compute derivative mapping factors
  for (int i=0; i<=(mx+1); i++)
    {
      double tmp = 1.0;
      facp->set(i,1, tmp );
      for (int k=2; k<=kmax; k++)
	{
	  tmp = tmp*(2.0/dx->get(i));
	  facp->set(i,k, tmp );
	}
      
      tmp = 1.0;
      facm->set(i,1, tmp );
      for (int k=2; k<=kmax; k++)
	{
	  tmp = tmp*(-2.0/dx->get(i));
	  facm->set(i,k, tmp );
	}
    }

  // initial flux fix-ups
  for (int i=1; i<=(mx+1); i++)
    {
      ftally_lft->set(i, 0 );
      ftally_rgt->set(i, 0 );

      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    Flft->set(i,m,k, 0.0 );
	    Frgt->set(i,m,k, 0.0 );
	  }
    }
}

void DogAsynch1d::check_if_done(int i)
// Check to see if input element is at final time
// If   NO: do nothing
// If  YES: remove input element from "elements_to_be_updated" list
{
  double t = time->get(i);
  if (fabs(t-tend)<=1.0e-15)
    {
      // switch current element with last one in the "elements_to_be_updated" list      
      if (curr_to_be_updated==num_to_be_updated)
	{
	  curr_to_be_updated = 1;
	}
      else
	{
	  int a = elements_to_be_updated->get(curr_to_be_updated);
	  elements_to_be_updated->set( curr_to_be_updated, elements_to_be_updated->get(num_to_be_updated) );
	  elements_to_be_updated->set( num_to_be_updated, a );
	}

      // decrement the number of elements that still need to be updated
      num_to_be_updated = num_to_be_updated - 1;      
    }
}

bool DogAsynch1d::is_updatable(int i) const
// Returns 1:  if element i is at time <= its neighbors
// Returns 0:  otherwise
{
  double tl = time->get(i-1);
  double tm = time->get(i);
  double tr = time->get(i+1);
  bool ready;

  if ( (tm<=tl) && (tm<=tr) )
    {  ready = 1;  }
  else
    {  ready = 0;  }

  return ready;
}


int DogAsynch1d::get_next_updatable_element()
// Find the next element that should be updated
{
  // check if done
  if (num_to_be_updated == 0)
    {
      return 0;
    }
  else
    {
      // check if current element is updatable
      int current = elements_to_be_updated->get(curr_to_be_updated);
      if(is_updatable(current))
	{
	  return current;
	}
      else
	{
	  // find next updatable element
	  curr_to_be_updated = curr_to_be_updated%num_to_be_updated + 1;
	  current = elements_to_be_updated->get(curr_to_be_updated);
	  
	  while(!is_updatable(current))
	    {
	      curr_to_be_updated = curr_to_be_updated%num_to_be_updated + 1;
	      current = elements_to_be_updated->get(curr_to_be_updated);
	    }

	  return current;
	}
    }
}


void DogAsynch1d::LxW_update(int i)
// Lax-Wendroff method to take solution from  t = tstart  and  t = tend
{
  // Compute max wave speed in current element
  double u = aux->get(i,1,1);

  // Raise velocity to powers
  for (int s=1; s<=(2*kmax-1); s++)
    {  us[s] = us[s-1]*u;  } 

  // Calculate time-step
  dt_max = tend - time->get(i);
  double dt_nice = Min( time->get(i+1)-time->get(i), time->get(i-1)-time->get(i) );
  double CFL_nice = fabs(u)*dt_nice/dx->get(i);
  if (CFL_nice>=(0.95*CFL) && CFL_nice<=CFL_MAX)
    {  dt->set(i, dt_nice );  }
  else
    {  dt->set(i, Min(CFL*(dx->get(i)/fabs(u)), dt_max) );  }

  // Raise dt to powers  
  double dti = dt->get(i);
  for (int s=1; s<=(2*kmax-1); s++)
    {  dts[s] = dts[s-1]*dti;  }

  // Construct time derivatives and store in T
  ConstructMRL(i,u);

  // Update solution
  UpdateSoln(i);

  // Set boundary conditions
  if (i==1 || i==mx) 
    { SetBndValues(i); }

  // Conservation fix
  //ConservationFix(i,u);
  
  // Check if current element is at the final time
  check_if_done(i);
}

void DogAsynch1d::ConstructMRL(int i, double u)
// Construct time derivatives and store in T
{
  // Compute internal part of update
  dTensor2 Qin(meqn,kmax);
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	Qin.set(m,k, q->get(i,m,k) );
      }
  SetMmat(i,u,Qin);

  // Compute boundary fluxes
  EdgeFlux(i,u,"left ");
  EdgeFlux(i,u,"right");

  // Compute left and right part of update
  SetLmat(i);
  SetRmat(i);
}


void DogAsynch1d::EdgeFlux(int i, double u, string LR)
// Get left edge values for element i
{
  dTensor2 DQl(meqn,kmax),DQr(meqn,kmax);
  dTensor2  Ql(meqn,kmax), Qr(meqn,kmax);
  double told,tnew,a,a2,a3,a4,a5;  

  // Ql and Qr
  SetQstates(i,Ql,Qr,LR);
    
  // Compute spatial derivatives of Ql and Qr
  ComputeDerivatives(Ql,Qr,DQl,DQr);

  // Solve Riemann problem
  RiemannSolve(i,u,DQl,DQr,LR);
}

// Compute left part of update
void DogAsynch1d::SetLmat(int i)
{
  // left flux
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	double tmp = 0.0;

	for (int s=k; s<=(kmax+k-1); s++)
	  {
	    tmp = tmp + us[s]*dts[s]/factorial[s]*Ltmp->get(m,s-k+1);
	  }

	Frgt->set(i,m,k, Frgt->get(i,m,k) + tmp );
      }

  // Lmat
  for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
      {
	double tmp = 0.0;
	for (int k=1; k<=kmax; k++)
	  {
	    tmp = tmp + facp->get(i,k)*dphiL->get(ell,k)*Frgt->get(i,m,k);
	  }
	Lmat->set(m,ell, tmp/dx->get(i) );
      }
}

// Compute right part of update
void DogAsynch1d::SetRmat(int i)
{
  // right flux
  for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
      {
	double tmp = 0.0;

	for (int s=k; s<=(kmax+k-1); s++)
	  {
	    tmp = tmp + us[s]*dts[s]/factorial[s]*Rtmp->get(m,s-k+1);
	  }

	Flft->set(i+1,m,k, Flft->get(i+1,m,k) + tmp );
      }

  
  // Rmat
  for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
      {
	double tmp = 0.0;
	for (int k=1; k<=kmax; k++)
	  {
	    tmp = tmp + facp->get(i,k)*dphiR->get(ell,k)*Flft->get(i+1,m,k);
	  }
	Rmat->set(m,ell, -tmp/dx->get(i) );
      }
}

// Compute internal part of update
void DogAsynch1d::SetMmat(int i,double u,dTensor2& Qin)
{
  switch(kmax)
    {
    case 2:
      for (int m=1; m<=meqn; m++)
	{
	  Mmat_s->set(m,2,1, sq3*Qin.get(m,1) );
	}
      break;

    case 3:
      for (int m=1; m<=meqn; m++)
	{
	  Mmat_s->set(m,2,1,  sq3*Qin.get(m,1) );
	  Mmat_s->set(m,3,1,  sq3*sq5*Qin.get(m,2) );	 
	  Mmat_s->set(m,3,2,  3.0*sq5*Qin.get(m,1) );
	}
      break;

    case 4:
      for (int m=1; m<=meqn; m++)
	{
	  Mmat_s->set(m,2,1,  sq3*Qin.get(m,1) );
	  Mmat_s->set(m,3,1,  sq3*sq5*Qin.get(m,2) );	 
	  Mmat_s->set(m,3,2,  3.0*sq5*Qin.get(m,1) );
	  Mmat_s->set(m,4,1,  sq5*sq7*Qin.get(m,3)+sq7*Qin.get(m,1) );
	  Mmat_s->set(m,4,2,  5.0*sq7*Qin.get(m,2)*sq3 );
	  Mmat_s->set(m,4,3,  15.0*sq7*Qin.get(m,1) );
	}
      break;

    case 5:
      for (int m=1; m<=meqn; m++)
	{
	  Mmat_s->set(m,2,1,  sq3*Qin.get(m,1) );
	  Mmat_s->set(m,3,1,  sq3*sq5*Qin.get(m,2) );	 
	  Mmat_s->set(m,3,2,  3.0*sq5*Qin.get(m,1) );
	  Mmat_s->set(m,4,1,  sq5*sq7*Qin.get(m,3)+sq7*Qin.get(m,1) );
	  Mmat_s->set(m,4,2,  5.0*sq7*Qin.get(m,2)*sq3 );
	  Mmat_s->set(m,4,3,  15.0*sq7*Qin.get(m,1) );
	  Mmat_s->set(m,5,1,  3.0*sq7*Qin.get(m,4)+3.0*sq3*Qin.get(m,2) );	  
	  Mmat_s->set(m,5,2,  21.0*sq5*Qin.get(m,3)+30.0*Qin.get(m,1) );
	  Mmat_s->set(m,5,3,  105.0*Qin.get(m,2)*sq3 );	  
	  Mmat_s->set(m,5,4,  315.0*Qin.get(m,1) );
	}
      break;
    }

  for (int m=1; m<=meqn; m++)
    for (int k=2; k<=kmax; k++)
      {
	double tmp = 0.0;
	for (int s=1; s<=(k-1); s++)
	  {
	    tmp = tmp + dts[s]/factorial[s]*us[s]*facp->get(i,s+1)*Mmat_s->get(m,k,s);
	  }
	Mmat->set(m,k, tmp );
      }
}

void DogAsynch1d::UpdateSoln(int i)
// Update solution
{
  for (int m=1; m<=meqn; m++)
    for (int ell=1; ell<=kmax; ell++)
      {
	double tmp = q->get(i,m,ell);
	qold->set(i,m,ell, tmp );

	tmp = tmp + (Mmat->get(m,ell) + Rmat->get(m,ell) + Lmat->get(m,ell));
	q->set(i,m,ell, tmp );
      }

  time_old->set(i, time->get(i) );
  time->set(i, time->get(i) + dt->get(i) );
}

void DogAsynch1d::ConservationFix(int i, double u)
{
  double tl_new = time->get(i-1);
  double tm_new = time->get(i);
  double tr_new = time->get(i+1);

  double tl_old = time_old->get(i-1);
  double tm_old = time_old->get(i);
  double tr_old = time_old->get(i+1);

  /*
  // interface i-1/2 is ready for fixup
  if (fabs(tm_new-tl_new)<=1.0e-15 && fabs(tm_old-tl_old)>=1.0e-15)
    {      
      cout << endl << endl << " i = " << i;
      cout << "   ftally_rgt->get(i) = " << ftally_rgt->get(i);
      cout << "   ftally_lft->get(i) = " << ftally_lft->get(i);
      cout << endl;

      // Correct fluxes
      if (ftally_rgt->get(i)==1 && ftally_lft->get(i)>1)
	{
	FluxCorrect(i-1,u,"right");
	}
      else if (ftally_rgt->get(i)>1 && ftally_lft->get(i)==1)
	{
	FluxCorrect(i,u,"left ");
	}

      // Reset tallies
      ftally_rgt->set(i, 0 );
      ftally_lft->set(i, 0 );

      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    flux_lft->set(i,m,k, 0.0 );
	    flux_rgt->set(i,m,k, 0.0 );
	  }
    }
  */

  // interface i+1/2 is ready for fixup
  if (fabs(tm_new-tr_new)<=1.0e-15 && fabs(tm_old-tr_old)>=1.0e-15)
    {
      cout << endl << endl << " i+1 = " << i+1;
      cout << "   ftally_rgt->get(i+1) = " << ftally_rgt->get(i+1);
      cout << "   ftally_lft->get(i+1) = " << ftally_lft->get(i+1);
      cout << endl;

      // Correct fluxes
      if (ftally_lft->get(i+1)==1 && ftally_rgt->get(i+1)>1)
	{
	  //FluxCorrect(i,u,"right");
	  exit(1);
	}
      else if (ftally_lft->get(i+1)>1 && ftally_rgt->get(i+1)==1)
	{	  
	  //FluxCorrect(i+1,u,"left ");
	}

      // Reset tallies
      ftally_rgt->set(i+1, 0 );
      ftally_lft->set(i+1, 0 );

      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    Flft->set(i+1,m,k, 0.0 );
	    Frgt->set(i+1,m,k, 0.0 );
	  }
    }

}


// Correct flux at interface
void DogAsynch1d::FluxCorrect(int i, double u, string LR)
{
  cout << " i = " << i << "   LR = " << LR << endl;

  if (LR=="right")
    {
       
    }
  else if (LR=="left ")
    { 
      /*     
      // Correct time derivatives  
      for (int m=1; m<=meqn; m++)    
	for (int ell=1; ell<=kmax; ell++)
	  {
	    double tmp = 0.0;
	    for (int k=1; k<=kmax; k++)
	      {
		us = us*u;
		
		tmp = tmp - us[s]*fac->get(i,k)*dphiL->get(ell,k)*
		  ( flux_lft->get(i,m,k) - flux_rgt->get(i,m,k));
	      }	    
	     q->set(i,m,ell, q->get(i,m,ell) + tmp );
	  }
      */
    }

  if (i==1 || i==mx)
    {  SetBndValues(i);  }

}


void DogAsynch1d::SetBndValues(int i)
// Set boundary conditions if i is either at left or right end
{
  // Left ghost cell
  if (i==mx)
    {
      // q
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    qold->set(0,m,k, qold->get(mx,m,k) );
	    q->set(0,m,k, q->get(mx,m,k) );

	    for (int s=1; s<=(kmax-1); s++)
	      {  Tmat->set(0,m,k,s, Tmat->get(mx,m,k,s) );  }
	  }
      
      // aux
      for (int m=1; m<=maux; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    auxold->set(0,m,k, auxold->get(mx,m,k) );
	    aux->set(0,m,k, aux->get(mx,m,k) );
	  }    
      
      // time
      time_old->set(0, time_old->get(mx) );
      time->set(0, time->get(mx) );   

      // dx
      dx->set(0, dx->get(mx) );

      // ftally
      ftally_lft->set(1, ftally_lft->get(mx+1) );

      // flux 
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    Flft->set(1,m,k, Flft->get(mx+1,m,k) );
	  }
    }
  // Right ghost cell    
  else if (i==1)
    {
      // q
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {         
	    qold->set(mx+1,m,k, qold->get(1,m,k) );
	    q->set(mx+1,m,k, q->get(1,m,k) );
	    
	    for (int s=1; s<=(kmax-1); s++)
	      {  Tmat->set(mx+1,m,k,s, Tmat->get(1,m,k,s) );  }
	  }
      
      // aux
      for (int m=1; m<=maux; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    auxold->set(mx+1,m,k, auxold->get(1,m,k) );
	    aux->set(mx+1,m,k, aux->get(1,m,k) );
	  }
      
      // time
      time_old->set(mx+1, time_old->get(1) );
      time->set(mx+1, time->get(1) );

      // dx
      dx->set(mx+1, dx->get(1) );

      // ftally
      ftally_rgt->set(mx+1, ftally_rgt->get(1) );

      // flux 
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    Frgt->set(mx+1,m,k, Frgt->get(1,m,k) );
	  }
    }

}
  
void DogAsynch1d::SetQstates(int i, dTensor2& Ql, dTensor2& Qr, string LR)
// Set left and right states at interface, interpolate in time if necessary
{
  if (LR=="left ")
    {
      double told = time_old->get(i-1);
      double tnew = time->get(i-1);
      double a;
      double dt1 = tnew-told;
      double dt2 = dt1*dt1;
      double dt3 = dt2*dt1;
      double dt4 = dt3*dt1;
      if (fabs(dt1)<=1.0e-15)
	{ a = 0.0; }
      else
	{ a = (time->get(i)-told)/dt1; }
      double a2 = a*a;
      double a3 = a2*a;
      double a4 = a3*a;
      double a5 = a4*a;

      switch(kmax)
	{
	case 1:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold = qold->get(i-1,m,k);
		double Qnew = q->get(i-1,m,k);

		Ql.set(m,k, (1.0-a)*Qold + a*Qnew );
		Qr.set(m,k, q->get(i,m,k) );
	      }
	  break;
	  
	case 2:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		//		double Qold = qold->get(i-1,m,k);
		//		double Qnew = q->get(i-1,m,k);
		//		double q_t  = Tmat->get(i-1,m,k,1);
		
		//		Ql.set(m,k, (1.0-a2)*Qold + a2*Qnew + a*dt1*(1.0-a)*q_t );
		//		Qr.set(m,k, q->get(i,m,k) );

		double Qold = qold->get(i-1,m,k);
		double Qnew = q->get(i-1,m,k);

		Ql.set(m,k, (1.0-a)*Qold + a*Qnew );
		Qr.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 3:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold = qold->get(i-1,m,k);
		double Qnew = q->get(i-1,m,k);
		double q_t  = Tmat->get(i-1,m,k,1);
		double q_tt = Tmat->get(i-1,m,k,2);
		
		Ql.set(m,k, (1.0-a3)*Qold + a3*Qnew +  dt1*a*(1.0-a2)*q_t 
		       + 0.5*a2*dt2*(1.0-a)*q_tt );
		Qr.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 4:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold  = qold->get(i-1,m,k);
		double Qnew  = q->get(i-1,m,k);
		double q_t   = Tmat->get(i-1,m,k,1);
		double q_tt  = Tmat->get(i-1,m,k,2);
		double q_ttt = Tmat->get(i-1,m,k,3);
		
		Ql.set(m,k, (1.0-a4)*Qold + a4*Qnew + a*dt1*(1.0-a3)*q_t 
		       + 0.5*a2*dt2*(1.0-a2)*q_tt + a3*dt3/6.0*(1.0-a)*q_ttt );
		Qr.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 5:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold   = qold->get(i-1,m,k);
		double Qnew   = q->get(i-1,m,k);
		double q_t    = Tmat->get(i-1,m,k,1);
		double q_tt   = Tmat->get(i-1,m,k,2);
		double q_ttt  = Tmat->get(i-1,m,k,3);
		double q_tttt = Tmat->get(i-1,m,k,4);
		
		Ql.set(m,k, (1.0-a5)*Qold + a5*Qnew + a*dt1*(1.0-a4)*q_t 
		       + 0.5*a2*dt2*(1.0-a3)*q_tt + a3*dt3/6.0*(1.0-a2)*q_ttt 
		       + a4*dt4/24.0*(1.0-a)*q_tttt );
		Qr.set(m,k, q->get(i,m,k) );
	      }	  
	  break;

	default:
	  cout << " Error in DogAsynch1d::EdgeFlux, kmax = " << kmax << " is not supported." << endl;
	  exit(1);
	  break;
	}
    }
  else if(LR=="right")
    {       
      double told = time_old->get(i+1);
      double tnew = time->get(i+1);
      double a;
      double dt1 = tnew-told;
      double dt2 = dt1*dt1;
      double dt3 = dt2*dt1;
      double dt4 = dt3*dt1;
      if (fabs(dt1)<=1.0e-15)
	{ a = 0.0; }
      else
	{ a = (time->get(i)-told)/dt1; }
      double a2 = a*a;
      double a3 = a2*a;
      double a4 = a3*a;
      double a5 = a4*a;

      switch(kmax)
	{
	case 1:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold = qold->get(i+1,m,k);
		double Qnew = q->get(i+1,m,k);
		
		Qr.set(m,k, (1.0-a)*Qold + a*Qnew );
		Ql.set(m,k, q->get(i,m,k) );
	      }
	  break;
	  
	case 2:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold = qold->get(i+1,m,k);
		double Qnew = q->get(i+1,m,k);
		double q_t  = Tmat->get(i+1,m,k,1);
		
		Qr.set(m,k, (1.0-a2)*Qold + a2*Qnew + a*dt1*(1.0-a)*q_t );
		Ql.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 3:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold = qold->get(i+1,m,k);
		double Qnew = q->get(i+1,m,k);
		double q_t  = Tmat->get(i+1,m,k,1);
		double q_tt = Tmat->get(i+1,m,k,2);
		
		Qr.set(m,k, (1.0-a3)*Qold + a3*Qnew +  dt1*a*(1.0-a2)*q_t 
		       + 0.5*a2*dt2*(1.0-a)*q_tt );
		Ql.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 4:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold  = qold->get(i+1,m,k);
		double Qnew  = q->get(i+1,m,k);
		double q_t   = Tmat->get(i+1,m,k,1);
		double q_tt  = Tmat->get(i+1,m,k,2);
		double q_ttt = Tmat->get(i+1,m,k,3);
		
		Qr.set(m,k, (1.0-a4)*Qold + a4*Qnew + a*dt1*(1.0-a3)*q_t 
		       + 0.5*a2*dt2*(1.0-a2)*q_tt + a3*dt3/6.0*(1.0-a)*q_ttt );
		Ql.set(m,k, q->get(i,m,k) );
	      }
	  break;

	case 5:
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		double Qold   = qold->get(i+1,m,k);
		double Qnew   = q->get(i+1,m,k);
		double q_t    = Tmat->get(i+1,m,k,1);
		double q_tt   = Tmat->get(i+1,m,k,2);
		double q_ttt  = Tmat->get(i+1,m,k,3);
		double q_tttt = Tmat->get(i+1,m,k,4);
		
		Qr.set(m,k, (1.0-a5)*Qold + a5*Qnew + a*dt1*(1.0-a4)*q_t 
		       + 0.5*a2*dt2*(1.0-a3)*q_tt + a3*dt3/6.0*(1.0-a2)*q_ttt 
		       + a4*dt4/24.0*(1.0-a)*q_tttt );
		Ql.set(m,k, q->get(i,m,k) );
	      }	  
	  break;

	default:
	  cout << " Error in DogAsynch1d::EdgeFlux, kmax = " << kmax << " is not supported." << endl;
	  exit(1);
	  break;
	}
    }
}

// Compute spatial derivatives of Ql and Qr
void DogAsynch1d::ComputeDerivatives(dTensor2& Ql,dTensor2& Qr,dTensor2& DQl,dTensor2& DQr)
{
  switch(kmax)
    {
    case 1:  // first-order in space

      for (int m=1; m<=meqn; m++)
	{
	  DQr.set(m,1, Qr.get(m,1) );
	  
	  DQl.set(m,1, Ql.get(m,1) );
	}

      break;
      
    case 2:  // second-order in space

      for (int m=1; m<=meqn; m++)
	{
	  DQr.set(m,1, Qr.get(m,1) - sq3*Qr.get(m,2) );
	  DQr.set(m,2, sq3*Qr.get(m,2) );
	  
	  DQl.set(m,1, Ql.get(m,1) + sq3*Ql.get(m,2) );
	  DQl.set(m,2, sq3*Ql.get(m,2) );
	}

      break;

    case 3:  // third-order in space
      
      for (int m=1; m<=meqn; m++)
	{
	  DQr.set(m,1, Qr.get(m,1) - sq3*Qr.get(m,2) + sq5*Qr.get(m,3) );
	  DQr.set(m,2, sq3*Qr.get(m,2) - 3.0*sq5*Qr.get(m,3) );
	  DQr.set(m,3, 3.0*sq5*Qr.get(m,3) );
	  
	  DQl.set(m,1, Ql.get(m,1) + sq3*Ql.get(m,2) + sq5*Ql.get(m,3) );
	  DQl.set(m,2, sq3*Ql.get(m,2) + 3.0*sq5*Ql.get(m,3) );
	  DQl.set(m,3, 3.0*sq5*Ql.get(m,3) );
	}

      break;

    case 4:  // fourth-order in space

      for (int m=1; m<=meqn; m++)
	{
	  DQr.set(m,1, Qr.get(m,1) - sq3*Qr.get(m,2) + sq5*Qr.get(m,3) - sq7*Qr.get(m,4) );
	  DQr.set(m,2, sq3*Qr.get(m,2) - 3.0*sq5*Qr.get(m,3) + 6.0*sq7*Qr.get(m,4) );
	  DQr.set(m,3, 3.0*sq5*Qr.get(m,3) - 15.0*sq7*Qr.get(m,4) );
	  DQr.set(m,4, 15.0*sq7*Qr.get(m,4) );	  
	  
	  DQl.set(m,1, Ql.get(m,1) + sq3*Ql.get(m,2) + sq5*Ql.get(m,3) + sq7*Ql.get(m,4) );
	  DQl.set(m,2, sq3*Ql.get(m,2) + 3.0*sq5*Ql.get(m,3) + 6.0*sq7*Ql.get(m,4) );
	  DQl.set(m,3, 3.0*sq5*Ql.get(m,3) + 15.0*sq7*Ql.get(m,4) );
	  DQl.set(m,4, 15.0*sq7*Ql.get(m,4) );
	}

      break;
  
    case 5:  // fifth-order in space

      for (int m=1; m<=meqn; m++)
	{
	  DQr.set(m,1, Qr.get(m,1) - sq3*Qr.get(m,2) + sq5*Qr.get(m,3) - sq7*Qr.get(m,4) + 3.0*Qr.get(m,5) );
	  DQr.set(m,2, sq3*Qr.get(m,2) - 3.0*sq5*Qr.get(m,3) + 6.0*sq7*Qr.get(m,4) - 30.0*Qr.get(m,5) );
	  DQr.set(m,3, 3.0*sq5*Qr.get(m,3) - 15.0*sq7*Qr.get(m,4) + 135.0*Qr.get(m,5) );
	  DQr.set(m,4, 15.0*sq7*Qr.get(m,4) - 315.0*Qr.get(m,5) );
	  DQr.set(m,5, 315.0*Qr.get(m,5) );
	  
	  DQl.set(m,1, Ql.get(m,1) + sq3*Ql.get(m,2) + sq5*Ql.get(m,3) + sq7*Ql.get(m,4) + 3.0*Ql.get(m,5) );
	  DQl.set(m,2, sq3*Ql.get(m,2) + 3.0*sq5*Ql.get(m,3) + 6.0*sq7*Ql.get(m,4) + 30.0*Ql.get(m,5) );
	  DQl.set(m,3, 3.0*sq5*Ql.get(m,3) + 15.0*sq7*Ql.get(m,4) + 135.0*Ql.get(m,5) );
	  DQl.set(m,4, 15.0*sq7*Ql.get(m,4) + 315.0*Ql.get(m,5) );
	  DQl.set(m,5, 315.0*Ql.get(m,5) );
	}

      break;

    default:
      cout << " ERROR, kmax = " << kmax << " is not yet implemented." << endl;
      exit(1);
    }
}


// Solve Riemann problem
void DogAsynch1d::RiemannSolve(int i,double u,dTensor2& DQl,dTensor2& DQr,string LR)
{
  if (LR=="left ")
    {
      if (u>=0.0)
	{
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {
		Ltmp->set(m,k, facm->get(i-1,k)*DQl.get(m,k) );  
	      }	    
	}
      else
	{
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      { 
		Ltmp->set(m,k, facm->get(i,k)*DQr.get(m,k) );  
	      }
	}
    }
  else if (LR=="right")
    {
      if (u>=0.0)
	{
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {	
		Rtmp->set(m,k, facm->get(i,k)*DQl.get(m,k) );  
	      }  
	}
      else
	{
	  for (int m=1; m<=meqn; m++)
	    for (int k=1; k<=kmax; k++)
	      {  
		Rtmp->set(m,k, facm->get(i+1,k)*DQr.get(m,k) );  
	      }
	}	
    }
}
