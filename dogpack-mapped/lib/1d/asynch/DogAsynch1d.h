#ifndef DOGASYNCH1D_H
#define DOGASYNCH1D_H

#include "tensors.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
using namespace std;

// 1D DogAsynch1d Object ----------------------------------------
class DogAsynch1d
{
 public:
    // Constructor
    // POST: Creates a mesh
    DogAsynch1d(int inmx, int inmeqn, int inmaux, int inkmax, int inmbc);
    
    // Copy constructor
    // POST: New mesh created with size and contents same as another_mesh
    DogAsynch1d(const DogAsynch1d& another_mesh);
    
    // Destructor
    // POST: mesh no longer exists
    ~DogAsynch1d();
    
    // Ouput all mesh, q, and aux information
    void Output(string outputdir,double t,int nframe);

    // Ouput q and aux information
    void OutputSoln(string outputdir,double t,int nframe);

    // Ouput mesh information
    void OutputMesh(string outputdir);

    // Create uniform mesh
    void InputMesh(double xlow, double xhigh);    
    
    // Input all mesh information from file
    void InputMesh(string inputdir);

    // Returns "CFL"
    const double& get_CFL() const;

    // Sets "CFL"
    void set_CFL(double in_CFL);

    // Returns "CFL_MAX"
    const double& get_CFL_MAX() const;

    // Sets "CFL_MX"
    void set_CFL_MAX(double in_CFL_MAX);

    // Returns "mx"
    const int& get_mx() const;

    // Returns "meqn"
    const int& get_meqn() const;

    // Returns "maux"
    const int& get_maux() const;

    // Returns "kmax"
    const int& get_kmax() const;
    
    // Returns "mbc"
    const int& get_mbc() const;

    // Returns "tstart"
    const double& get_tstart() const;
    
    // Returns "tend"
    const double& get_tend() const;   

    // Returns "num_to_be_updated"
    const int& get_num_to_be_updated() const;
    
    // Set the correct values of dx once node values have been set
    void dx_compute();    

    // ------------------------------------
    // Asynchronous time-stepping functions
    // ------------------------------------

    // Returns 1:  if element i is at time <= its neighbors
    // Returns 0:  otherwise
    bool is_updatable(int i) const;

    // Initialize state to time advance from tstart to tend
    void initialize_for_time_stepping(double intstart, double intend);
    
    // Find the next element that should be updated
    int get_next_updatable_element();

    // Check to see if input element is at final time
    // If   NO: do nothing
    // If  YES: remove input element from "elements_to_be_updated" list
    void check_if_done(int i);

    // Set boundary conditions if i is either at left or right end
    void SetBndValues(int i);

    // Lax-Wendroff method for updating element
    void LxW_update(int i);

    // Compute edge fluxes for element i
    void EdgeFlux(int i, double u, string LR);

    // Ql and Qr
    void SetQstates(int i, dTensor2& Ql, dTensor2& Qr, string LR);
    
    // Compute spatial derivatives of Ql and Qr
    void ComputeDerivatives(dTensor2& Ql,dTensor2& Qr,dTensor2& DQl,dTensor2& DQr);

    // Solve Riemann problem
    void RiemannSolve(int i,double u,dTensor2& DQl,dTensor2& Dqr,string LR);
    
    // Construct time derivatives and store in T
    void ConstructMRL(int i, double u);

    // Compute internal part of update
    void SetMmat(int i,double u,dTensor2& Qin);

    // Compute left part of udpate
    void SetLmat(int i);
    
    // Compute right part of update
    void SetRmat(int i);
    
    // Update solution
    void UpdateSoln(int i);
    
    // Conservation fix
    void ConservationFix(int i, double u);

    // Correct flux at interface
    void FluxCorrect(int i, double u, string LR);

    // Compute the time step, dt, that will
    //       give the slowest moving element a courant
    //       number of CFL.
    double get_max_dt_step();

    // Compute the time step, dt, that will
    //       give the fastest moving element a courant
    //       number of CFL.
    double get_min_dt_step();

    // -------------------
    // Important variables
    // -------------------

    // q: solution on the mesh
    dTensorBC3* q;      // (mx,meqn,kmax,mbc);

    // qold: previous time step solution on the mesh
    dTensorBC3* qold;      // (mx,meqn,kmax,mbc);

    // aux: auxiliary array on the mesh
    dTensorBC3* aux;    // (mx,maux,kmax,mbc);

    // auxold: auxiliary array on the mesh
    dTensorBC3* auxold;    // (mx,maux,kmax,mbc);
    
    // node: list of x coordinates of nodes
    dTensorBC1* node;   // (mx+1);

    // element lengths
    dTensorBC1* dx;       // (mx)

    // element times
    dTensorBC1* time;     // (mx)
    dTensorBC1* time_old; // (mx)
    
    // element time-steps
    dTensorBC1* dt;       // (mx)

    // Correction matrices
    dTensor2* Mmat;     // (meqn,kmax)
    dTensor2* Rmat;     // (meqn,kmax)
    dTensor2* Lmat;     // (meqn,kmax)

    // Time derivative matrix
    dTensorBC4* Tmat;     // (mx,meqn,1..kmax,1..kmax-1)

    // Stored fluxes at interfaces --> this is need for conservation fix-ups
    dTensor3* Flft;       // (mx+1,meqn,kmax)
    dTensor3* Frgt;       // (mx+1,meqn,kmax)
    iTensor1* ftally_lft; // (mx+1)
    iTensor1* ftally_rgt; // (mx+1)

 private:
    // Basic constants
    double CFL;     // CFL number
    double CFL_MAX; // max allowed CFL number
    int mx;      // Number of total elements in mesh
    int meqn;    // Number of conservation laws
    int maux;    // Number of aux arrays
    int kmax;    // Number of Legendre coefficients
    int mbc;     // Number of ghost cells on each end
    double tstart,tend,dt_max;  // Starting and ending times for time-stepping
    int num_to_be_updated;      // Number of elements that still need to be updated to tend
    int curr_to_be_updated;     // Pointer to current element that still needs to be updated
    iTensor1* elements_to_be_updated;   // List that keeps track of elements that
                                        // have not yet reached the final time.

    // vector of (2*u/dx)^(k-1) for each element i
    dTensorBC2* facp;   // (mx,kmax)

    // vector of (-2*u/dx)^(k-1) for each element i
    dTensorBC2* facm;   // (mx,kmax)

    // Helpers for correction matrices
    dTensor3* Mmat_s;   // (meqn,kmax,kmax-1);
    dTensor2* Ltmp;     // (meqn,kmax)
    dTensor2* Rtmp;     // (meqn,kmax)

    // matrix of phi derivatives at xi=-1 (left) and xi=1 (right)
    dTensor2* dphiR;    // (l=5,s+1=1..5)
    dTensor2* dphiL;    // (l=5,s+1=1..5)

    // vector of u^s
    double* us;

    // vector of dt^s
    double* dts;
};
// ---------------------------------------------------------

#endif
