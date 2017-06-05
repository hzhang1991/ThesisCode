#ifndef _RUNDOGPACK_UNST_H_
#define _RUNDOGPACK_UNST_H_

// Necessary header files for RunDogpack_Unst:
#include "dogdefs.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "DogStateUnst2.h"
#include "DogParamsUnst2.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "ext_time.h" /* for get_utime and timeval_diff */

// ------------------------------------------------------------
// Function definitions
void RunMeshCopyScript(const string& outputdir);

//void ParseArguments(int argc,char**argv,string& outputdir);
//void GetParams_Unst(string inputdir,string& time_stepping_method);
void L2Project_Unst(
    const dTensor2* vel_vec,
    const int istart, 
        const int iend, 
        const int QuadOrder, 
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, 
        const dTensor3* auxin, 
        dTensor3* fout, 
    void (*Func)(
        const dTensor2* vel_vec,
        const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor2&));

void QinitFuncWrapper (const dTensor2* vel_vec, const dTensor2&,const dTensor2&,const dTensor2&,dTensor2&);
void AuxFuncWrapper   (const dTensor2* vel_vec, const dTensor2&,const dTensor2&,const dTensor2&,dTensor2&);

void AfterQinit_Unst(const mesh& Mesh, dTensor3& aux, dTensor3& q);
void Output_Unst(const mesh& Mesh, const dTensor3& aux,
        const dTensor3& q, double t, int nframe, 
        string outputdir);

void SetEdgeData_Unst(int NumEdges, int morder, int kmax,
        const mesh& Mesh, edge_data_Unst& EdgeData);
void SetEdgeData_Unst(const mesh& Mesh, 
        int NumQuadPoints, 
        int NumBasisOrder, 
        edge_data_Unst& EdgeData);
void DogSolveRK_Unst(
    const dTensor2* vel_vec,
    const mesh& Mesh,const edge_data_Unst& EdgeData,
        dTensor3& aux,dTensor3& qold, dTensor3& qnew,		       
        const double tstart,const double tend, DogStateUnst2& dogStateUnst2 );
void DogSolveUser_Unst(const mesh& Mesh,const edge_data_Unst& EdgeData,
        dTensor3& aux,dTensor3& qold,dTensor3& qnew,
        const double tstart,const double tend, 
        const string outputdir);
// ------------------------------------------------------------

#endif
