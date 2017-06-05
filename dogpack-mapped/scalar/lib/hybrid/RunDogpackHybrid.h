#ifndef _RUNDOGPACKHYBRID_H_
#define _RUNDOGPACKHYBRID_H_

// Necessary header files for RunDogpackHybrid:
#include "dogdefs.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "DogStateHybrid.h"
#include "DogParamsUnst2.h"
#include "DogParamsCart2.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "ext_time.h" /* for get_utime and timeval_diff */

// ------------------------------------------------------------
// Function definitions
void RunMeshCopyScript(const string& outputdir);

//void ParseArguments(int argc,char**argv,string& outputdir);
//void GetParams_Unst(string inputdir,string& time_stepping_method);
void L2Project_Unst(const int istart, 
        const int iend, 
        const int QuadOrder, 
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, 
        const dTensor3* auxin, 
        dTensor3* fout, 
        void (*Func)(const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&));

void QinitFuncWrapper  (const dTensor2* vel_vec,const dTensor2&,const dTensor2&,dTensor2&);
void AuxFuncWrapper    (const dTensor2* vel_vec,const dTensor2&,const dTensor2&,dTensor2&);

void AfterQinit_Unst(const mesh& Mesh, dTensorBC5& q);
void Output_Unst(const mesh& Mesh, const dTensor3& aux,
        const dTensor3& q, double t, int nframe, 
        string outputdir);
void Output_Hybrid(const mesh& Mesh, 
        const dTensorBC5& q, double t, int nframe, 
        string outputdir);
//  void ConSoln_Unst(const mesh& Mesh, 
//          const dTensor3& aux, const dTensor3& q, 
//          double t, string outputdir);

void SetEdgeData_Unst(int NumEdges, int morder, int kmax,
        const mesh& Mesh, edge_data_Unst& EdgeData);
void SetEdgeData_Unst(const mesh& Mesh, 
        int NumQuadPoints, 
        int NumBasisOrder, 
        edge_data_Unst& EdgeData);

// Driver for the 2D, unstructured solver:
void DogSolveHybrid(const mesh& Mesh,const edge_data_Unst& EdgeData,
        dTensorBC5& q, 
        const double tstart,const double tend );

// Driver for the 2D, structured, semi-lagrangian solver:
// TODO
// ------------------------------------------------------------

#endif
