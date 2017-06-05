#ifndef _RUN_DOGPACK_H_
#define _RUN_DOGPACK_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "IniDocument.h"
#include "DogParamsCart1.h"

// If we want to use DogSolver from the top-level library, this needs to be
// written:
// #include "DogState1d.h"  

using namespace std;


// ------------------------------------------------------------
// Functions use in RunDogpack.cpp
void GridSetup(int,double,double,dTensor2&,dTensor1&);
void L2Project(int mopt, int istart, int iend,
        const dTensor2& node,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
        dTensorBC3& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, 
            const dTensor2&, dTensor2&));
void Output(const dTensor2&,const dTensorBC3&,const dTensorBC3&,
        double,int,string);
void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
void ConSoln(const int method[], const dTensor2& node, const dTensorBC3& aux,
        const dTensorBC3& q, double t, string outputdir);
void DogSolveRK(const dTensor2&,const dTensor1&,
        dTensorBC3&,dTensorBC3&,dTensorBC3&,
        dTensorBC1&,double,double,int,const int[],
        double[],const double[],string);
void DogSolveSDC(const dTensor2&,const dTensor1&,
        dTensorBC3&,dTensorBC3&,dTensorBC3&,
        dTensorBC1&,double,double,int,const int[],
        double[],const double[],string);
void DogSolveLxW(const dTensor2&,const dTensor1&,
        dTensorBC3&,dTensorBC3&,dTensorBC3&,
        dTensorBC1&,double,double,int,const int[],
        double[],const double[],string);
void DogSolveUser(const dTensor2&,const dTensor1&,
        dTensorBC3&,dTensorBC3&,dTensorBC3&,
        dTensorBC1&,double,double,int,const int[],
        double[],const double[],string);
void InitApp(IniDocument& ini_doc);
// ------------------------------------------------------------

#endif
