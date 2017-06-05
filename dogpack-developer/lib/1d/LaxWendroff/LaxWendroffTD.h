#ifndef _LaxWendroffTD_H_
#define _LaxWendroffTD_H_

    void EvalRiemmanData(const int& i, const dTensorBC3& q, const dTensorBC3& aux,  
            const dTensorBC3& LxWF, dTensor1& Ql, dTensor1& Qr, 
            dTensor1& LxFFluxl, dTensor1& LxFFluxr,
            dTensor1& Auxl, dTensor1& Auxr);

    double RiemannSolveLxW(const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
            const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
            dTensor1& Fl, dTensor1& Fr,
            void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
                const dTensor1&,double&,double&));

    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);

    void L2Project(int,int,int,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&,dTensorBC3&,
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

    void L2ProjectLxW(const int istart, const int iend, 
            const dTensor2& node, const dTensorBC3& qin, const dTensorBC3& auxin,  
            dTensorBC3& IntF1, dTensorBC3& IntF2, dTensorBC3& IntF3,
            dTensorBC3& Flux, dTensorBC3& LxWFlux2, dTensorBC3& LxWFlux3,
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
            void (*DFunc)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor3&),
            void (*D2Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor4&));

    void L2ProjectLxWTD(const int istart, const int iend, 
            const dTensor2& node, const dTensorBC3& qin, const dTensorBC3& auxin,  
            dTensorBC3& IntF1, dTensorBC3& IntF2, 
            dTensorBC3& Flux, dTensorBC3& LxWFlux2, 
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
            void (*DFunc)(const dTensor1&, const dTensor2&, const dTensor2&,
            dTensor3&) );

    // Flux function and its derivatives:
    void FluxFunc   (const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void DFluxFunc  (const dTensor1&,const dTensor2&,const dTensor2&,dTensor3&);
    void D2FluxFunc (const dTensor1&,const dTensor2&,const dTensor2&,dTensor4&);

    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void LstarExtra(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);
    void CopyQ(const dTensorBC3&,dTensorBC3&);

#endif
