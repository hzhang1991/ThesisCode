class dTensor1;
class dTensor2;
class dTensorBC4;

void ApplyLimiterKrivodonova(
  dTensorBC4& aux, dTensorBC4& q,
  void (*ProjectRightEig)(int,const dTensor1&,
    const dTensor1&,const dTensor2&, dTensor2&),
  void (*ProjectLeftEig)(int,const dTensor1&,
    const dTensor1&,const dTensor2&, dTensor2&));
