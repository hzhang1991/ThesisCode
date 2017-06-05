#ifndef _APPLY_POS_LIMITER_H__ 
#define _APPLY_POS_LIMITER_H__

void SetPositivePoints  (const int& space_order, dTensor2& spts);
void SetLegendreAtPoints(const int& space_order, const dTensor2& spts, dTensor2& phi);

#endif
