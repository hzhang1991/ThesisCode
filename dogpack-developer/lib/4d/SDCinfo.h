#ifndef _SDCINFO_H_
#define _SDCINFO_H_

// SDC information
class dTensor1;
struct SDCinfo
{
  int num_iter;
  int max_iter;
  int num_stage;
  int timevec[10];

  bool is_final_iter()
  {
    return (num_iter==max_iter);
  }
};

#endif
