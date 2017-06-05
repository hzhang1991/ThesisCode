#include "tensors.h"
#include "constants.h"

// declare inner methods in-line to avoid function call overhead
//
namespace OutSliceCart3Inline
{
  inline void slice_q_local_1(const int kmax2d,
			      const dTensor1& qin,
			      dTensor1& qout)
  {
    switch(kmax2d)
      {
      case 10: // kmax=20
	qout.set(10, qin.get(19) );
	qout.set(9,  qin.get(18) );
	qout.set(8,  qin.get(13) );
	qout.set(7,  qin.get(11) );
	qout.set(6,  qin.get(9)  );
	qout.set(5,  qin.get(8)  );
	qout.set(4,  qin.get(5)  );
	qout.set(3,  qin.get(3) - 0.5*sq5*qin.get(16) );
	qout.set(2,  qin.get(2) - 0.5*sq5*qin.get(15) );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(10) );
	break;
	
      case 6: // kmax = 10
	qout.set(6,  qin.get(9)  );
	qout.set(5,  qin.get(8)  );
	qout.set(4,  qin.get(5)  );
	qout.set(3,  qin.get(3)  );
	qout.set(2,  qin.get(2)  );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(10) );
	break;
	
      case 3: // kmax = 4
	qout.set(3,  qin.get(3)  );
	qout.set(2,  qin.get(2)  );
	qout.set(1,  qin.get(1)  );
	break;
	
      case 1: // kmax = 1
	qout.set(1,  qin.get(1)  );
	break;
	
      default:
	// should never get here
	eprintf("Error in slice_q_local, kmax2d must be 1,3,6,10.  kmax = %d\n",kmax2d);
      }
  }
  
  inline void slice_q_local_2(const int kmax2d,
			      const dTensor1& qin,
			      dTensor1& qout)
  {
    switch(kmax2d)
      {
      case 10: // kmax=20
	qout.set(10, qin.get(20) );
	qout.set(9,  qin.get(18) );
	qout.set(8,  qin.get(15) );
	qout.set(7,  qin.get(12) );
	qout.set(6,  qin.get(10) );
	qout.set(5,  qin.get(8)  );
	qout.set(4,  qin.get(6)  );
	qout.set(3,  qin.get(4) - 0.5*sq5*qin.get(14) );
	qout.set(2,  qin.get(2) - 0.5*sq5*qin.get(13) );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(9)  );
	break;
	
      case 6: // kmax = 10
	qout.set(6,  qin.get(10) );
	qout.set(5,  qin.get(8)  );
	qout.set(4,  qin.get(6)  );
	qout.set(3,  qin.get(4)  );
	qout.set(2,  qin.get(2)  );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(9)  );
	break;
	
      case 3: // kmax = 4
	qout.set(3,  qin.get(4)  );
	qout.set(2,  qin.get(2)  );
	qout.set(1,  qin.get(1)  );
	break;
	
      case 1: // kmax = 1
	qout.set(1,  qin.get(1)  );
	break;
	
      default:
	// should never get here
	eprintf("Error in slice_q_local, kmax2d must be 1,3,6,10.  kmax = %d\n",kmax2d);
      }
  }
  
  inline void slice_q_local_3(const int kmax2d,
			      const dTensor1& qin,
			      dTensor1& qout)
  {
    switch(kmax2d)
      {
      case 10: // kmax=20
	qout.set(10, qin.get(20) );
	qout.set(9,  qin.get(19) );
	qout.set(8,  qin.get(16) );
	qout.set(7,  qin.get(14) );
	qout.set(6,  qin.get(10) );
	qout.set(5,  qin.get(9)  );
	qout.set(4,  qin.get(7)  );
	qout.set(3,  qin.get(4) - 0.5*sq5*qin.get(12) );
	qout.set(2,  qin.get(3) - 0.5*sq5*qin.get(11) );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(8)  );
	break;
	
      case 6: // kmax = 10
	qout.set(6,  qin.get(10) );
	qout.set(5,  qin.get(9)  );
	qout.set(4,  qin.get(7)  );
	qout.set(3,  qin.get(4)  );
	qout.set(2,  qin.get(3)  );
	qout.set(1,  qin.get(1) - 0.5*sq5*qin.get(8)  );
	break;
	
      case 3: // kmax = 4
	qout.set(3,  qin.get(4)  );
	qout.set(2,  qin.get(3)  );
	qout.set(1,  qin.get(1)  );
	break;
	
      case 1: // kmax = 1
	qout.set(1,  qin.get(1)  );
	break;
	
      default:
	// should never get here
	eprintf("Error in slice_q_local, kmax2d must be 1,3,6,10.  kmax = %d\n",kmax2d);
      }
  }
  
}
