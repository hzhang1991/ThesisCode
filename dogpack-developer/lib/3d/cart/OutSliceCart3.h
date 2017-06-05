#ifndef _OUTSLICECART3_H_
#define _OUTSLICECART3_H_

// Including this header file allows access to all the 3D paramters described in
// the [outslice] section of any 3D parameters.ini file.
//
// The singleton outSliceCart3 allows access to these parameters.  For example,
// after including OutSliceCart3.h, one may accesses numslices, the number of slices
// via:
//
//      outSliceCart3.get_numslices()
//

#include <cstring>
#include "IniDocument.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
#include "assert.h"
#include "debug.h"
#include "dog_io.h"
#include "dog_str.h"
#include "tensors.h"
#include "OutSliceCart3Inline.h"

class IniDocument;

struct OutSliceCart3
{
public:
  
  // Constructor
  OutSliceCart3();

  // Destructor
  ~OutSliceCart3();

  // initialize (i.e. parse and save) all the parameters from [outslice]
  // section of parameters.ini file.  This routine remembers if all these
  // parameters have been previously saved.
  void init();
  
  // Accessors for plotting routines:
  void write_qhelp_slice(const char* filename);
  void write_plotting_help_files(const char* outputdir);
  
  // Prints a list of the parameters to standard output:
  void reportParameters();
  void checkParameters();
  
  // Accessors:
  int get_numslices() const;
  int get_slicekind(int idx) const;
  int get_sliceidx(int idx) const;
  bool get_is_initialized() const;
  int get_output3D() const;

  // Output function
  void write_output(const int noutput,
		    const double t,
		    const dTensorBC5& q,
		    const dTensorBC5& aux);

  void slice_q(const int slicenumber,
	       const dTensorBC5* q,
	       dTensor4*& qslice);

  void OutSliceExtra(const int noutput,
		     const double t);

private:
  bool is_initialized;
  int numslices;
  int* slicekind;
  int* sliceidx;
  int output3D;
};

extern OutSliceCart3 outSliceCart3;

#endif
