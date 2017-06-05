#ifndef _OUTSLICECart4_H_
#define _OUTSLICECart4_H_

// Including this header file allows access to all the 3D paramters described in
// the [outslice] section of any 3D parameters.ini file.
//
// The singleton outSliceCart4 allows access to these parameters.  For example,
// after including OutSliceCart4.h, one may accesses numslices, the number of slices
// via:
//
//      outSliceCart4.get_numslices()
//

#include <cstring>
#include "IniDocument.h"
#include "DogParams.h"
#include "DogParamsCart4.h"
#include "assert.h"
#include "debug.h"
#include "dog_io.h"
#include "dog_str.h"
#include "tensors.h"
#include "OutSliceCart4Inline.h"

class IniDocument;

struct OutSliceCart4
{
    public:

        // Constructor
        OutSliceCart4();

        // Destructor
        ~OutSliceCart4();

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
        int get_output4D() const;

        // Output function
        void write_output(const int noutput,
                const double t,
                const dTensorBC6& q,
                const dTensorBC6& aux);

        void slice_q(const int slicenumber,
                const dTensorBC6* q, dTensor4*& qslice);

        void OutSliceExtra(const int noutput,
                const double t);

    private:
        bool is_initialized;
        int numslices;
        int* slicekind;
        int* sliceidx;
        int output4D;
};

extern OutSliceCart4 outSliceCart4;

#endif
