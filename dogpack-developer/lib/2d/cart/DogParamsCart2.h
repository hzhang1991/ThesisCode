#ifndef _DOGPARAMSCART2_H_
#define _DOGPARAMSCART2_H_

// Including this header file allows access to all the 2D paramters described in
// the [grid] section of any 2D parameters.ini file.
//
// The singleton dogParamsCart2 allows access to these parameters.  For example,
// after including DogParamsCart2.h, one may accesses mx, the number of grid cells
// via:
//
//      dogParamsCart2.get_mx()
//
// See DogParamsCart2.cpp to see how all this data gets parsed.

#include <cstring>
#include "IniDocument.h"
#include "DogParams.h"
#include "assert.h"
#include "debug.h"
#include "dog_io.h"
#include "dog_str.h"

class IniDocument;

struct DogParamsCart2
{

    // IO routines.  (for printing to qhelp, etc.)
    public:

        // Accessors for plotting routines:
        void write_qhelp(const char* outputdir,int plot_mx,int plot_my); //deprecated
        void write_output_parameters_ini(const char* filename, int plot_mx, int plot_my);
        void write_plotting_help_files( int plot_mx, int plot_my, const char* outputdir);
        bool get_is_initialized(){return is_initialized;}

        // Prints a list of the parameters to standard output:
        void reportParameters();

        // Constructor
        DogParamsCart2(): plot_rx(NULL), plot_ry(NULL), is_initialized(false){nframe=-1;}

        // iniitialize (i.e. parse and save) all the parameters from [grid]
        // section of parameters.ini file.  This routine remembers if all these
        // parameters have been previously saved.
        void init();

    //
    // Accessors:
    //
    public:

        const int*    get_plot_rx() const{ return plot_rx;}
        const int*    get_plot_ry() const{ return plot_ry;}
        int get_plot_mx(int idx) const;
        int get_plot_my(int idx) const;
        int    get_mx()    const{ return  mx;   }
        int    get_my()    const{ return  my;   }
        int    get_mbc()   const{ return  mbc;  }
        double get_xlow()  const{ return  xlow; }
        double get_xhigh() const{ return  xhigh;}
        double get_ylow()  const{ return  ylow; }
        double get_yhigh() const{ return  yhigh;}

        double get_dx()    const{ return  dx;   }
        double get_dy()    const{ return  dy;   }
        double get_prim_vol() const{ return prim_vol;}

        double get_xc(int i) const{ return xlow + (double(i)-0.5)*dx;}
        double get_yc(int j) const{ return ylow + (double(j)-0.5)*dy;}
        double get_xl(int i) const{ return xlow + (double(i-1))*dx;}
        double get_yl(int j) const{ return ylow + (double(j-1))*dy;}
        double get_xlength() const{ return xhigh - xlow;}
        double get_ylength() const{ return yhigh - ylow;}

    //
    // Methods for changing user supplied parameters:
    //
    public:
        void reset_mx(int mx_in);
        void reset_my(int my_in);
        void set_xlims(double,double);
        void set_ylims(double,double);

        int get_nframe() const{ return nframe; }
        void set_nframe(const int in) { nframe=in; }

    //
    // Variables:
    //
    private:
        
        int mx;
        int my;
        int mbc;
        double xlow;
        double xhigh;
        double ylow;
        double yhigh;
        int nframe;

        // calculated (derived) quantities
        double dx;
        double dy;
        double prim_vol;

        //
        int* plot_rx;
        int* plot_ry;

        //
        bool is_initialized;

    //
    // Methods:
    //
    private:
        void update_dx();
        void update_dy();
        void checkParameters();
        void setDerivedParameters();
};

extern DogParamsCart2 dogParamsCart2;

#endif
