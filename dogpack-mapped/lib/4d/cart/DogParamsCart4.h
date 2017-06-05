#ifndef _DOGPARAMSCART4_H_
#define _DOGPARAMSCART4_H_

// Including this header file allows access to all the 3D paramters described in
// the [grid] section of any 3D parameters.ini file.
//
// The singleton dogParamsCart4 allows access to these parameters.  For example,
// after including DogParamsCart4.h, one may accesses mx, the number of grid cells
// via:
//
//      dogParamsCart4.get_mx()
//

#include <cstring>
#include "IniDocument.h"
#include "DogParams.h"
#include "assert.h"
#include "debug.h"
#include "dog_io.h"
#include "dog_str.h"

class IniDocument;

struct DogParamsCart4
{
    // IO routines.  (for printing to qhelp, etc.)
    public:

        // Accessors for plotting routines:
        void write_plotting_help_files( int plot_mx, int plot_my, int plot_mz, int plot_mw, const char* outputdir);
        void append_qhelp(const char* filename, int plot_mx, int plot_my, int plot_mz, int plot_mw);
        void append_output_parameters_ini(const char* filename, int plot_mx, int plot_my, int plot_mz, int plot_mw);
        bool get_is_initialized(){return is_initialized;}

        // Prints a list of the parameters to standard output:
        void reportParameters();

        // Constructor
        DogParamsCart4(): plot_rx(NULL), plot_ry(NULL), plot_rz(NULL), is_initialized(false){nframe=-1;}

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
        const int*    get_plot_rz() const{ return plot_rz;}

        int get_plot_mx(int idx) const;
        int get_plot_my(int idx) const;
        int get_plot_mz(int idx) const;
        int get_plot_mw(int idx) const;
        int    get_mx()    const{ return  mx;   }
        int    get_my()    const{ return  my;   }
        int    get_mz()    const{ return  mz;   }
        int    get_mw()    const{ return  mw;   }
        int    get_mbc()   const{ return  mbc;  }
        double get_xlow()  const{ return  xlow; }
        double get_xhigh() const{ return  xhigh;}
        double get_ylow()  const{ return  ylow; }
        double get_yhigh() const{ return  yhigh;}
        double get_zlow()  const{ return  zlow; }
        double get_zhigh() const{ return  zhigh;}
        double get_wlow()  const{ return  wlow; }
        double get_whigh() const{ return  whigh;}

        double get_dx()    const{ return  dx;   }
        double get_dy()    const{ return  dy;   }
        double get_dz()    const{ return  dz;   }
        double get_dw()    const{ return  dw;   }
        double get_prim_vol() const{ return prim_vol;}

        double get_xc(int i) const{ return xlow + (double(i)-0.5)*dx;}
        double get_yc(int j) const{ return ylow + (double(j)-0.5)*dy;}
        double get_zc(int k) const{ return zlow + (double(k)-0.5)*dz;}
        double get_wc(int k) const{ return wlow + (double(k)-0.5)*dw;}
        double get_xl(int i) const{ return xlow + (double(i-1))*dx;}
        double get_yl(int j) const{ return ylow + (double(j-1))*dy;}
        double get_zl(int k) const{ return zlow + (double(k-1))*dz;}
        double get_wl(int k) const{ return wlow + (double(k-1))*dw;}
        double get_xlength() const{ return xhigh - xlow;}
        double get_ylength() const{ return yhigh - ylow;}
        double get_zlength() const{ return zhigh - zlow;}
        double get_wlength() const{ return whigh - wlow;}

        //
        // Methods for changing user supplied parameters:
        //
    public:
        void reset_mx(int mx_in);
        void reset_my(int my_in);
        void reset_mz(int mz_in);
        void reset_mw(int mw_in);
        void set_xlims(double,double);
        void set_ylims(double,double);
        void set_zlims(double,double);
        void set_wlims(double,double);

        int get_nframe() const{ return nframe; }
        void set_nframe(const int in) { nframe=in; }

        //
        // Variables:
        //
    private:

        bool is_initialized;
        int mx,my,mz,mw;
        int mbc;
        double xlow,xhigh;
        double ylow,yhigh;
        double zlow,zhigh;
        double wlow,whigh;
        int nframe;

        // calculated (derived) quantities
        double dx,dy,dz,dw;
        double prim_vol;

        //
        int* plot_rx;
        int* plot_ry;
        int* plot_rz;
        int* plot_rw;

        //
        // Methods:
        //
    private:
        void update_dx();
        void update_dy();
        void update_dz();
        void update_dw();
        void checkParameters();
        void setDerivedParameters();
};

extern DogParamsCart4 dogParamsCart4;

#endif
