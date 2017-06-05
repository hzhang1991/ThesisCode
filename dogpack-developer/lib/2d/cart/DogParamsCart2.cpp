#include "DogParamsCart2.h"

using std::string;
using std::vector;

// This file parses and provides options to access the [grid] section of the 
// parameters.ini file

DogParamsCart2 dogParamsCart2;

int DogParamsCart2::get_plot_mx(int idx) const 
{
    assert_divides(plot_rx[idx],mx);
    return mx/plot_rx[idx];
}

int DogParamsCart2::get_plot_my(int idx) const 
{
    assert_divides(plot_ry[idx],my);
    return my/plot_ry[idx];
}

void DogParamsCart2::update_dx()
{
  dx = (xhigh-xlow)/double(mx);
}
void DogParamsCart2::update_dy()
{
  dy = (yhigh-ylow)/double(my);
}

void DogParamsCart2::reset_mx(int mx_in)
{
    mx = mx_in;
    update_dx();
    //update_plot_mx();
    Wprintf("reset mx to %d\n", get_mx());
}
void DogParamsCart2::reset_my(int my_in)
{
    my = my_in;
    update_dy();
    //update_plot_my();
    Wprintf("reset my to %d\n", get_my());
}

void DogParamsCart2::set_xlims(double xlow_in, double xhigh_in)
{
    xlow = xlow_in;
    xhigh = xhigh_in;
    update_dx();
}

void DogParamsCart2::set_ylims(double ylow_in, double yhigh_in)
{
    ylow = ylow_in;
    yhigh = yhigh_in;
    update_dy();
}

static void not_set_error(string varname)
{
    eprintf("In parameters.ini, section [grid],"
            " you must set the parameter %s\n", varname.c_str());
}

static bool invalid_value(const char* varname, const char* val)
{
    eprintf("invalid value: %s = [%s]", varname, val);
    return false;
}

// It does not make sense to define defaults for
// these options, so we do not use the defaults
// mechanism used e.g. in DogParams::init
void DogParamsCart2::init()
{
    if(is_initialized) {
        return;
    }

    string section_label = "grid";
    IniDocument::Section& grid = ini_doc[section_label];

    // options for which we provide non-empty-string defaults
    vector<string> option_names_list;
    option_names_list.push_back("mx"   );
    option_names_list.push_back("my"   );
    option_names_list.push_back("mbc"  );
    option_names_list.push_back("xlow" );
    option_names_list.push_back("xhigh");
    option_names_list.push_back("ylow" );
    option_names_list.push_back("yhigh");

    // get defaults from e.g. $DOGPACK/config/dogParams_defaults.ini
    void get_defaults(
            IniDocument::Section& params,
            const string& section_label,
            const vector<string>& option_names_list);
    get_defaults      (grid, section_label, option_names_list);
    bool verify_options_set(
            const IniDocument::Section& params,
            const string& section_name,
            const vector<string>& option_names_list);
    verify_options_set(grid, section_label, option_names_list);

    const char* s_mx    = grid["mx"   ].c_str();
    const char* s_my    = grid["my"   ].c_str();
    const char* s_mbc   = grid["mbc"  ].c_str();
    const char* s_xlow  = grid["xlow" ].c_str();
    const char* s_xhigh = grid["xhigh"].c_str();
    const char* s_ylow  = grid["ylow" ].c_str();
    const char* s_yhigh = grid["yhigh"].c_str();

    // convert strings to numbers;
    // populate class with parameter data.
    //
    sscanf(s_mx   ,"%d" , &mx    ) || invalid_value("mx"    , s_mx    );
    sscanf(s_my   ,"%d" , &my    ) || invalid_value("my"    , s_my    );
    sscanf(s_mbc  ,"%d" , &mbc   ) || invalid_value("mbc"   , s_mbc   );
    sscanf(s_xlow ,"%lf", &xlow  ) || invalid_value("xlow"  , s_xlow  );
    sscanf(s_xhigh,"%lf", &xhigh ) || invalid_value("xhigh" , s_xhigh );
    sscanf(s_ylow ,"%lf", &ylow  ) || invalid_value("ylow"  , s_ylow  );
    sscanf(s_yhigh,"%lf", &yhigh ) || invalid_value("yhigh" , s_yhigh );

    // read plot mesh sizes out of parameters plot_[mr][xy]
    const char* s_plot_rx = grid["plot_rx"].c_str();
    const char* s_plot_ry = grid["plot_ry"].c_str();
    if(dogParams.get_how_many_plot_resolutions()!=0)
    {
        if(strlen(s_plot_rx)!=0)
        {
            int numel;
            numel = new_array_from_str(s_plot_rx,plot_rx,1);
            if(numel<0) eprintf("invalid syntax: plot_rx = %s", s_plot_rx);
            if(numel!=dogParams.get_how_many_plot_resolutions())
                eprintf("plot_rx=%s does not match up with nout_per_plot",
                        s_plot_rx);
            //
            numel = new_array_from_str(s_plot_ry,plot_ry,1);
            if(numel<0) eprintf("invalid syntax: plot_ry = %s", s_plot_ry);
            assert_printf(numel==dogParams.get_how_many_plot_resolutions(),
                    "plot_ry does not match up with nout_per_plot");
            for(int i=1;i<dogParams.get_how_many_plot_resolutions();i++)
            {
                assert_divides(plot_rx[i],mx);
                assert_divides(plot_ry[i],my);
            }
        }
        // if plot_rx and plot_ry are not defined,
        // then assume that plot_mx and plot_my are defined
        // and convert them to the equivalent plot_rx and plot_ry.
        // (should we bother to support this?)
        else
        {
            assert_eq(strlen(s_plot_rx),0);
            assert_eq(strlen(s_plot_ry),0);
            const char* s_plot_mx = grid["plot_mx"].c_str();
            const char* s_plot_my = grid["plot_my"].c_str();
            if(strlen(s_plot_mx)==0 || strlen(s_plot_my)==0)
            {
                eprintf("if nout_per_plot is nonempty "
                        "must set either plot_rx and plot_ry or plot_mx and plot_my");
            }
            int numel;
            numel = new_array_from_str(s_plot_mx,plot_rx,1);
            if(numel<0) eprintf("invalid syntax: plot_mx = %s", s_plot_mx);
            assert_printf(numel==dogParams.get_how_many_plot_resolutions(),
                    "plot_mx does not match up with nout_per_plot");
            //
            numel = new_array_from_str(s_plot_my,plot_ry,1);
            if(numel<0) eprintf("invalid syntax: plot_my = %s", s_plot_my);
            assert_printf(numel==dogParams.get_how_many_plot_resolutions(),
                    "plot_my does not match up with nout_per_plot");
            for(int i=1;i<=dogParams.get_how_many_plot_resolutions();i++)
            {
                assert_divides(plot_rx[i],mx);
                assert_divides(plot_ry[i],my);
                plot_rx[i]=mx/plot_rx[i];
                plot_ry[i]=my/plot_ry[i];
            }
        }
        plot_rx[0]=1;
        plot_ry[0]=1;
    }

    set_xlims(xlow, xhigh);
    set_ylims(ylow, yhigh);

    checkParameters();
    setDerivedParameters();
    reportParameters();

    is_initialized = true;
}

void DogParamsCart2::reportParameters()
{
    printf(  "   === parameters from [grid] ===");
    printf("\n   mx   : %d", mx   );
    printf("\n   my   : %d", my   );
    printf("\n   mbc  : %d", mbc  );
    printf("\n   xlow : %f", xlow );
    printf("\n   xhigh: %f", xhigh);
    printf("\n   ylow : %f", ylow );
    printf("\n   yhigh: %f", yhigh);
    printf("\n   --- parameters derived from [grid] ---");
    printf("\n   dx   : %f ", dx);
    printf("\n   dy   : %f ", dy);
    printf("\n\n");
    const int num = dogParams.get_how_many_plot_resolutions();
    if (num>0)
      {
	printf("   === parameters from [grid.plot] ===");    
	printf("\n   plot_rx = "); fprint_array(stdout,plot_rx,1,num,',',"");
	printf("\n   plot_ry = "); fprint_array(stdout,plot_ry,1,num,',',"");
	printf("\n\n");
      }
}

void DogParamsCart2::setDerivedParameters()
{
    update_dx();
    update_dy();
    prim_vol = dx*dy;
}

void DogParamsCart2::checkParameters()
{
    if(mx <= 0) eprintf("ERROR: mx=%d must be positive.\n",mx);
    if(my <= 0) eprintf("ERROR: my=%d must be positive.\n",my);
    if(mbc < 0) eprintf("ERROR: mbc=%d must be nonnegative.\n",mbc);
    if((xhigh-xlow) <= 0) {
        eprintf("ERROR: xhigh=%d should be greater than"
                " xlow=%d\n", xhigh, xlow);
    }
    if((yhigh-ylow) <= 0) {
        eprintf("ERROR: yhigh=%d should be greater than"
                " ylow=%d\n", yhigh, ylow);
    }
}

// data put into qhelp.dat (which is then used by plotting routines)
void DogParamsCart2::write_qhelp(const char* filename, int plot_mx, int plot_my)
{
    dogParams.write_qhelp(filename);
    FILE* file = fopen(filename,"a");
    fprintf(file,"%16d : mx\n",plot_mx);
    fprintf(file,"%16d : my\n",plot_my);
    fprintf(file,"%16.8e : xlow\n", xlow);
    fprintf(file,"%16.8e : xhigh\n",xhigh);
    fprintf(file,"%16.8e : ylow\n", ylow);
    fprintf(file,"%16.8e : yhigh\n",yhigh);
    fclose(file);
}

// the new way
void DogParamsCart2::write_output_parameters_ini(const char* filename,
        int plot_mx, int plot_my)
{
    const double plot_dx = (xhigh-xlow)/plot_mx;
    const double plot_dy = (yhigh-ylow)/plot_my;
    dogParams.write_output_parameters_ini(filename);
    FILE* file = fopen(filename,"a");
    fprintf(file,"[grid]\n");
    fprintf(file,"; mesh used to plot\n");
    fprintf(file,"xlow    = %.16g\n",xlow);
    fprintf(file,"ylow    = %.16g\n",ylow);
    fprintf(file,"xhigh   = %.16g\n",xhigh);
    fprintf(file,"yhigh   = %.16g\n",yhigh);
    fprintf(file,"plot_mx = %d\n",plot_mx);
    fprintf(file,"plot_my = %d\n",plot_my);
    fprintf(file,"plot_dx = %.16g\n",plot_dx);
    fprintf(file,"plot_dy = %.16g\n",plot_dy);
    fprintf(file,"; actual computational mesh\n");
    fprintf(file,"mx      = %d\n",mx);
    fprintf(file,"my      = %d\n",my);
    fprintf(file,"dx      = %.16g\n",dx);
    fprintf(file,"dy      = %.16g\n",dy);
    fprintf(file,"datafmt = %d\n",dogParams.get_datafmt());
    fclose(file);
}

void DogParamsCart2::write_plotting_help_files(
        int plot_mx, int plot_my, const char* outputdir)
{
    // create deprecated qhelp.dat file
    string qhelp = string(outputdir)+"/qhelp.dat";
    write_qhelp(qhelp.c_str(),plot_mx,plot_my);

    string filename = string(outputdir)+"/out_parameters.ini";
    write_output_parameters_ini(filename.c_str(),plot_mx,plot_my);
}
