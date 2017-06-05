#include <assert.h>
#include "DogParams.h"
#include "dog_ini.h"
#include "constants.h"
#include <stdlib.h>
#include <float.h> // for DBL_MAX
#include <cstring>
#include <string>
#include "dog_str.h"
#include "dog_io.h"
// for derr; eliminate using eprintf
#include "ErrorStream.h"

// -------------------------------------------------------------------------- //
// Is it possible to move these static function declarations elsewhere? (-DS)
static SplittingType::Enum get_splitting_type(const char* s_splitting);
static const char* print(SplittingType::Enum splitting);
static bool invalid_value(const char* varname, const char* val);
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Singleton used throughout the code for access to [dogParams]
DogParams dogParams;  // see DogParams.h, allocated as: extern DogParams dogParams;
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
//
// This function, ini() does the bulk work of scanning the section labeled
// [dogParams] from the parameters.ini file,
// and saving the information located there inside a singleton of this struct.
//
// Inside this section, there are parameters that other parts of the code
// know what their values are. For example, meqn stays constant through the
// simulation, but may parts need this.
//
// There are values that are derived from input parameters, and they too are
// saved inside this singleton struct.  For example, kmax and dA (cell area
// for Structured code), are common parameters that are needed in many places
// in the code, yet are not strictly specified in the parameters.ini file.
//
// This function computes these derived parameters, and sets conveneience
// accessors to other parameters.
//
// -------------------------------------------------------------------------- //
void DogParams::init()
{
    if(is_initialized) {
        return;
    }
    string section_label = "dogParams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    // list of the options with user-configurable defaults
    // (options not listed here are empty strings by default)
    vector<string> option_names_list;
    option_names_list.push_back("ndims"               );
    option_names_list.push_back("mesh_type"           );
    option_names_list.push_back("nout"                );
    option_names_list.push_back("tfinal"              );
    option_names_list.push_back("dtv(1)"              );
    option_names_list.push_back("dtv(2)"              );
    option_names_list.push_back("cflv(1)"             );
    option_names_list.push_back("cflv(2)"             );
    option_names_list.push_back("nv"                  );
    option_names_list.push_back("time_stepping_method");
    option_names_list.push_back("limiter_method"      );
    option_names_list.push_back("space_order"         );
    option_names_list.push_back("time_order"          );
    option_names_list.push_back("use_limiter"         );
    option_names_list.push_back("verbosity"           );
    option_names_list.push_back("mcapa"               );
    option_names_list.push_back("maux"                );
    option_names_list.push_back("source_term"         );
    option_names_list.push_back("flux_term"           );
    option_names_list.push_back("splitting"           );
    option_names_list.push_back("meqn"                );
    option_names_list.push_back("mrestart"            ); // deprecated
    option_names_list.push_back("nstart"              );
    option_names_list.push_back("nrestart"            );
    option_names_list.push_back("maintained_restart"  );
    option_names_list.push_back("report_frame_idx"    );
    option_names_list.push_back("datafmt"             );
    option_names_list.push_back("nout_per_restart"    );
    option_names_list.push_back("num_subintervals"    );
    option_names_list.push_back("use_divfree"         );
    option_names_list.push_back("ic_quad_order"       );
    //option_names_list.push_back("withPyClawPlotting"  );
    //option_names_list.push_back("nout_per_plot"       ); //default: empty string
    //option_names_list.push_back("which_compnt_divfree"); //default: empty string

    // get defaults from e.g. $DOGPACK/config/dogParams_defaults.ini
    get_defaults      (ini_sec, section_label, option_names_list);
    verify_options_set(ini_sec, section_label, option_names_list);

    const char* s_ndims                = ini_sec["ndims"               ].c_str();
    const char* s_mesh_type            = ini_sec["mesh_type"           ].c_str();
    const char* s_nout                 = ini_sec["nout"                ].c_str();
    const char* s_tfinal               = ini_sec["tfinal"              ].c_str();
    const char* s_dtv_1                = ini_sec["dtv(1)"              ].c_str();
    const char* s_dtv_2                = ini_sec["dtv(2)"              ].c_str();
    const char* s_cflv_1               = ini_sec["cflv(1)"             ].c_str();
    const char* s_cflv_2               = ini_sec["cflv(2)"             ].c_str();
    const char* s_nv                   = ini_sec["nv"                  ].c_str();
    const char* s_time_stepping_method = ini_sec["time_stepping_method"].c_str();
    const char* s_limiter_method       = ini_sec["limiter_method"      ].c_str();
    const char* s_space_order          = ini_sec["space_order"         ].c_str();
    const char* s_time_order           = ini_sec["time_order"          ].c_str();
    const char* s_use_limiter          = ini_sec["use_limiter"         ].c_str();
    const char* s_verbosity            = ini_sec["verbosity"           ].c_str();
    const char* s_mcapa                = ini_sec["mcapa"               ].c_str();
    const char* s_maux                 = ini_sec["maux"                ].c_str();
    const char* s_source_term          = ini_sec["source_term"         ].c_str();
    const char* s_flux_term            = ini_sec["flux_term"           ].c_str();
    const char* s_splitting            = ini_sec["splitting"           ].c_str();
    const char* s_meqn                 = ini_sec["meqn"                ].c_str();
    const char* s_mrestart             = ini_sec["mrestart"            ].c_str();
    const char* s_nstart               = ini_sec["nstart"              ].c_str();
    const char* s_nrestart             = ini_sec["nrestart"            ].c_str();
    const char* s_maintained_restart   = ini_sec["maintained_restart"  ].c_str();
    const char* s_report_frame_idx     = ini_sec["report_frame_idx"    ].c_str();
    const char* s_datafmt              = ini_sec["datafmt"             ].c_str();
    const char* s_nout_per_restart     = ini_sec["nout_per_restart"    ].c_str();
    const char* s_nout_per_plot        = ini_sec["nout_per_plot"       ].c_str();
    const char* s_num_subintervals     = ini_sec["num_subintervals"    ].c_str();
    const char* s_use_divfree          = ini_sec["use_divfree"         ].c_str();
    const char* s_which_compnt_divfree = ini_sec["which_compnt_divfree"].c_str();
    const char* s_ic_quad_order        = ini_sec["ic_quad_order"       ].c_str();
    //const char* s_withPyClawPlotting   = ini_sec["withPyClawPlotting"  ].c_str();

    int& space_order = method[1];
    int& time_order  = method[2];
    int& use_limiter = method[3];
    int& verbosity   = method[4];
    int& mcapa       = method[5];
    int& maux        = method[6];
    int& source_term = method[7];

    int datafmt_in;
    int nstart;
    int mrestart;

    // populate class with parameter data.
    sscanf(s_ndims      ,"%d" ,&ndims      )|| invalid_value("ndims"      , s_ndims      );
    sscanf(s_nout       ,"%d" ,&nout       )|| invalid_value("nout"       , s_nout       );
    sscanf(s_tfinal     ,"%lf",&tfinal     )|| invalid_value("tfinal"     , s_tfinal     );
    sscanf(s_dtv_1      ,"%lf",&dtv[1]     )|| invalid_value("dtv(1)"     , s_dtv_1      );
    sscanf(s_dtv_2      ,"%lf",&dtv[2]     )|| invalid_value("dtv(2)"     , s_dtv_2      );
    sscanf(s_cflv_1     ,"%lf",&cflv[1]    )|| invalid_value("cflv(1)"    , s_cflv_1     );
    sscanf(s_cflv_2     ,"%lf",&cflv[2]    )|| invalid_value("cflv(2)"    , s_cflv_2     );
    sscanf(s_nv         ,"%d" ,&nv         )|| invalid_value("nv"         , s_nv         );
    sscanf(s_space_order,"%d" ,&space_order)|| invalid_value("space_order", s_space_order);
    sscanf(s_time_order ,"%d" ,&time_order )|| invalid_value("time_order" , s_time_order );
    sscanf(s_use_limiter,"%d" ,&use_limiter)|| invalid_value("use_limiter", s_use_limiter);
    sscanf(s_verbosity  ,"%d" ,&verbosity  )|| invalid_value("verbosity"  , s_verbosity  );
    sscanf(s_mcapa      ,"%d" ,&mcapa      )|| invalid_value("mcapa"      , s_mcapa      );
    sscanf(s_maux       ,"%d" ,&maux       )|| invalid_value("maux"       , s_maux       );
    sscanf(s_source_term,"%d" ,&source_term)|| invalid_value("source_term", s_source_term);
    sscanf(s_meqn       ,"%d" ,&meqn       )|| invalid_value("meqn"       , s_meqn       );
    sscanf(s_mrestart   ,"%d" ,&mrestart   )|| invalid_value("mrestart"   , s_mrestart   );
    sscanf(s_nstart     ,"%d" ,&nstart     )|| invalid_value("nstart"     , s_nstart     );
    sscanf(s_nrestart   ,"%d" ,&nrestart   )|| invalid_value("nrestart"   , s_nrestart   );
    sscanf(s_maintained_restart,"%d" ,&maintained_restart)
        || invalid_value("maintained_restart", s_maintained_restart);
    sscanf(s_report_frame_idx,"%d" ,&report_frame_idx)
        || invalid_value("report_frame_idx", s_report_frame_idx);
    sscanf(s_datafmt    ,"%d" ,&datafmt_in )|| invalid_value("datafmt"    , s_datafmt    );
    sscanf(s_nout_per_restart,"%d", &nout_per_restart)
        || invalid_value("nout_per_restart" , s_nout_per_restart);
    sscanf(s_num_subintervals,"%d", &num_subintervals)
        || invalid_value("num_subintervals" , s_num_subintervals);
    sscanf(s_use_divfree,"%d", &use_divfree)|| invalid_value("use_divfree", s_use_divfree);
    //sscanf(s_withPyClawPlotting,"%d", &withPyClawPlotting)
    //  || invalid_value("withPyClawPlotting", s_withPyClawPlotting);
    sscanf(s_ic_quad_order,"%d", &ic_quad_order)|| invalid_value("ic_quad_order", s_ic_quad_order);

    datafmt   = (DataFmt)datafmt_in;
    mesh_type = strdup(s_mesh_type);
    time_stepping_method = strdup(s_time_stepping_method);

    // initialize and save the necessary space for the limiter (this will be
    // freed later)
    limiter_method = strdup(s_limiter_method);

    // flux_term
    if(!strcmp(s_flux_term,"0"))
        flux_term = false;
    else if(!strcmp(s_flux_term,"1"))
        flux_term = true;
    else
        invalid_value_error(s_flux_term);

    splitting = ::get_splitting_type(s_splitting);

    // parse lists of numbers
    //
    how_many_vectors_divfree = new_array_from_str(
            s_which_compnt_divfree,which_compnt_divfree,1,',');
    assert_printf(how_many_vectors_divfree>=0, "invalid syntax: %s = %s\n",
            "which_compnt_divfree", s_which_compnt_divfree);
    //
    how_many_plot_resolutions = new_array_from_str(
            s_nout_per_plot,nout_per_plot,1,',');
    assert_printf(how_many_plot_resolutions>=0,"invalid syntax: %s = %s\n",
            "nout_per_plot",s_nout_per_plot);
    if(how_many_plot_resolutions==0 && nout_per_restart!=1)
    {
        how_many_plot_resolutions = 1;
        delete nout_per_plot;
        nout_per_plot = new int[2];
        nout_per_plot[1]=1;
    }
    if(nout_per_plot) nout_per_plot[0]=nout_per_restart;

    // reset parameters based on couplings and remappings among parameters
    //
    if(dtv[1] > dtv[2])
    {
        Wprintf("reducing dtv[1] from %24.16e\n\t"
                "      to dtv[2] =    %24.16e",
                dtv[1], dtv[2]);
        dtv[1] = dtv[2];
    }
    if(how_many_vectors_divfree==0)
    {
        use_divfree = 0;
    }
    if(mrestart == 0) // deprecated
    {
        // for backwards compatibility
        nstart = -1;
    }
    if(nrestart >= 0 && nstart >=0)
    {
        dprint1(nrestart); dprint1(nstart);
        eprintf("cannot set both nrestart and nstart.")
    }
    if(nstart >=0)
    {
        assert_divides(nout_per_restart,nstart);
        nrestart = nstart/nout_per_restart;
    }

    resetParameters();
    setDerivedParameters();

    checkParameters();
    reportParameters();

    // remember if this struct has been initialized
    is_initialized = true;
}


// -------------------------------------------------------------------------- //
// Destructor.  There are only a handful of variables that are actually
// allocated.
DogParams::~DogParams()
{
    free(mesh_type);
    free(time_stepping_method);
    free(limiter_method);
    delete [] nout_per_plot;
    delete [] which_compnt_divfree;
    delete [] generic_components;
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Methods used for defining derived parameters
// -------------------------------------------------------------------------- //
void DogParams::set_kmax()
{
    switch(ndims)
    {
        case 1:
            kmax = get_space_order();
            break;
        case 2:
            kmax = (get_space_order()*(get_space_order()+1))/2;
            break;
        case 3:
            kmax = (get_space_order()*(get_space_order()+1)*(get_space_order()+2))/6;
            break;
        case 4:
            kmax = (get_space_order()*(get_space_order()+1)*(get_space_order()+2)*( get_space_order()+3 ))/24;
            break;

            // For future reference, in dimension d, a scheme with order M needs
            // kmax = nchoosek( M+d-1, d ) polynomials, where nchoosek is the
            // binomial coefficient.
            // Purely split schemes require slightly more unkowns.

        default:
            unsupported_value_error(ndims);
    }
}

// -------------------------------------------------------------------------- //
void DogParams::set_kmax_divfree()
{
    switch(ndims)
    {
        case 1:
            kmax_divfree = -1;
            break;
        case 2:
            kmax_divfree = get_space_order()*(get_space_order()+3)/2;
            break;
        case 3:
            kmax_divfree = -1;
            break;
        case 4:
            kmax_divfree = -1;
            break;
        default:
            unsupported_value_error(ndims);
    }
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
void DogParams::setDerivedParameters()
{
    set_kmax();
    set_kmax_divfree();
    // set generic_components lookup table
    // (and how_many_generic_components)
    //
    generic_components = new int[meqn+1];
    if(use_divfree)
    {
        for(int i=0;i<=meqn;i++) generic_components[i]=-1;
        int generic_component=0;
        for(int m=1;m<=meqn;m++)
        {
            bool is_generic = true;
            // check if m is not generic
            for(int i=1;i<=how_many_vectors_divfree;i++)
            {
                if(m==which_compnt_divfree[i] || m==which_compnt_divfree[i]+1)
                    is_generic=false;
            }
            if(is_generic) generic_components[++generic_component]=m;
        }
        how_many_generic_components = generic_component;
    }
    else
    {
        for(int i=0;i<=meqn;i++) generic_components[i]=i;
        how_many_generic_components = meqn;
    }
    // dprint(how_many_generic_components);

    frame_interval = tfinal/nout;
    if(num_subintervals > EPSILON)
        frame_subinterval = frame_interval/num_subintervals;
    else if(num_subintervals < -EPSILON) // infinity
        frame_subinterval = 0;
    else
        frame_subinterval = DBL_MAX; //numeric_limits<double>::max();
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Error message to give the user feedback about what the offending parameters
// did wrong.
static bool invalid_value(const char* varname, const char* val)
{
    eprintf("invalid value: %s = [%s]", varname, val);
    return false;
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
void DogParams::resetParameters()
{
    if (str_eq(time_stepping_method, "SDC") && get_time_order()==1)
    {
        free(time_stepping_method);
        time_stepping_method = strdup("Runge-Kutta");
        printf(
                "   NOTE:  1st order SDC is the same as 1st order Runge-Kutta. \n"
                "          The program will use time_stepping_method = %s\n\n",
                time_stepping_method);
    }
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
void DogParams::checkParameters()
{
    assert_divides(nout_per_restart,nout);

    // check cflv
    if (cflv[2] >= cflv[1])
    {
        derr << "desired CFL must be strictly smaller than max CFL\n";
    }

    int& space_order = method[1];
    int& time_order  = method[2];
    int& use_limiter = method[3];
    int& verbosity   = method[4];
    int& mcapa       = method[5];
    int& maux        = method[6];
    int& source_term = method[7];

    // check mcapa
    if (mcapa > maux)
    {
        derr << "mcapa cannot be larger than maux " << endl;
    }

    // check meqn
    if (meqn < 1)
    {
        derr << "meqn must be at least 1 " << endl;
    }

    if (use_limiter<0 || use_limiter>1)
    {
        derr << "use_limiter must be set either 0 or 1: use_limiter = " << use_limiter << endl;
    }

    // Extra checks based on which dimension is being solved.
    switch(ndims)
    {
        case 1:
            checkParameters1();
            break;
        case 2:
            checkParameters2();
            break;
        case 3:
            checkParameters3();
            break;
        case 4:
            checkParameters4();
            break;
        default:
            derr << " ERROR: number of dimensions must be set to 1, 2, 3 or 4.\n\n";
    }
}

// checks based on what has been implemented for one dimension
//
void DogParams::checkParameters1()
{
    if(get_mrestart())
    {
        derr << "DogParams: restarting from data has not yet "
            "been implemented for one dimension." << endl;
    }

    // check space_order
    if (get_space_order()<1 || get_space_order()>6 )
    {
        derr << "incorrect spatial accuracy," 
            << " must have space_order = 1, 2, 3, 4, 5, or 6" << endl;
    }  

    // check time_stepping_method
    if (!str_eq(time_stepping_method, "Runge-Kutta") &&
            !str_eq(time_stepping_method, "SDC") &&
            !str_eq(time_stepping_method, "Lax-Wendroff") &&
            !str_eq(time_stepping_method, "User-Defined") &&
            !str_eq(time_stepping_method, "LxW_synch") &&
            !str_eq(time_stepping_method, "LxW_asynch"))
    {
        derr << "the time-stepping method " << time_stepping_method 
            << " has not been implemented. " << endl << endl
            << " The available options are " << endl
            << "   1. Runge-Kutta " << endl
            << "   2. SDC " << endl
            << "   3. Lax-Wendroff " << endl
            << "   4. User-Defined " << endl;
    }

    // checks per time_stepping_method
    if (str_eq(time_stepping_method, "Runge-Kutta"))
    {
        if (get_time_order()<0 || get_time_order()>5)
        {
            derr << "Runge-Kutta must have time_order = 1, 2, 3, 4, or 5 " << endl;
        }
    }
    else if (str_eq(time_stepping_method, "SDC"))
    {
        if (get_time_order()<0 || get_time_order()>6)
        {
            derr << "SDC must have time_order = 1, 2, 3, 4, 5, or 6 " << endl;
        }
    }
    else if (str_eq(time_stepping_method, "Lax-Wendroff"))
    {
        if (get_time_order()<0 || (get_time_order()>5 ) )
        {
            derr << "Lax-Wendroff must have time_order = 1, 2, 3, 4, or 5 " << endl;
        }
    }

    // check limiter method
    if( !str_eq(limiter_method, "moment") &&
        !str_eq(limiter_method, "viscosity") &&
        !str_eq(limiter_method, "positive") &&
        !str_eq(limiter_method, "relax"))
    {
        derr << "the limiter method " << limiter_method
            << " has not been implemented. " << endl << endl
            << " The available options are " << endl
            << "   1. moment " << endl
            << "   2. viscosity " << endl
            << "   3. positive " << endl
            << "   4. relax " << endl;
    }


    // this option is not currently implemented in 1D (TODO - JR)
    ic_quad_order = get_space_order();

}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// checks based on what has been implemented for two dimensions
//
void DogParams::checkParameters2()
{
    // check space_order
    if (get_space_order()<1 || get_space_order()>5 )
    {
        derr << "invalid spatial accuracy,"
            << " must have space_order = 1, 2, 3, 4, or 5 " << endl;
    }  

    // check time_stepping_method
    if (!str_eq(time_stepping_method, "Runge-Kutta") &&
            !str_eq(time_stepping_method, "SDC") &&
            !str_eq(time_stepping_method, "Lax-Wendroff") && 
            !str_eq(time_stepping_method, "User-Defined") )
    {
        derr << "the time-stepping method " << time_stepping_method 
            << " has not been implemented. " << endl << endl
            << " The available options are " << endl
            << "   1. Runge-Kutta " << endl
            << "   2. SDC " << endl
            << "   3. Lax-Wendroff " << endl
            << "   4. User-Defined " << endl;
    }

    // check time_order per time_stepping_method
    if (str_eq(time_stepping_method, "Runge-Kutta") )
    {
        if (get_time_order()<0 || get_time_order()>5)
        {
            derr << "Runge-Kutta must have time_order = 1, 2, 3, 4, or 5 " << endl;
        }
    }
    else if (str_eq(time_stepping_method, "SDC"))
    {
        if (get_time_order()<0 || get_time_order()>5)
        {
            derr << "SDC must have time_order = 1, 2, 3, 4, or 5 " << endl;
        }
    }
    else if (str_eq(time_stepping_method, "Lax-Wendroff"))
    {
        if (get_time_order()<0 || (get_time_order()>3 ) )
        {
            derr << "Lax-Wendroff must have time_order = 1, 2, or 3 " << endl;
        }
    }

    // check limiter method
    if (!str_eq(limiter_method, "moment") &&
            !str_eq(limiter_method, "viscosity") &&
            !str_eq(limiter_method, "positive"))
    {
        derr << "the limiter method " << limiter_method
            << " has not been implemented. " << endl << endl
            << " The available options are " << endl
            << "   1. moment " << endl
            << "   2. viscosity " << endl
            << "   3. positive " << endl;
    }

    // Quadrature order for initial data projection
    if ( ic_quad_order == -100)
    {
        ic_quad_order = get_space_order();
    }
    else if ( ic_quad_order < 1 || ic_quad_order > 20 )
    {
        derr << " ic_quad_order must be an integer from 1 to 20,  ic_quad_order = " << ic_quad_order << endl;
    }

    if (ic_quad_order < get_space_order())
    {
        derr << " ic_quad_order must be >= space_order, ic_quad_order = " << ic_quad_order << ", space_order = " 
            << get_space_order() << endl;
    }

}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// checks based on what has been implemented for two dimensions
//
void DogParams::checkParameters3()
{
    // check space_order
    if (get_space_order()<1 || get_space_order()>4 )
    {
        derr << "invalid spatial accuracy,"
            << " must have space_order = 1, 2, 3, or 4 " << endl;
    }  

    // check time_stepping_method
    if (!str_eq(time_stepping_method, "Runge-Kutta"))// &&
        //      !str_eq(time_stepping_method, "SDC") &&
        //      !str_eq(time_stepping_method, "Lax-Wendroff") && 
        //      !str_eq(time_stepping_method, "User-Defined") )
    {
        derr << "the time-stepping method " << time_stepping_method 
            << " has not been implemented. " << endl << endl
            << " The only currently available option is Runge-Kutta." << endl;
        //       << " The available options are " << endl
        //       << "   1. Runge-Kutta " << endl
        //       << "   2. SDC " << endl
        //       << "   3. Lax-Wendroff " << endl
        //       << "   4. User-Defined " << endl;
    }

    // check time_order per time_stepping_method
    if (str_eq(time_stepping_method, "Runge-Kutta") )
    {
        if (get_time_order()<0 || get_time_order()>4)
        {
            derr << "Runge-Kutta must have time_order = 1, 2, 3, or 4 " << endl;
        }
    }
    /*
       else if (str_eq(time_stepping_method, "SDC"))
       {
       if (get_time_order()<0 || get_time_order()>5)
       {
       derr << "SDC must have time_order = 1, 2, 3, 4, or 5 " << endl;
       }
       }
       else if (str_eq(time_stepping_method, "Lax-Wendroff"))
       {
       if (get_time_order()<0 || (get_time_order()>3 ) )
       {
       derr << "Lax-Wendroff must have time_order = 1, 2, or 3 " << endl;
       }
       }
     */

    // check limiter method
    if (!str_eq(limiter_method, "moment") )
    {
        derr << "the limiter method " << limiter_method
            << " has not been implemented. " << endl << endl
            << " The only currently available option is moment." << endl;
        //            << " The available options are " << endl
        //            << "   1. moment " << endl
        //            << "   2. viscosity " << endl
        //            << "   3. positive " << endl;
    }

    // Quadrature order for initial data projection
    if ( ic_quad_order == -100)
    {
        ic_quad_order = get_space_order();
    }
    else if ( ic_quad_order < 1 || ic_quad_order > 20 )
    {
        derr << " ic_quad_order must be an integer from 1 to 20,  ic_quad_order = " << ic_quad_order << endl;
    }

    if (ic_quad_order < get_space_order())
    {
        derr << " ic_quad_order must be >= space_order, ic_quad_order = " << ic_quad_order << ", space_order = " 
            << get_space_order() << endl;
    }

}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// checks based on what has been implemented for three dimensions (6/9/2014)
//
void DogParams::checkParameters4()
{
    // check space_order
    if (get_space_order()<1 || get_space_order()>4 )
    {
        // TODO
        derr << "invalid spatial accuracy,"
            << " must have space_order = 1, 2, 3, or 4 " << endl;
    }  

    // check time_stepping_method
    if (!str_eq(time_stepping_method, "Runge-Kutta"))// &&
        //      !str_eq(time_stepping_method, "SDC") &&
        //      !str_eq(time_stepping_method, "Lax-Wendroff") && 
        //      !str_eq(time_stepping_method, "User-Defined") )
    {
        derr << "the time-stepping method " << time_stepping_method 
            << " has not been implemented. " << endl << endl
            << " The only currently available option is Runge-Kutta." << endl;
    }

    // check time_order per time_stepping_method
    if (str_eq(time_stepping_method, "Runge-Kutta") )
    {
        if (get_time_order()<0 || get_time_order()>4)
        {
            derr << "Runge-Kutta must have time_order = 1, 2, 3, or 4 " << endl;
        }
    }

    // check limiter method
    if (!str_eq(limiter_method, "moment"))// &&
    {
        derr << "the limiter method " << limiter_method
            << " has not been implemented. " << endl << endl
            << " The only currently available option is moment." << endl;
    }

    // Quadrature order for initial data projection
    if ( ic_quad_order == -100)
    {
        ic_quad_order = get_space_order();
    }
    else if ( ic_quad_order < 1 || ic_quad_order > 20 )
    {
        derr << " ic_quad_order must be an integer from 1 to 20,  ic_quad_order = " << ic_quad_order << endl;
    }

    if (ic_quad_order < get_space_order())
    {
        derr << " ic_quad_order must be >= space_order, ic_quad_order = " << ic_quad_order << ", space_order = " 
            << get_space_order() << endl;
    }

}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Print a friendly message to the screen describing a summary of what
// parameters were read in.
//
void DogParams::reportParameters()
{

    // Output parameters to screen
    printf(
            "   === parameters from [dogParams] ===\n"
            "                      Mesh Type:  %s\n"
            "   Number of Spatial Dimensions:  %d\n"
            "            Number of Equations:  %d\n"
            "      Order of Accuracy in Time:  %d\n"
            "     Order of Accuracy in Space:  %d\n"
            "           Time-Stepping Method:  %s\n"
            "                 Limiter Method:  %s\n"
            "                       Limiters:  %s\n"
            "         desired Courant number:  %f\n",
            mesh_type,ndims,meqn,
            get_time_order(),get_space_order(),
            time_stepping_method,limiter_method,
            get_use_limiter()?"yes":"no",cflv[2]);
    if(!get_source_term()||!get_flux_term()||!splitting)
    {

        printf(
                "                    source_term:  %d\n"
                "                      flux_term:  %d\n" 
                "                      splitting:  %s\n",
                get_source_term(),
                get_flux_term(),
                print(splitting));
    }
    if(nrestart>=0)
        printf(
                "                       nrestart:  %d\n",
                nrestart);
    if(maintained_restart>=0)
        printf(
                "             maintained_restart:  %d\n",
                maintained_restart);
    if(report_frame_idx>=0)
        printf(
                "               report_frame_idx:  %d\n",
                report_frame_idx);
    if(use_divfree){
        printf(
                "                    use_divfree:  %d\n"
                "       how_many_vectors_divfree:  %d\n"
                "           which_compnt_divfree:  ",
                use_divfree,how_many_vectors_divfree);
        fprint_array(stdout,which_compnt_divfree,1,how_many_vectors_divfree);
        printf(
                "    how_many_generic_components:  %d\n"
                "             generic_components:  ",
                how_many_generic_components);
        fprint_array(stdout,generic_components,1,how_many_generic_components);
    }
    printf("\n");
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
void DogParams::write_qhelp(const char* filename)
{
    FILE* file = fopen(filename,"w"); 
    fprintf(file,"%16d : ndims\n", ndims);
    fprintf(file,"%16s : mesh_type\n", mesh_type);    
    fprintf(file,"%16d : meqn\n", meqn);
    fprintf(file,"%16d : maux\n", get_maux());
    fprintf(file,"%16d : nout\n", nout);
    fprintf(file,"%16d : space_order\n", get_space_order());
    fprintf(file,"%16d : datafmt\n", get_datafmt());
    fclose(file);
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
void DogParams::write_output_parameters_ini(const char* filename)
{
    FILE* file = fopen(filename,"w");
    if( file == NULL )
    {
        printf("Error in opening filename: %s\n", filename );
        printf("Terminating program\n");
        exit(1);
    }
    fprintf(file,"[dogParams]\n");
    fprintf(file,"mesh_type = %s\n", mesh_type);
    fprintf(file,"meqn = %d\n", meqn);
    fprintf(file,"maux = %d\n", get_maux());
    fprintf(file,"nout = %d\n", nout);
    fprintf(file,"space_order = %d\n", get_space_order());
    fprintf(file,"datafmt = %d\n",get_datafmt());
    fclose(file);
}
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Limiter options:
// -------------------------------------------------------------------------- //
bool DogParams::using_positive_limiter()
{
    return (get_use_limiter()==1 && str_eq(get_limiter_method(),"positive"));
}

bool DogParams::using_moment_limiter()
{
    return (get_use_limiter()==1 && str_eq(get_limiter_method(),"moment"));
}

bool DogParams::using_viscosity_limiter()
{
    return (get_use_limiter()==1 && str_eq(get_limiter_method(),"viscosity"));
}

bool DogParams::using_relax_limiter()
{
    return (get_use_limiter()==1 && str_eq(get_limiter_method(),"relax"));
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// ??? I don't know what this part is for ... (-DS)
// -------------------------------------------------------------------------- //
static const char* print(SplittingType::Enum splitting)
{
    using namespace SplittingType;
    switch(splitting)
    {
        // These appear to be a list of splitting options. Can we please document
        // them?  (-DS)
        case none:
            return "none";
        case sf:
            return "sf";
        case fs:
            return "fs";
        case sfs:
            return "sfs";
        case fsf:
            return "fsf";
        default:
            unsupported_value_error(splitting);
    }
    return "error";
}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
//
// TODO - see above comment.  These need documentation. -DS
//
// -------------------------------------------------------------------------- //
static SplittingType::Enum get_splitting_type(const char* s_splitting)
{
    using namespace SplittingType;
    if(!strcmp(s_splitting,"none"))
        return none;
    if(!strcmp(s_splitting,"sf"))
        return sf;
    if(!strcmp(s_splitting,"fs"))
        return fs;
    if(!strcmp(s_splitting,"sfs"))
        return sfs;
    if(!strcmp(s_splitting,"fsf"))
        return fsf;
    eprintf("unsupported value for parameter splitting: %s", s_splitting);
    return none;
}
// -------------------------------------------------------------------------- //
