#include "DogParamsCart3.h"

using namespace std;

// This file parses and provides options to access the [grid] section of the 
// parameters.ini file

DogParamsCart3 dogParamsCart3;

int DogParamsCart3::get_plot_mx(int idx) const 
{
  assert_divides(plot_rx[idx],mx);
  return mx/plot_rx[idx];
}

int DogParamsCart3::get_plot_my(int idx) const 
{
  assert_divides(plot_ry[idx],my);
  return my/plot_ry[idx];
}

int DogParamsCart3::get_plot_mz(int idx) const 
{
  assert_divides(plot_rz[idx],mz);
  return mz/plot_rz[idx];
}

void DogParamsCart3::update_dx()
{
  dx = (xhigh-xlow)/double(mx);
}

void DogParamsCart3::update_dy()
{
  dy = (yhigh-ylow)/double(my);
}

void DogParamsCart3::update_dz()
{
  dz = (zhigh-zlow)/double(mz);
}

void DogParamsCart3::reset_mx(int mx_in)
{
  mx = mx_in;
  update_dx();
  Wprintf("reset mx to %d\n", get_mx());
}

void DogParamsCart3::reset_my(int my_in)
{
    my = my_in;
    update_dy();
    Wprintf("reset my to %d\n", get_my());
}

void DogParamsCart3::reset_mz(int mz_in)
{
    mz = mz_in;
    update_dz();
    Wprintf("reset mz to %d\n", get_mz());
}

void DogParamsCart3::set_xlims(double xlow_in, double xhigh_in)
{
  xlow = xlow_in;
  xhigh = xhigh_in;
  update_dx();
}

void DogParamsCart3::set_ylims(double ylow_in, double yhigh_in)
{
  ylow = ylow_in;
  yhigh = yhigh_in;
  update_dy();
}

void DogParamsCart3::set_zlims(double zlow_in, double zhigh_in)
{
  zlow = zlow_in;
  zhigh = zhigh_in;
  update_dz();
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
void DogParamsCart3::init()
{
  if(is_initialized) 
    {
      return;
    }
  
  string section_label = "grid";
  IniDocument::Section& grid = ini_doc[section_label];
  
  // options for which we provide non-empty-string defaults
  vector<string> option_names_list;
  option_names_list.push_back("mx"   );
  option_names_list.push_back("my"   );
  option_names_list.push_back("mz"   );
  option_names_list.push_back("mbc"  );
  option_names_list.push_back("xlow" );
  option_names_list.push_back("xhigh");
  option_names_list.push_back("ylow" );
  option_names_list.push_back("yhigh");
  option_names_list.push_back("zlow" );
  option_names_list.push_back("zhigh");
  
  // get defaults from e.g. $DOGPACK/config/dogParams_defaults.ini
  void get_defaults(IniDocument::Section& params,
		    const string& section_label,
		    const vector<string>& option_names_list);
  get_defaults(grid, section_label, option_names_list);
  bool verify_options_set(const IniDocument::Section& params,
			  const string& section_name,
			  const vector<string>& option_names_list);
  verify_options_set(grid, section_label, option_names_list);
  
  const char* s_mx    = grid["mx"   ].c_str();
  const char* s_my    = grid["my"   ].c_str();
  const char* s_mz    = grid["mz"   ].c_str();
  const char* s_mbc   = grid["mbc"  ].c_str();
  const char* s_xlow  = grid["xlow" ].c_str();
  const char* s_xhigh = grid["xhigh"].c_str();
  const char* s_ylow  = grid["ylow" ].c_str();
  const char* s_yhigh = grid["yhigh"].c_str();
  const char* s_zlow  = grid["zlow" ].c_str();
  const char* s_zhigh = grid["zhigh"].c_str();
  
  // convert strings to numbers;
  // populate class with parameter data.
  //
  sscanf(s_mx   ,"%d" , &mx    ) || invalid_value("mx"    , s_mx    );
  sscanf(s_my   ,"%d" , &my    ) || invalid_value("my"    , s_my    );
  sscanf(s_mz   ,"%d" , &mz    ) || invalid_value("mz"    , s_mz    );
  sscanf(s_mbc  ,"%d" , &mbc   ) || invalid_value("mbc"   , s_mbc   );
  sscanf(s_xlow ,"%lf", &xlow  ) || invalid_value("xlow"  , s_xlow  );
  sscanf(s_xhigh,"%lf", &xhigh ) || invalid_value("xhigh" , s_xhigh );
  sscanf(s_ylow ,"%lf", &ylow  ) || invalid_value("ylow"  , s_ylow  );
  sscanf(s_yhigh,"%lf", &yhigh ) || invalid_value("yhigh" , s_yhigh );
  sscanf(s_zlow ,"%lf", &zlow  ) || invalid_value("zlow"  , s_zlow  );
  sscanf(s_zhigh,"%lf", &zhigh ) || invalid_value("zhigh" , s_zhigh );
  
  // read plot mesh sizes out of parameters plot_[mr][xy]
  const char* s_plot_rx = grid["plot_rx"].c_str();
  const char* s_plot_ry = grid["plot_ry"].c_str();
  const char* s_plot_rz = grid["plot_rz"].c_str();
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
	  //
	  numel = new_array_from_str(s_plot_rz,plot_rz,1);
	  if(numel<0) eprintf("invalid syntax: plot_rz = %s", s_plot_rz);
	  assert_printf(numel==dogParams.get_how_many_plot_resolutions(),
			"plot_rz does not match up with nout_per_plot");
	  //
	  for(int i=1;i<dogParams.get_how_many_plot_resolutions();i++)
            {
	      assert_divides(plot_rx[i],mx);
	      assert_divides(plot_ry[i],my);
	      assert_divides(plot_rz[i],mz);
            }
        }
      // if plot_rx, plot_ry, and plot_rz are not defined,
      // then assume that plot_mx, plot_my, plot_mz are defined
      // and convert them to the equivalent plot_rx, plot_ry, and plot_rz.
      // (should we bother to support this?)
      else
        {
	  assert_eq(strlen(s_plot_rx),0);
	  assert_eq(strlen(s_plot_ry),0);
	  assert_eq(strlen(s_plot_rz),0);
	  const char* s_plot_mx = grid["plot_mx"].c_str();
	  const char* s_plot_my = grid["plot_my"].c_str();
	  const char* s_plot_mz = grid["plot_mz"].c_str();
	  if(strlen(s_plot_mx)==0 || strlen(s_plot_my)==0 || strlen(s_plot_mz)==0)
            {
	      eprintf("if nout_per_plot is nonempty "
		      "must set either plot_rx, plot_ry, plot_rz or plot_mx, plot_my, and plot_mz");
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
	  //
	  numel = new_array_from_str(s_plot_mz,plot_rz,1);
	  if(numel<0) eprintf("invalid syntax: plot_mz = %s", s_plot_mz);
	  assert_printf(numel==dogParams.get_how_many_plot_resolutions(),
			"plot_mz does not match up with nout_per_plot");
	  //
	  for(int i=1;i<=dogParams.get_how_many_plot_resolutions();i++)
            {
	      assert_divides(plot_rx[i],mx);
	      assert_divides(plot_ry[i],my);
	      assert_divides(plot_rz[i],mz);
	      plot_rx[i]=my/plot_rx[i];
	      plot_ry[i]=my/plot_ry[i];
	      plot_rz[i]=mz/plot_rz[i];
            }
        }
      plot_rx[0]=1;
      plot_ry[0]=1;
      plot_rz[0]=1;
    }
  
  set_xlims(xlow, xhigh);
  set_ylims(ylow, yhigh);
  set_zlims(zlow, zhigh);
  
  checkParameters();
  setDerivedParameters();
  reportParameters();
  
  is_initialized = true;
}

void DogParamsCart3::reportParameters()
{
  printf(  "   === parameters from [grid] ===");
  printf("\n   mx   : %d", mx   );
  printf("\n   my   : %d", my   );
  printf("\n   mz   : %d", mz   );
  printf("\n   mbc  : %d", mbc  );
  printf("\n   xlow : %f", xlow );
  printf("\n   xhigh: %f", xhigh);
  printf("\n   ylow : %f", ylow );
  printf("\n   yhigh: %f", yhigh);
  printf("\n   zlow : %f", zlow );
  printf("\n   zhigh: %f", zhigh);
  printf("\n   --- parameters derived from [grid] ---");
  printf("\n   dx   : %f ", dx);
  printf("\n   dy   : %f ", dy);
  printf("\n   dz   : %f ", dz);
  printf("\n\n");
  const int num = dogParams.get_how_many_plot_resolutions();
  if (num>0)
    {
      printf("   === parameters from [grid.plot] ===");
      printf("\n   plot_rx = "); fprint_array(stdout,plot_rx,1,num,',',"");
      printf("\n   plot_ry = "); fprint_array(stdout,plot_ry,1,num,',',"");
      printf("\n   plot_rz = "); fprint_array(stdout,plot_rz,1,num,',',"");
      printf("\n\n");
    }
}

void DogParamsCart3::setDerivedParameters()
{
  update_dx();
  update_dy();
  update_dz();
  prim_vol = dx*dy*dz;
}

void DogParamsCart3::checkParameters()
{
  if(mx <= 0) eprintf("ERROR: mx=%d must be positive.\n",mx);
  if(my <= 0) eprintf("ERROR: my=%d must be positive.\n",my);
  if(mz <= 0) eprintf("ERROR: mz=%d must be positive.\n",my);
  if(mbc < 0) eprintf("ERROR: mbc=%d must be nonnegative.\n",mbc);

  if((xhigh-xlow) <= 0) {
    eprintf("ERROR: xhigh=%d should be greater than"
	    " xlow=%d\n", xhigh, xlow);
  }
  if((yhigh-ylow) <= 0) {
    eprintf("ERROR: yhigh=%d should be greater than"
	    " ylow=%d\n", yhigh, ylow);
  }
  if((zhigh-zlow) <= 0) {
    eprintf("ERROR: zhigh=%d should be greater than"
	    " zlow=%d\n", zhigh, zlow);
  }
}

// deprecated
void DogParamsCart3::write_qhelp(const char* filename, int plot_mx, int plot_my, int plot_mz)
{
  dogParams.write_qhelp(filename);
  FILE* file = fopen(filename,"a");
  fprintf(file,"%16d : mx\n",plot_mx);
  fprintf(file,"%16d : my\n",plot_my);
  fprintf(file,"%16d : mz\n",plot_mz);
  fprintf(file,"%16.8e : xlow\n", xlow);
  fprintf(file,"%16.8e : xhigh\n",xhigh);
  fprintf(file,"%16.8e : ylow\n", ylow);
  fprintf(file,"%16.8e : yhigh\n",yhigh);
  fprintf(file,"%16.8e : zlow\n", zlow);
  fprintf(file,"%16.8e : zhigh\n",zhigh);
  fclose(file);
}

// the new way
void DogParamsCart3::write_output_parameters_ini(const char* filename,
						 int plot_mx, 
						 int plot_my,
						 int plot_mz)
{
  const double plot_dx = (xhigh-xlow)/plot_mx;
  const double plot_dy = (yhigh-ylow)/plot_my;
  const double plot_dz = (zhigh-zlow)/plot_mz;
  dogParams.write_output_parameters_ini(filename);
  FILE* file = fopen(filename,"a");
  fprintf(file,"[grid]\n");
  fprintf(file,"; mesh used to plot\n");
  fprintf(file,"xlow    = %.16g\n",xlow);
  fprintf(file,"ylow    = %.16g\n",ylow);
  fprintf(file,"zlow    = %.16g\n",zlow);
  fprintf(file,"xhigh   = %.16g\n",xhigh);
  fprintf(file,"yhigh   = %.16g\n",yhigh);
  fprintf(file,"zhigh   = %.16g\n",zhigh);
  fprintf(file,"plot_mx = %d\n",plot_mx);
  fprintf(file,"plot_my = %d\n",plot_my);
  fprintf(file,"plot_mz = %d\n",plot_mz);
  fprintf(file,"plot_dx = %.16g\n",plot_dx);
  fprintf(file,"plot_dy = %.16g\n",plot_dy);
  fprintf(file,"plot_dz = %.16g\n",plot_dz);
  fprintf(file,"; actual computational mesh\n");
  fprintf(file,"mx      = %d\n",mx);
  fprintf(file,"my      = %d\n",my);
  fprintf(file,"mz      = %d\n",mz);
  fprintf(file,"dx      = %.16g\n",dx);
  fprintf(file,"dy      = %.16g\n",dy);
  fprintf(file,"dz      = %.16g\n",dz);
  fprintf(file,"datafmt = %d\n",dogParams.get_datafmt());
  fclose(file);
}

void DogParamsCart3::write_plotting_help_files(int plot_mx, 
					       int plot_my, 
					       int plot_mz,
					       const char* outputdir)
{
  string qhelp = string(outputdir)+"/qhelp.dat";
  write_qhelp(qhelp.c_str(),plot_mx,plot_my,plot_mz);
  
  string filename = string(outputdir)+"/out_parameters.ini";
  write_output_parameters_ini(filename.c_str(),plot_mx,plot_my,plot_mz);
}
