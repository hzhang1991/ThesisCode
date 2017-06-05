#include <stdlib.h> // for exit()
#include <sstream>
#include "DogParamsCart1.h"
#include "IniDocument.h"
#include "debug.h"
using namespace std;

DogParamsCart1 dogParamsCart1;

void DogParamsCart1::set_xlims(double xlow_in, double xhigh_in)
{
    xlow = xlow_in;
    xhigh = xhigh_in;
    dx = (xhigh-xlow)/mx;
}

static void not_set_error(string varname)
{
  eprintf("In parameters.ini, section [grid],"
    " you must set the parameter %s\n", varname.c_str());
}

static bool invalid_value(string varname, string val)
{
  eprintf("invalid value: %s = [%s]", varname.c_str(), val.c_str());
  return false;
}

// It does not make sense to define defaults for
// these options, so we do not use the defaults
// mechanism used e.g. in DogParams::init
void DogParamsCart1::init(IniDocument& ini_doc)
{
  if(is_initialized) {
    return;
  }

  IniDocument::Section& grid = ini_doc["grid"];
  //
  string s_mx    = grid["mx"   ];
  string s_mbc   = grid["mbc"  ];
  string s_xlow  = grid["xlow" ];
  string s_xhigh = grid["xhigh"];
  string s_read_grid = grid["read_grid"];
  //
  // set defaults (should we read these in from a file?)
  //
  if(s_mx      .empty()) not_set_error("mx");
  if(s_mbc     .empty()) not_set_error("mbc");
  if(s_xlow    .empty()) not_set_error("xlow");
  if(s_xhigh   .empty()) not_set_error("xhigh");
  if(s_read_grid.empty())
    { read_grid = 0; }
  else
    { 
      istringstream is_read_grid (s_read_grid );
      (is_read_grid  >> read_grid );
    }

  // convert strings to numbers
  //
  istringstream is_mx    (s_mx    );
  istringstream is_mbc   (s_mbc   );
  istringstream is_xlow  (s_xlow  );
  istringstream is_xhigh (s_xhigh );

  // populate class with parameter data.
  //
  (is_mx     >> mx    ) || invalid_value("mx"    , s_mx    );
  (is_mbc    >> mbc   ) || invalid_value("mbc"   , s_mbc   );
  (is_xlow   >> xlow  ) || invalid_value("xlow"  , s_xlow  );
  (is_xhigh  >> xhigh ) || invalid_value("xhigh" , s_xhigh );

  set_xlims(xlow, xhigh);

  checkParameters();
  setDerivedParameters();
  reportParameters();

  is_initialized = true;
}

void DogParamsCart1::reportParameters()
{
  printf(  "   === parameters from [grid] ===");
  printf("\n   mx   : %d", mx   );
  printf("\n   mbc  : %d", mbc  );
  printf("\n   xlow : %f", xlow );
  printf("\n   xhigh: %f", xhigh);
  printf("\n   --- parameters derived from [grid] ---");
  printf("\n   dx   : %f ", dx);
  printf("\n");
}

void DogParamsCart1::setDerivedParameters()
{
    dx = (xhigh-xlow)/mx;
}

void DogParamsCart1::checkParameters()
{
    if(mx <= 0) eprintf("ERROR: mx=%d must be positive.\n",mx);
    if(mbc < 0) eprintf("ERROR: mbc=%d must be nonnegative.\n",mbc);
    if((xhigh-xlow) <= 0) {
        eprintf("ERROR: xhigh=%d should be greater than"
                " xlow=%d\n", xhigh, xlow);
    }
}
