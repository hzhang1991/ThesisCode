#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "EulerParams.h"

void EulerParams::init(IniDocument& ini_doc)
{
  string section_label = "eulerParams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [mhdparams]
  //
  string s_gamma             = ini_sec["gamma"         ];
  string s_x0                = ini_sec["x0"            ];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_gamma           (s_gamma         );
  istringstream is_x0              (s_x0            );

  //
  // Store information in global parameter values
  //
  gamma = 1.6666666666666667;
  x0 = 1.0/6.0;

  is_gamma >> gamma;  
  is_x0 >> x0;

  //
  // Output values to screen
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   gamma             :  " << gamma           << endl;
  cout << "   x0                :  " << x0              << endl;
  cout << endl;
}

void EulerParams::write_eulerhelp(const char* filename)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"%24.16e\n", gamma);
  fclose(file);
}
