#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "VlasovParams.h"

VlasovParams vlasovParams;

VlasovParams::~VlasovParams()
{

}

void VlasovParams::init()
{

printf("Initializing vlasovparams\n");

    if( is_initialized )
    { return; }

    string section_label = "vlasovparams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [vlasovparams]
    //
    string s_rho0              = ini_sec["rho0"];
    string s_area              = ini_sec["area"];
//  string s_temp              = ini_sec["temp"];
//  string s_eta               = ini_sec["eta"];
//  string s_omega0            = ini_sec["omega0"];
//  string s_n0                = ini_sec["n0"];
//  string s_vth2              = ini_sec["vth2"];
//  string s_a                 = ini_sec["a"];
//  string s_R2                = ini_sec["R2"];

    //
    // Convert read-in strings to numerical values
    //
    istringstream is_rho0            (s_rho0         );
    istringstream is_area            (s_area         );
//  istringstream is_temp            (s_temp         );
//  istringstream is_eta             (s_eta          );
//  istringstream is_omega0          (s_omega0       );
//  istringstream is_n0              (s_n0           );
//  istringstream is_vth2            (s_vth2         );
//  istringstream is_a               (s_a            );
//  istringstream is_R2              (s_R2           );

    //
    // Store information in global parameter values
    //

    // default values
    rho0   = 1.0;
    area   = 1.0;
//  temp   = 1.;

    // default values for the beam problem:
//  eta      = 0.25;
//  omega0   = 8.0;
//  n0       = 2.0*(1-eta*eta)*omega0*omega0;
//  vth2     = 1.0;
//  a        = 1.0;
//  R2       = 0.25;


//  is_area   >> area;  
//  is_temp   >> temp;  
//  is_rho0   >> rho0;  

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   rho_0             :  " << rho0          << endl;
    cout << "   area              :  " << area          << endl;
//  cout << "   temp              :  " << temp          << endl;
//  cout << "These parameters are not used in the code anywhere ... " << endl;
//  cout << "   eta               :  " << eta           << endl;
//  cout << "   omega0            :  " << omega0        << endl;
//  cout << "   n0                :  " << n0            << endl;
//  cout << "   vth2              :  " << vth2          << endl;
//  cout << "   a                 :  " << a             << endl;
//  cout << "   R2                :  " << R2            << endl;
    cout << endl;

    is_initialized = true;

}
