#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "AcousticParams.h"

void AcousticParams::init(IniDocument& ini_doc)
{

    string section_label = "acousticparams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [acousticparams]
    //
    string s_c  = ini_sec["c"];

    //
    // Convert read-in strings to numerical values
    //
    istringstream is_c ( s_c );

    //
    // Store information in global parameter values
    //

    // default values
    c = 1.0;
    set_c( c );

    is_c  >> c;  

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   c        :  " << c << endl;
    cout << endl;
}

const double& AcousticParams::get_c(void)
{ return c; }

void AcousticParams::set_c(double cin)
{ c = cin; }
