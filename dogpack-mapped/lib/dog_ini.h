#ifndef DOG_INI_H
#define DOG_INI_H
#include <vector>
#include "IniDocument.h"
using namespace std;

// set defaults from defaults file
// ($DOGPACK/config/$(section_label)_defaults.ini)
void get_defaults(
  IniDocument::Section& params,
  const string& section_label,
  const vector<string>& option_names_list);

bool verify_options_set(
  const IniDocument::Section& params,
  const string& section_name,
  const vector<string>& option_names_list);

#endif
