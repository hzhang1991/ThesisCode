#include <stdlib.h>
#include "debug.h"
#include "dog_ini.h"

// access defaults from a file that other applications
// (e.g. matlab) might also want to read
static void get_defaults_helper(
        IniDocument::Section& params,
        const vector<string>& option_names_list,
        const string& section_label,
        const string& defaults_option_label,
        const string& defaults_file_in)
{
    string defaults_file = defaults_file_in;
    // search and replace $DOGPACK with value
    string S_DOGPACK("$DOGPACK");
    size_t found = defaults_file.find(S_DOGPACK);
    if(found!=string::npos)
    {
        // do we need the assignment here?:
        defaults_file = defaults_file.replace(found,
                S_DOGPACK.length(), string(getenv("DOGPACK")));
    }
    //ifstream default_ini_ifstream(defaults_file.c_str(), ios::in);
    IniDocument defaults_doc;
    //defaults_doc.fromInputStream(default_ini_ifstream);
    defaults_doc.initFromFile(defaults_file);
    IniDocument::Section& defaults_params = defaults_doc[section_label];

    vector<string>::iterator io;
    // have to cast option_names_list as non-const to iterate over it
    for(io = ((vector<string>&)option_names_list).begin();
            io != ((vector<string>&)option_names_list).end(); io++)
    {
        if(params[(*io)].empty())
        {
            // if the parameter is not set use the default.
            params.set((*io),defaults_params[(*io)]);
        }
    }
    // call self recursively in case the defaults file has defaults
    string new_defaults_file
        = defaults_params[defaults_option_label];
    if(!new_defaults_file.empty())
    {
        if(new_defaults_file == defaults_file_in)
        {
            eprintf("defaults file [%s] references itself\n", defaults_file_in.c_str());
        }
        get_defaults_helper(
                params,
                option_names_list,
                section_label,
                defaults_option_label,
                new_defaults_file);
    }
    return;
}

void get_defaults(
        IniDocument::Section& params,
        const string& section_label,
        const vector<string>& option_names_list)
{
    string defaults_option_label = "defaults_file";
    //string defaults_option_label = "dogParams_defaults";
    string defaults_file = params[defaults_option_label];
    // the default defaults file
    if(defaults_file.empty())
        defaults_file = "$DOGPACK/config/ini_defaults/"+section_label+".ini";
    get_defaults_helper(
            params,
            option_names_list,
            section_label,
            defaults_option_label,
            defaults_file);
}

bool verify_options_set(
        const IniDocument::Section& params,
        const string& section_name,
        const vector<string>& option_names_list)
{
    vector<string>::iterator io;
    for(io = ((vector<string>&)option_names_list).begin();
            io != ((vector<string>&)option_names_list).end(); io++)
    {
        const string& str = ((IniDocument::Section&)params)[(*io)];
        if(str.empty() ||
                str=="must_set")
        {
            eprintf("In parameters.ini, section [%s]"
                    " you must set the parameter %s\n",
                    section_name.c_str(),
                    (*io).c_str());
            return false;
        }
    }
    return true;
}

