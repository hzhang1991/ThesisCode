#ifndef _IniDocument_h_
#define _IniDocument_h_

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "debug.h"
#include "getdelim.h"

// The purpose of this class is to define an interface to the parameters.ini
// document which is required for each application in DoGPack.  The singleton
// instance named ini_doc provides the required interface.
// 
// Inside the IniDocument class, each ini file is described as a list of 
// "Section" objects, each of which is a list of options objects.
//
// For examples of how to interact with this instance, see the following:
//
// [dogParams] is parsed in $(DOGPACK)/lib/1d/main_global.cpp;
// [dogParams] is parsed in $(DOGPACK)/lib/2d/main_global.cpp;
// [grid]      is parsed in $(DOGPACK)/lib/2d/struct/DogParamsCart2.cpp
//
// See also: $DOGPACK/lib/dog_ini.cpp
class IniDocument
{

    public:
        // each line in a section defines an option, e.g.
        // mynumber = 23.4e-6 ; 
        // mystring = "string\" with escaped quote and newline\\n" ; this is a comment
        class Option
        {
            public:
                std::string label;
                std::string value;
                std::string line;
            public:
                Option(
                        const std::string& label_in,
                        const std::string& value_in,
                        const std::string& line_in)
                    : label(label_in),
                    value(value_in),
                    line(line_in){}
                Option(){}
                ~Option(){}
                std::string getLabel(){return label;}
                std::string getValue(){return value;}
        };

        // a collection of options preceded by a header line such as:
        // [section_name]
        class Section
        {
            public:
                std::string name;
                std::string line;
                std::vector<class Option> option_list;
            public:
                Section(const std::string& name_in, const std::string& line_in):
                    name(name_in), line(line_in){}
                ~Section(){}
                std::string getSectionName(){return name;}
                void pushOption(
                        std::string label,
                        std::string value,
                        std::string line);
                void set(const std::string &label, const std::string &value);
                std::string operator[](const std::string &label);
                Section& operator=(const Section &rhs);
        };

    private:
        bool is_initialized_;
        std::vector<class Section> section_list;
        Section* currentSection;
    private:
        Section* createSection(std::string name, std::string comment);
        bool parse_line(const char* line);
    public:
        bool is_initialized(){return is_initialized_;}
        IniDocument():
            currentSection(0),
            is_initialized_(false){}
        ~IniDocument(){}
        void initFromFile(const std::string& iniFileName);
        Section & operator[](const std::string&);
};

// singleton global class member:
extern IniDocument ini_doc;

#endif // _IniDocument_h_
