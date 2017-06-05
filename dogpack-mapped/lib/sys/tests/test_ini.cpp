#include "IniDocument.h"
//#include "ini.h"
#include <iostream>
using namespace std;

void test_ini(string inputfile)
{
  IniDocument iniDocument;
  iniDocument.initFromFile(inputfile);
  IniDocument::Section section1 = iniDocument["section1"];
  cout << "dtv(2)=[" << section1["dtv(2)"] << "]" << endl;
  cout << "method(2)=[" << section1["method(2)"] << "]" << endl;
  cout << "matlab_addpath=[" << section1["matlab_addpath"] << "]" << endl;
  IniDocument::Section section2 = iniDocument["section2"];
  cout << "mx=[" << section2["mx"] << "]" << endl;
  cout << "mystring=[" << section2["mystring"] << "]" << endl;
  cout << "inidoc[\"section1\"][\"time_stepping_method\"] = [" 
       << iniDocument["section1"]["time_stepping_method"] << "]" << endl;
  cout << "unquoted_list_of_numbers{3}=[" << section2["unquoted_list_of_numbers{3}"] << "]" << endl;
}

int main(int argc, char* argv[])
{
  string inputfile = "test_ini.ini";
  test_ini(inputfile);
  return 0;
}

