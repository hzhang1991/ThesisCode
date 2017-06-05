#ifndef IniDocument_h
#define IniDocument_h
#include <string>
#include <vector>

// class StrList
// {
//   char** list;
//  public:
//   StrList(int size) {
//     list = new char*[size];
//     for(int i=0;i<size;i++) list[i]=0;
//   }
//   // can use to set.
//   char* operator[](const int i) {
//     return list[i];
//   }
// }

// represent ini file as a list of Section objects
// each of which is a list of option objects.
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

extern IniDocument ini_doc;
#endif // IniDocument_h
