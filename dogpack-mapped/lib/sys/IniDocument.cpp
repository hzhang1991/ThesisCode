#include "IniDocument.h"

using namespace std;

// declare global instance of ini_doc
IniDocument ini_doc;

// === code to parse lines ===

enum LineType{ERROR=0,EMPTY,OPTION,SECTION_LABEL};

// character classes

static bool is_paren_character(char c)
{
  switch(c)
  {
    case '(':
    case ')':
    case '<':
    case '>':
    case '[':
    case ']':
    case '{':
    case '}':
      return true;
  }
  return false;
}

static bool is_separator_character(char c)
{
  switch(c)
  {
    case '/':
    case '\\':
    case '|':
      return true;
  }
  return false;
}

static bool is_punctuation_character(char c)
{
  switch(c)
  {
    case ',':
    case '.':
    case ':':
    case ';':
    case '-':
    case '_':
      return true;
  }
  return false;
}

static bool is_whitespace_character(char c)
{
  switch(c)
  {
    case ' ':
    case '\t':
    case '\n':
    case '\r':
      return true;
  }
  return false;
}

static bool is_comment_character(const char c) 
{
  switch(c)
  {
    case ';':
    case '#':
    case '%':
      return true;
  }
  return false;
}

static bool is_quote_character(char c)
{
  switch(c)
  {
    case '"':
    case '\'':
    case '`':
      return true;
  }
  return false;
}

static bool is_alphabetic_character(char c)
{
  if(('A' <= c && c <= 'Z') ||
     ('a' <= c && c <= 'z') ||
     c == '_')
  {
    return true;
  }
  return false;
}

static bool is_numeric_character(char c)
{
  if('0' <= c && c <= '9') return true;
  return false;
}

static bool is_alphanumeric_character(char c)
{
  if(is_alphabetic_character(c)) return true;
  if(is_numeric_character(c)) return true;
  return false;
}

static bool is_valid_name_character(char c)
{
  if(is_whitespace_character(c)) return false;
  if(is_quote_character(c)) return false;
  switch(c)
  {
    case '=':
    case '[':
    case ']':
      return false;
  }
  return true;
  // we could shrink this to the following
  if(is_alphanumeric_character(c)) return true;
  if(is_paren_character(c)) return true;
  if(is_punctuation_character(c)) return true;
  return false;
}

inline static bool is_valid_section_name_character(char c)
{
  return is_valid_name_character(c);
}

inline static bool is_valid_option_name_character(char c)
{
  return is_valid_name_character(c);
}

static void step_over_whitespace(const char*& ptr)
{
  while(*ptr && is_whitespace_character(*ptr)) ptr++;
}

static LineType get_lineType(const char* ptr)
{
  step_over_whitespace(ptr);
  if(!*ptr || is_comment_character(*ptr))
    return EMPTY;
  if(*ptr == '[')
    return SECTION_LABEL;
  if(is_alphabetic_character(*ptr))
    return OPTION;
  Wprintf("Could not evaluate line type: ptr=%s\n",ptr);
  return ERROR;
}

static bool get_section_name(const char* ptr, string& section_name)
{
  assert(*ptr=='[');
  ptr++;
  // ptr should now point to the beginning of the string.
  if(!is_alphabetic_character(*ptr)){
    Wprintf("first character is not alphabetic\n")
    return false;
  }
  size_t i=0;
  while(ptr[i] && is_valid_section_name_character(ptr[i])) i++;
  if(ptr[i]!=']'){
    Wprintf("first nonalphanumeric character must be ']'\n");
    return false;
  }
  section_name = string(ptr, i);
  dprintf3("section_name=%s\n",section_name.c_str());
  return true;
}

static void get_option_name(const char*& ptr, string& var_name)
{
  size_t i=0;
  // while(ptr[i] && is_alphanumeric_character(ptr[i])) i++;
  while(ptr[i] && is_valid_option_name_character(ptr[i])) i++;
  var_name = string(ptr, i);
  ptr+=i;
}

// copy characters from ptr to values skipping escaping backslashes
// until we encounter closing double quote.
// we expect the first character to be a double quote
bool get_value_from_quoted_string(const char* ptr, string& value)
{
  value = ptr;
  assert(*ptr=='"');
  *ptr++;
  int i=0;
  int j=0;
  for(;ptr[i] && ptr[i] != '"';i++)
  {
    // can escape any character
    if(ptr[i]=='\\'){
      i++;
    }
    value[j++] = ptr[i];
  }
  if(!ptr[i]){
    Wprintf("missing closing quote\n");
    return false;
  }
  value.resize(j);
  // confirm that any trailing non-whitespace begins
  // with a comment character.
  ++ptr+=i;
  step_over_whitespace(ptr);
  if(*ptr && !is_comment_character(*ptr)) 
  {
    dprintf3("ptr=%s\n", ptr)
    Wprintf("invalid comment character: %c\n", *ptr);
    return false;
  }
  return true;
}

static bool is_valid_unquoted_character(char c)
{
  // // expansive list: any
  // // nonspace noncomment character;
  // if(is_whitespace_character(c)) return false;
  // if(is_comment_character(c)) return false;
  // if(is_quote_character(c)) return false;
  // return true;
  if(is_alphanumeric_character(c)) return true;
  switch(c)
  {
    // at least should support numbers and lists thereof.
    case '.':
    case ',':
    case '-':
    case '+':
      return true;
  }
  return false;
}


static bool get_value_from_unquoted_string(const char* ptr, string& value)
{
  value = ptr;
  if(!*ptr)
  {
    Wprintf("option is missing value.\n");
    return false;
  }
  const char* str=ptr;
  int i=0;
  for(;ptr[i] && is_valid_unquoted_character(ptr[i]);i++)
  {
    value[i] = ptr[i];
  }
  value.resize(i);
  ptr+=i;
  step_over_whitespace(ptr);
  if(*ptr && !is_comment_character(*ptr)) 
  {
    Wprintf("invalid comment character: %c\n", *ptr);
    return false;
  }
  return true;
}


static bool get_option_name_and_value(
  const char* ptr, string& option_name, string& value)
{
  get_option_name(ptr, option_name);
  dprintf3("%s\n",ptr);
  // skip over whitespace preceeding equals sign
  step_over_whitespace(ptr);
  if(*ptr++!='='){
    Wprintf("missing equals sign.\n");
    return false;
  }
  // skip over whitespace trailing equals sign
  step_over_whitespace(ptr);
  if(!*ptr){
    Wprintf("no value given for option.\n");
    return false;
  }
  if(*ptr=='"'){
    return get_value_from_quoted_string(ptr, value);
  } else {
    return get_value_from_unquoted_string(ptr, value);
  }
  return false;
}

static ssize_t my_getline(char** lineptr, size_t* n, FILE* stream)
{
  ssize_t ret = getdelim (lineptr, n, '\n', stream);
  dprintf3("my_getline ret=%d\n",ret);
  return ret;
}

bool IniDocument::parse_line( const char* line)
{
  const char* ptr=line;
  step_over_whitespace(ptr);
  LineType lineType = get_lineType(ptr);
  // printf("parsing line: %s\n",line);
  switch(lineType)
  {
    case SECTION_LABEL:
     {
      string section_name;
      bool success = get_section_name(ptr, section_name);
      if(!success){
        Wprintf("Could not parse section name.\n");
        return false;
      }
      // create new section
      currentSection = IniDocument::createSection(section_name, line);
      dprintf3("new section name is [%s]\n",section_name.c_str());
     }
      break;
    case OPTION:
     {
      string option_name;
      string value;
      bool success = get_option_name_and_value(ptr, option_name, value);
      dprintf3("option '%s' has value '%s'\n",
        option_name.c_str(),value.c_str());
      if(!success) return false;
      currentSection->pushOption(option_name, value, line);
     }
      break;
    case EMPTY:
      break;
    case ERROR:
      return false;
      break;
    default:
      eprintf("Invalid line type: %d\n", lineType);
  }
  return true;
}

void IniDocument::initFromFile(const std::string& ini_filename)
{
  size_t bytes_read;
  size_t nbytes = 100;
  char *line;

  FILE *fp;
  long len;
  char *buf;
  fp=fopen(ini_filename.c_str(),"rb");
  if(!fp){
    eprintf("Could not open file %s\n",ini_filename.c_str());
  }

  line = (char *) malloc (nbytes + 1);
  // parse lines into ini document
  //
  // create default introductory section
  currentSection = IniDocument::createSection(
    "defaultSection", string());
  rewind(fp); //fseek(fp,0,SEEK_SET); // rewind the file.
  while(bytes_read = my_getline (&line, &nbytes, fp), bytes_read!=-1)
  {
    bool success = parse_line(line);
    if(!success) eprintf("could not parse line: %s\n", line);
  }

  is_initialized_=true;
  // clean up
  fclose(fp);
  free(line);
}

void IniDocument::Section::pushOption(string label, string value, string line)
{
  this->option_list.push_back(Option(label, value, line));
}

void IniDocument::Section::set(const std::string &label, const std::string &value)
{
  vector<class IniDocument::Option>::iterator i;
  for(i = this->option_list.begin(); i != this->option_list.end(); i++)
  {
    if( (*i).label == label)
    {
      (*i).value = value;
      return;
    }
  }
  /* the label does not exist, so add it */
  pushOption(label, value, string());
  return;
}

IniDocument::Section & IniDocument::Section::operator=(
  const IniDocument::Section &rhs)
{
  this->name = rhs.name;
  this->line = rhs.line;
  this->option_list = rhs.option_list;
  return *this;
}

string IniDocument::Section::operator[](const std::string& str)
{
  vector<class IniDocument::Option>::iterator i;
  for(i = this->option_list.begin(); i != this->option_list.end(); i++)
  {
    if( (*i).label == str) {
      return (*i).value;
    }
  }
  // if it doesn't exist then return an empty string
  return string();
}

// return the section with a given name
class IniDocument::Section & IniDocument::operator[] (const string& name_in)
{
  vector<IniDocument::Section>::iterator i;
  for(i = this->section_list.begin(); i != this->section_list.end(); i++)
  {
    if( (*i).name == name_in)
    {
      return *i;
    }
  }
  // if it doesn't exist then return the first (the default)
  i = this->section_list.begin();
  return *i;
}

IniDocument::Section* IniDocument::createSection(string name, string line)
{
  this->section_list.push_back(Section(name, line));
  return &(this->section_list.back());
}

