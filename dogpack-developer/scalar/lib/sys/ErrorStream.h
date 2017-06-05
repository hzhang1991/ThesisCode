#ifndef ERROR_STREAM_H
#define ERROR_STREAM_H
#include <sstream>
using namespace std;
#ifndef DEBUG_H
#include "debug.h"
#endif

//////////////////////////////////////////
// for user documentation see debug.cpp //
//////////////////////////////////////////

class ErrorStream
{
private:
  char *func;
  char *file;
  int line;
  std::stringstream sstream;
public:
  // assumes that file_ will not be modified!
  ErrorStream(const char *func_, const char *file_, int line_) :
    func((char*)func_), file((char*)file_), line(line_) {};
  ~ErrorStream() { eprintf_fileLine(func, file, line, "%s", sstream.str().c_str()); }
  std::ostream & stream() { return sstream; }
private:
  ErrorStream(const ErrorStream &);
  ErrorStream & operator=(const ErrorStream &);
};
#define derr ErrorStream(__func__, __FILE__, __LINE__).stream()
#endif
