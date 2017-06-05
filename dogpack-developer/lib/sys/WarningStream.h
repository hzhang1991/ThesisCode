#include <sstream>
using namespace std;
#ifndef DEBUG_H
#include "debug.h"
#endif

// Should this even exist?
//
class WarningStream
{
private:
  char *func;
  char *file;
  int line;
  std::stringstream sstream;
public:
  // assumes that file_ will not be modified!
  WarningStream(const char *func_, char *file_, int line_) :
    func((char*)func_), file(file_), line(line_) {};
  ~WarningStream() { Wprintf_fileLine(func, file, line, "%s", sstream.str().c_str()); }
  std::ostream & stream() { return sstream; }
private:
  WarningStream(const WarningStream &);
  WarningStream & operator=(const WarningStream &);
};
#define dwarn WarningStream(__func__, __FILE__, __LINE__).stream()
