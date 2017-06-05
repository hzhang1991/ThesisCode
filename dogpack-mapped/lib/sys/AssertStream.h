#if !defined(AssertStream_h)
#define AssertStream_h

//////////////////////////////////////////
// for user documentation see debug.cpp //
//////////////////////////////////////////

using namespace std;

#include "assert.h"
#include "VoidStream.h"

// name chosen for compatibility with system assert.h
#if (defined(NDEBUG) || MAX_ASSERT_LEVEL < 1)
  // list of what we are supporting
  #define assert_err1(e) if(0) voidStream
  #define assert_err2(e) if(0) voidStream
  #define assert_err3(e) if(0) voidStream
  #define assert_err(e) if(0) voidStream
#else // ifndef NDEBUG

#include <iostream>
#include <sstream>
class AssertStream
{
public:
  // assumes that file_ and func_will not be modified!
  AssertStream(const char* expr_, const char * func_, char * file_, int line_) :
    expr((char*)expr_), func((char *)func_), file((char*)file_), line(line_) { }
  ~AssertStream() {
    dassert_printf_fileLine(expr, file, line, func, "%s", sstream.str().c_str());
  }
  std::ostream & stream() { return sstream; }
private:
  bool val;
  char *expr;
  char *file;
  char *func;
  int line;
  std::stringstream sstream;
private:
  AssertStream & operator=(const AssertStream &);
};

// When the anonymous AssertStream object goes out of scope
// (when control leaves the line) the destructor is called and
// the contents are written.
#if(MAX_ASSERT_LEVEL>=1)
  #define assert_err1(e) AssertStream(#e, __func__, __FILE__, __LINE__).stream()
#else
  #define assert_err1(e) if(0) voidStream
#endif
#if(MAX_ASSERT_LEVEL>=2)
  #define assert_err2(e) AssertStream(#e, __func__, __FILE__, __LINE__).stream()
  #define assert_err assert_err2
#else
  #define assert_err2(e) if(0) voidStream
  #define assert_err assert_err2
#endif
#if(MAX_ASSERT_LEVEL>=3)
  #define assert_err3(e) AssertStream(#e, __func__, __FILE__, __LINE__).stream()
#else
  #define assert_err3(e) if(0) voidStream
#endif

#endif // NDEBUG

#endif // AssertStream_h
