#if !defined(DEBUGSTREAM_H_)
#define DEBUGSTREAM_H_

//////////////////////////////////////////
// for user documentation see debug.cpp //
//////////////////////////////////////////

using namespace std;

// Often the user wants to turn on debug in certain parts
// of the code.  For this the maximum debug level would be
// defined at compile time per file.
//
#ifndef DEBUG_H
#include "debug.h"
#endif
#ifdef DEBUG_ON_
  #include <iostream>
  #include <sstream>
#endif
#include "VoidStream.h"

// taken from http://stackoverflow.com/questions/1328568/custom-stream-manipulator-for-class
#ifdef DEBUG_ON_
class DebugStream
{
public:
  // assumes that file_ and func_will not be modified!
  DebugStream(short level_, const char * func_, const char * file_, int line_) :
    level(level_), func((char *)func_), file((char*)file_), line(line_) { }
  ~DebugStream() {
    dprintf_fileLine(level, func, file, line, "%s", sstream.str().c_str());
  }
  std::ostream & stream() { return sstream; }
private:
  short level;
  char *file;
  char *func;
  int line;
  std::stringstream sstream;
private:
  //DebugStream(const DebugStream &);
  DebugStream & operator=(const DebugStream &);
};
#endif // DEBUG_ON_

// When the anonymous DebugStream object goes out of scope
// (when control leaves the line) the destructor is called and
// the contents are written.
#if(MAX_DEBUG_LEVEL>=1)
  #define dout1 if (1<=DebugLevel::get()) \
    DebugStream(1, __func__, __FILE__, __LINE__).stream()
#else
  #define dout1 if(0) voidStream
#endif
#if(MAX_DEBUG_LEVEL>=2)
  #define dout2 if (2<=DebugLevel::get()) \
    DebugStream(2, __func__, __FILE__, __LINE__).stream()
#else
  #define dout2 if(0) voidStream
#endif
#if(MAX_DEBUG_LEVEL>=3)
  #define dout3 if (3<=DebugLevel::get()) \
    DebugStream(3, __func__, __FILE__, __LINE__).stream()
#else
  #define dout3 if(0) voidStream
#endif
#ifdef DEBUG_ON_
  #define doutn(n) if (n<=DebugLevel::get()) \
    DebugStream(n, __func__, __FILE__, __LINE__).stream()
#else
  #define doutn(n) if(0) voidStream
#endif

// debug level of dout is intended to be contextually defined;
// define a default debug level for dout
// (user is expected to override this definition)
#define dout dout2

#endif // DEBUGSTREAM_H_

