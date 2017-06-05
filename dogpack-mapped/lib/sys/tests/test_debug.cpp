//#include "stdio.h"
// override the default maximum debug level defined in debug.h
#define MAX_DEBUG_LEVEL 2
#include "debug.h"
//#include "DebugStream.h"
//#include "ErrorStream.h"
#include "DebugStreams.h"

// an old way of doing things that I abandoned.
//
//void test_debug_stream1()
//{
//  debug << "Debug" << std::endl;
//  debug.On() = false;
//  debug << "This line will not print" << std::endl;
//}

void test_debug_stream2()
{
  dout << "dout: DebugLevel::get():" << DebugLevel::get();
  doutn(2) << "doutn(2)";
  dout1 << "dout1: " << 1.1;
  dout2 << "dout2: " << 2.2;
  dout3 << "dout3: " << 3;
  dwarn << "warning: something unignorable happened.";
}

void test_debug_printf()
{
  //debug_level = 1;
  dprintf("something high-level happened: %s has value %d", "variable", 9);
  dprintf1("something high-level happened: %s has value %d", "variable", 1);
  dprintf2("small loop output: %s has value %d", "variable", 2);
  dprintf3("deeply nested loop output: %s has value %d", "variable", 3);
  dprintfn(4,"verbose nested output: %s has value %d", "variable", 4);
  Wprintf("something unignorable happened: %s has value %d", "variable", 0);
}

void test_error()
{
  int four=4;
  invalid_value_error(four);
  // uncomment to test the following line
  //derr << "derr " << -1 << endl;
  eprintf(
    "something fatal happened: \n"
    "\t%s has value %d\n", "variable", -1);
}

int main() {

  DebugLevel::set(3);
  test_debug_stream2();
  test_debug_printf();
  test_error();
}

