//#include "stdio.h"
// override the default maximum assert level defined in debug.h
#include <assert.h>
#include <string>
#include <cmath>
using namespace std;
declare_assert_errmsg(const string&,const string&);
// eliminate abort so that everything will get tested.
#define abort() (void(0))

void test_assert()
{
  int a=1;
  assert_printf (a<0,"assert_printf: a=%d",a);
  assert_printf3(a<0,"assert_printf3: a=%d",a);
  assert_printf2(a<0,"assert_printf2: a=%d",a);
  assert_printf1(a<0,"assert_printf1: a=%d",a);
  assert(a<0);
  assert3(a<0);
  assert2(a<0);
  assert1(a<0);
  assert_eq(a,3);
  assert_le(a,0);
  string mystr = "wow";
  const char* str1="hello";
  string str2 = "hello";
  assert_streq(str1,"world");
  assert_streq(str1,"hello");
  assert_streq(str2.c_str(),"world");
  assert_streq(str2.c_str(),"hello");
  //assert_eq(str2,"world");
  //assert_eq(str2,"hello");
}

#include "AssertStream.h"
#define abort() (void(0))

void test_assertStream()
{
  int a=1;
  int b=2;
  assert_err(a>0) << "assert_err: a=" << a;
}

int main() {
  assert_almost_eq(M_PI,M_PI+M_PI*1e-12);
  test_assert();
  test_assertStream();
}

