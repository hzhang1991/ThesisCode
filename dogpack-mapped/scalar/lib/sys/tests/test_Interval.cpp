#include "../Interval.h"

void test2()
{
  for(int i=-2;i<=2;i++)
  for(int j=-2;j<=2;j++)
  {
    Interval x(i,2);
    Interval y(j,2);
    x.print();
    printf(" * ");
    y.print();
    printf(" = ");
    Interval z = x*y;
    z.print();
    printf(" = ");
    x.multiply(y);
    x.print();
    printf("\n");
    assert_eq(x.min(),z.min());
    assert_eq(x.max(),z.max());
  }
  for(int i=-2;i<=2;i++)
  for(int j=-2;j<=2;j++)
  {
    Interval x(-2,i);
    Interval y(-2,j);
    x.print();
    printf(" * ");
    y.print();
    printf(" = ");
    Interval z = x*y;
    z.print();
    printf(" = ");
    x.multiply(y);
    x.print();
    printf("\n");
    assert_eq(x.min(),z.min());
    assert_eq(x.max(),z.max());
  }
}

void test1()
{
    Interval x(0,2);
    Interval y(-1,2);
    x.print();
    printf(" * ");
    y.print();
    printf(" = ");
    (x*y).print();
    printf(" = ");
    x.multiply(y);
    x.print();
    printf("\n");
}

void test3()
{
  Interval I(3.,4.);
  Interval J(7.,8.);
  Interval K=I+J;
  K.print();
  printf("\n");
  (I-J).print();
  printf("\n");
  (I*J).print();
  printf("\n");
  (I/J).print();
  printf("\n");
  I.reciprocal().print();
  printf("\n");
  Interval L=K-10.;
  L.print();
  printf("\n");
  printf("%d\n",(!(L>=0.)));
  printf("\n");
}

int main()
{
  Interval I(3.,4.);
  (-2.*I).print();
  printf("\n");
  test1();
  test2();
  test3();
}
