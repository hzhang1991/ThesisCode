#include "dog_math.h"
#include "dogdefs.h"

void test_minmod()
{
  double a[3];
  double b[3];
  a[2] = 2.;
  a[1] = 1.;
  a[0] = 0.;
  b[2] = -2.;
  b[1] = -1.;
  b[0] = -0.;
  for(int i=0;i<3;i++)
  {
    printf(
      " signbit(a[%d]) = %d"
      " signbit(b[%d]) = %d\n",
      i, signbit(a[i]),
      i, signbit(b[i])
    );
    for(int j=0;j<3;j++)
    {
      printf(
        " minmod(a[%d],b[%d]) = %f "
        " minmod(a[%d],a[%d]) = %f "
        " minmod(b[%d],b[%d]) = %f "
        " minmod(b[%d],a[%d]) = %f \n",
        i,j, minmod(a[i],b[j]),
        i,j, minmod(a[i],a[j]),
        i,j, minmod(b[i],b[j]),
        i,j, minmod(b[i],a[j])
      );
    }
  }
  for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
  for(int k=0;k<3;k++)
  {
      printf(
        " minmod(a[%d],a[%d],a[%d]) = %f "
        " minmod(b[%d],b[%d],b[%d]) = %f \n",
        i,j,k, minmod(a[i],a[j],a[k]),
        i,j,k, minmod(b[i],b[j],b[k]));
  }
  for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
  {
      printf(
        " minmod(a[ 0],a[%d],a[%d]) = %f "
        " minmod(a[%d],a[ 0],a[%d]) = %f "
        " minmod(a[%d],a[%d],a[ 0]) = %f "
        " minmod(b[ 0],b[%d],b[%d]) = %f "
        " minmod(b[%d],b[ 0],b[%d]) = %f "
        " minmod(b[%d],b[%d],b[ 0]) = %f \n",
        i,j, minmod(a[ 0],a[ i],a[ j]),
        i,j, minmod(a[ i],a[ 0],a[ j]),
        i,j, minmod(a[ i],a[ j],a[ 0]),
        i,j, minmod(b[ 0],b[ i],b[ j]),
        i,j, minmod(b[ i],b[ 0],b[ j]),
        i,j, minmod(b[ i],b[ j],b[ 0]));
  }
}

main()
{
  test_minmod();
}
