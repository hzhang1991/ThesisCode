#include "tensors.h"

void QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  hi the region of array a that is to be sorted

    int i = lo;
    int j = hi;
    double x=a.get( (i+j)/2 );
    double h;
    int itmp;

    //  partition
    while(i<=j) 
    {           
        while (a.get(i)<x) 
        {i++;} 

        while (a.get(j)>x) 
        {j--;}

        if (i<=j)
        {
            h=a.get(i);
            a.set(i,a.get(j));
            a.set(j,h);

            itmp = index.get(i);
            index.set(i, index.get(j) );
            index.set(j, itmp );

            i++; j--;
        }
    }

    //  recursion
    if (lo<j) QuickSort_double(a, index, lo, j);
    if (i<hi) QuickSort_double(a, index, i, hi);
}
#include "tensors.h"

void QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  hi the region of array a that is to be sorted

    int i = lo;
    int j = hi;
    int x=a.get( (i+j)/2 );
    int h;
    int itmp;

    //  partition
    while(i<=j) 
    {           
        while (a.get(i)<x) 
        {i++;} 

        while (a.get(j)>x) 
        {j--;}

        if (i<=j)
        {
            h=a.get(i);
            a.set(i,a.get(j));
            a.set(j,h);

            itmp = index.get(i);
            index.set(i, index.get(j) );
            index.set(j, itmp );

            i++; j--;
        }
    }

    //  recursion
    if (lo<j) QuickSort_int(a, index, lo, j);
    if (i<hi) QuickSort_int(a, index, i, hi);
}
