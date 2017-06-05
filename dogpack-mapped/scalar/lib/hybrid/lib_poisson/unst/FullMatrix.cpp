// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (FullMatrix.cpp)
// --------------------------------------------------------------------------

#include "FullMatrix.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
using namespace std;

inline int iMax(int a, int b)
{ return a>b ? a : b; }

FullMatrix::FullMatrix(int in_NumRows, int in_NumCols)
// Constructor
// POST: Creates a FullMatrix
{
  NumRows     = in_NumRows;
  NumCols     = in_NumCols;
  TotalNZ     = 0;
  MaxNZrow    = 0;
  IsSparsify  = false;

  NZrow      = new iTensor1(NumRows);         // number of non-zeros in each row
  Matrix     = new dTensor2(NumRows,NumCols); // store full matrix
  Index      = new iTensor2(NumRows,NumCols);
  Value      = new dTensor2(NumRows,NumCols);

  for (int i=1; i<=NumRows; i++)
    for (int j=1; j<=NumCols; j++)
      {
	Matrix->set(i,j, 0.0 );
	Index->set(i,j, 0 );
	Value->set(i,j, 0.0 );
      }
}

FullMatrix::FullMatrix(const FullMatrix& another)
// Copy constructor
// POST: New mesh created with size and contents same as amesh
{
  NumRows    = another.NumRows;
  NumCols    = another.NumCols;
  TotalNZ    = another.TotalNZ;
  MaxNZrow   = another.MaxNZrow;
  IsSparsify = another.IsSparsify;
  
  for (int i=1; i<=NumRows; i++)
    {
      NZrow->set(i, another.NZrow->get(i) );
    }

  for (int i=1; i<=NumRows; i++)
    for (int j=1; j<=NumCols; j++)
      {
	Matrix->set(i,j, another.Matrix->get(i,j) );
	Index->set(i,j, another.Index->get(i,j) );
	Value->set(i,j, another.Value->get(i,j) );
      }
}

FullMatrix::~FullMatrix()
// Destructor
// POST: FullMatrix no longer exists
{
}

void FullMatrix::Sparsify()
{
  int iMax(int,int);
  TotalNZ = 0;
  MaxNZrow = 0;
  
  for (int i=1; i<=NumRows; i++)
    {
      int ntmp = 0;
      for (int j=1; j<=NumCols; j++)
	{
	  double tmp = Matrix->get(i,j);
	  if (fabs(tmp)>1.0e-15)
	    {
	      TotalNZ = TotalNZ + 1;
	      ntmp = ntmp + 1;

	      Index->set(i,ntmp, j   );
	      Value->set(i,ntmp, tmp );
	    }
	}
      NZrow->set(i, ntmp );
      MaxNZrow = iMax(MaxNZrow,ntmp);
    }
  
  IsSparsify = true;
}

const double& FullMatrix::get(int i, int j) const
{
  return Matrix->get(i,j);
}

void FullMatrix::set(int i, int j, double input)
{
  Matrix->set(i,j, input );
}

const int& FullMatrix::get_NumRows() const
{
  return NumRows;
}
 
const int& FullMatrix::get_NumCols() const
{
  return NumCols;
}

const bool& FullMatrix::get_IsSparsify() const
{
  return IsSparsify;
}

const int& FullMatrix::get_MaxNZrow() const
{
  if (IsSparsify==false)
    {
      printf("\n");
      printf(" ERROR in FullMatrix.cpp: Invalid function call, ");
      printf("first you need to call Sparsify().\n");
      printf("\n");
      exit(1);
    }
  return MaxNZrow;
}

const int& FullMatrix::get_TotalNZ() const
{
  if (IsSparsify==false)
    {
      printf("\n");
      printf(" ERROR in FullMatrix.cpp: Invalid function call, ");
      printf("first you need to call Sparsify().\n");
      printf("\n");
      exit(1);
    }
  return TotalNZ;
}

const int& FullMatrix::get_NZrow(int i) const
{
  if (IsSparsify==false)
    {
      printf("\n");
      printf(" ERROR in FullMatrix.cpp: Invalid function call, ");
      printf("first you need to call Sparsify().\n");
      printf("\n");
      exit(1);
    }
  return NZrow->get(i);
}

const int& FullMatrix::get_Index(int i,int j) const
{
  if (IsSparsify==false)
    {
      printf("\n");
      printf(" ERROR in FullMatrix.cpp: Invalid function call, ");
      printf("first you need to call Sparsify().\n");
      printf("\n");
      exit(1);
    }
  return Index->get(i,j);
}

const double& FullMatrix::get_Value(int i,int j) const
{
  if (IsSparsify==false)
    {
      printf("\n");
      printf(" ERROR in FullMatrix.cpp: Invalid function call, ");
      printf("first you need to call Sparsify().\n");
      printf("\n");
      exit(1);
    }
  return Value->get(i,j);
}
