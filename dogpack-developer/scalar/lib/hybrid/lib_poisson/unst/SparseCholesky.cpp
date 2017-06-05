// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (SparseCholesky.cpp)
// --------------------------------------------------------------------------

#include "SparseCholesky.h"
using namespace std;

SparseCholesky::SparseCholesky(int size_in)
{
  size = size_in;
  num_nonzeros = 0;

  Rnz_row = new int[size+1];
  Rind = new int*[size+1];
  Rval = new double*[size+1];

  Lnz_row = new int[size+1];
  Lind = new int*[size+1];
  Lval = new double*[size+1];

  P = new int[size+1];
}

SparseCholesky::SparseCholesky(const SparseCholesky& amat)
{
  size = amat.size;
  num_nonzeros = amat.num_nonzeros;
  
  for (int i=1; i<=size; i++)
    {
      Rnz_row[i] = amat.Rnz_row[i];
      Lnz_row[i] = amat.Lnz_row[i];
      P[i] = amat.P[i];
    }

  for (int i=1; i<=size; i++)
    for (int j=1; j<=Rnz_row[i]; j++)
      {
	Rind[i][j] = amat.Rind[i][j];
	Rval[i][j] = amat.Rval[i][j];
      }

  for (int i=1; i<=size; i++)
    for (int j=1; j<=Lnz_row[i]; j++)
      {
	Lind[i][j] = amat.Lind[i][j];
	Lval[i][j] = amat.Lval[i][j];
      }
}

SparseCholesky::~SparseCholesky()
{
  delete[] Rnz_row;
  Rnz_row = NULL;

  delete[] Rind;
  Rind = NULL;

  delete[] Rval;
  Rval = NULL;

  delete[] Lnz_row;
  Lnz_row = NULL;

  delete[] Lind;
  Lind = NULL;

  delete[] Lval;
  Lval = NULL;

  delete[] P;
  P = NULL;
}


const double& SparseCholesky::get_Rval(int i, int j) const
// Returns jth non-zero entry on ith row
{
  return Rval[i][j];
}

const int& SparseCholesky::get_Rind(int i, int j) const
// Returns column index of jth non-zero entry on ith row
{
  return Rind[i][j];
}

const double& SparseCholesky::get_Lval(int i, int j) const
// Returns jth non-zero entry on ith row
{
  return Lval[i][j];
}

const int& SparseCholesky::get_Lind(int i, int j) const
// Returns column index of jth non-zero entry on ith row
{
  return Lind[i][j];
}

const int& SparseCholesky::get_Pval(int i) const
// Returns ith entry of permutation array
{
  return P[i];
}

void SparseCholesky::read(string inputdir)
// Read in sparse matrix
{
  string fname;

  // input: number of non-zeros in R
  fname = inputdir + "/matlab_output/Rnz.dat";
  ifstream read1(fname.c_str(), ios::in);
  num_nonzeros = 0;
  Rnz_row[0]=0;
  for (int i=1; i<=size; i++)
    {
      int tmp_int;

      read1 >> tmp_int;
      Rnz_row[i] = tmp_int;  

      num_nonzeros = num_nonzeros + tmp_int;

      Rind[i] = new int[tmp_int+1];
      Rval[i] = new double[tmp_int+1];
    }
  read1.close();
  read1.clear();
  
  // input: get all non-zeros in R (index and value)
  int* counter = new int[size+1];
  for (int i=0; i<=size; i++)
    {  counter[i] = 0;  }
  fname = inputdir + "/matlab_output/R.dat";
  ifstream read2(fname.c_str(), ios::in);
  for (int m=1; m<=num_nonzeros; m++)
    {
      int i,j;
      double tmp_double;
      read2 >> i >> j >> tmp_double;

      counter[i] = counter[i] + 1;

      Rind[i][counter[i]] = j;
      Rval[i][counter[i]] = tmp_double;
    }
  read2.close();
  read2.clear();

  // input: number of non-zeros in L
  fname = inputdir + "/matlab_output/Lnz.dat";
  ifstream read3(fname.c_str(), ios::in);
  int num_nonzeros_tmp = 0;
  Lnz_row[0]=0;
  for (int i=1; i<=size; i++)
    {
      int tmp_int;

      read3 >> tmp_int;
      Lnz_row[i] = tmp_int;  

      num_nonzeros_tmp = num_nonzeros_tmp + tmp_int;

      Lind[i] = new int[tmp_int+1];
      Lval[i] = new double[tmp_int+1];
    }
  read3.close();
  read3.clear();
  assert_eq(num_nonzeros,num_nonzeros_tmp);
  
  // input: get all non-zeros in L (index and value)
  for (int i=0; i<=size; i++)
    {  counter[i] = 0;  }
  fname = inputdir + "/matlab_output/L.dat";
  ifstream read4(fname.c_str(), ios::in);
  for (int m=1; m<=num_nonzeros; m++)
    {
      int i,j;
      double tmp_double;
      read4 >> i >> j >> tmp_double;

      counter[i] = counter[i] + 1;

      Lind[i][counter[i]] = j;
      Lval[i][counter[i]] = tmp_double;
    }
  read4.close();
  read4.clear();

  // input: get all permutation values
  fname = inputdir + "/matlab_output/P.dat";
  ifstream read5(fname.c_str(), ios::in);
  P[0]=0;
  for (int m=1; m<=size; m++)
    {
      int tmp_int;
      read5 >> tmp_int;

      P[m] = tmp_int;
    }
  read5.close();
  read5.clear();
}

void SparseCholesky::write(string outputdir)
// Write out sparse matrix
{
  string fname;

  // output: number of non-zeros in R
  fname = outputdir + "/matlab_output/Rnz_out.dat";
  ofstream write1(fname.c_str(), ios::out);
  for (int i=1; i<=size; i++)
    {
      write1 << setw(8) << Rnz_row[i] << endl;
    }
  write1.close();
  write1.clear();
  
  // output: write all non-zeros in R (index and value)
  fname = outputdir + "/matlab_output/R_out.dat";
  ofstream write2(fname.c_str(), ios::out);
  write2 << setprecision(16);
  for (int i=1; i<=size; i++)
    for (int k=1; k<=Rnz_row[i]; k++)
      {
	int j = Rind[i][k];
	write2 << setw(8) << i << setw(8) << j << setw(32) << scientific << Rval[i][k] << endl;
      }
  write2.close();
  write2.clear();

  // output: number of non-zeros in L
  fname = outputdir + "/matlab_output/Lnz_out.dat";
  ofstream write3(fname.c_str(), ios::out);
  for (int i=1; i<=size; i++)
    {
      write3 << setw(8) << Lnz_row[i] << endl;
    }
  write3.close();
  write3.clear();
  
  // output: write all non-zeros in L (index and value)
  fname = outputdir + "/matlab_output/L_out.dat";
  ofstream write4(fname.c_str(), ios::out);
  write4 << setprecision(16);
  for (int i=1; i<=size; i++)
    for (int k=1; k<=Lnz_row[i]; k++)
      {
	int j = Lind[i][k];
	write4 << setw(8) << i << setw(8) << j << setw(32) << scientific << Lval[i][k] << endl;
      }
  write4.close();
  write4.clear();

  // output: write all permutation values
  fname = outputdir + "/matlab_output/P_out.dat";
  ofstream write5(fname.c_str(), ios::out);
  P[0]=0;
  for (int m=1; m<=size; m++)
    {
      write5 << setw(8) << P[m] << endl;
    }
  write5.close();
  write5.clear();
}

void SparseCholesky::init(string outputdir)
{
  char command_str[1024];
  int numchars;
  int exit_status;
  string matlaboutputdir = outputdir + "/matlab_output";

  // script to call MATLAB routines if necessary
  numchars = snprintf(command_str,1024,
                      "${DOGPACK}/scripts/matlabcopy_script %s\n",
                      matlaboutputdir.c_str());
  //printf("%s\n",command_str);
  assert(numchars<1023);
  assert(numchars>0);
  exit_status = system(command_str);
}
    
const int& SparseCholesky::get_size() const
// Returns matrix size (i.e., number of rows = number of columns)
{
  return size;
}

const int& SparseCholesky::get_num_nonzeros() const
// Returns total number of nonzero entries 
//   0 <= num_nonzeros <= size*size
{
  return num_nonzeros;
}

const int& SparseCholesky::get_Rnz_row(int i) const
{
  return Rnz_row[i];
}

const int& SparseCholesky::get_Lnz_row(int i) const
{
  return Lnz_row[i];
}

const void SparseCholesky::ForwardSubs(const dTensor1& b, dTensor1& y) const
// Forward substitution using L = R^T
{
  for (int i=1; i<=size; i++)
    {
      int kend = Lnz_row[i];
      double denom = Lval[i][kend]; 
      double numer = b.get(P[i]);

      assert_eq(Lind[i][kend],i);
    
      for (int k=1; k<=(kend-1); k++)
	{
	  int j      = Lind[i][k];
	  double val = Lval[i][k];
	  numer = numer - val*y.get(j);
	}
      
      y.set(i, numer/denom );
    }
}

const void SparseCholesky::BackwardSubs(const dTensor1& y, dTensor1& x) const
// Backward substitution using R
{
  dTensor1 xtmp(size);

  for (int i=size; i>=1; i--)
    {
      double denom = Rval[i][1]; 
      double numer = y.get(i);

      assert_eq(Rind[i][1],i);
    
      for (int k=2; k<=Rnz_row[i]; k++)
	{
	  int j      = Rind[i][k];
	  double val = Rval[i][k];
	  numer = numer - val*xtmp.get(j);
	}
      
      xtmp.set(i, numer/denom );
    }

  for (int i=1; i<=size; i++)
    {
      x.set(P[i], xtmp.get(i) );
    }
  
}
