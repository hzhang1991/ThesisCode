#ifndef SPARSECHOLESKY_H
#define SPARSECHOLESKY_H

#include <new>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <assert.h>
#include "tensors.h"
using namespace std;

// SparseCholesky Object ----------------------------------------
class SparseCholesky
{
 public:
    SparseCholesky(int n);
    // Constructor
    // POST: Creates a matrix
    
    SparseCholesky(const SparseCholesky& amat);
    // Copy constructor
    // POST: New SparseCholesky created with size and contents same as another_matrix
    
    ~SparseCholesky();
    // Destructor
    // POST: SparseCholesky no longer exists

    const double& get_Rval(int i, int j) const;
    // Returns jth non-zero entry on ith row

    const int& get_Rind(int i, int j) const;
    // Returns column index of jth non-zero entry on ith row

    const double& get_Lval(int i, int j) const;
    // Returns jth non-zero entry on ith row

    const int& get_Lind(int i, int j) const;
    // Returns column index of jth non-zero entry on ith row
    
    const int& get_Pval(int i) const;
    // Returns ith entry of permutation array

    void read(string inputdir);
    // Read in sparse matrix

    void write(string outputdir);
    // Write out sparse matrix

    void init(string outputdir);
    // Create necessary files in order to be able to read sparse matrix
    
    const int& get_size() const;
    // Returns matrix size (i.e., number of rows = number of columns)

    const int& get_num_nonzeros() const;
    // Returns total number of nonzero entries 
    //   0 <= num_nonzeros <= size*size

    const int& get_Rnz_row(int i) const;
    // Returns number of nonzeros on row i 
    //   0 <= num_nonzeros <= size

    const int& get_Lnz_row(int i) const;
    // Returns number of nonzeros on row i 
    //   0 <= num_nonzeros <= size

    const void ForwardSubs(const dTensor1& b, dTensor1& y) const;
    // Forward substitution using L

    const void BackwardSubs(const dTensor1& y, dTensor1& x) const;
    // Backward substitution using R

 private:
    // Basic constants
    int size;          // Number of rows/columns
    int num_nonzeros;  // Number of nonzeros
    
    // Arrays
    int* Rnz_row;
    int** Rind;
    double** Rval;
    int* Lnz_row;
    int** Lind;
    double** Lval;
    int* P;

};
// ---------------------------------------------------------

#endif
