#ifndef _FULLMATRIX_H_
#define _FULLMATRIX_H_

#include "tensors.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
using namespace std;

// FullMatrix Object ----------------------------------------
class FullMatrix
{
    public:
        FullMatrix(int n, int m);
        // Constructor
        // POST: Creates a matrix

        FullMatrix(const FullMatrix& another_matrix);
        // Copy constructor
        // POST: New FullMatrix created with size and contents same as another_matrix

        ~FullMatrix();
        // Destructor
        // POST: FullMatrix no longer exists

        const double& get(int i, int j) const;
        // Returns (i,j) entry

        void set(int i, int j, double input);
        // Sets (i,j) entry to value input

        void Sparsify();
        // get sparse information from FullMatrix

        // Various get functions
        const int& get_NumRows() const;
        const int& get_NumCols() const;
        const int& get_TotalNZ() const;
        const int& get_MaxNZrow() const;
        const bool& get_IsSparsify() const;
        const int& get_NZrow(int) const;
        const int& get_Index(int,int) const;
        const double& get_Value(int,int) const;

    private:
        // Basic constants
        int NumRows;        // Number of rows
        int NumCols;        // Number of columns
        int TotalNZ;        // Total number of nonzeros
        int MaxNZrow;       // Maximum over all rows of number of non-zeros in each row
        bool IsSparsify;    // Check if matrix has been sparsified

        // Arrays
        iTensor1* NZrow; // number of non-zeros in each row
        iTensor2* Index; // store index of non-zero entry
        dTensor2* Value; // store value of non-zero entry
        dTensor2* Matrix; // store full matrix
};
// ---------------------------------------------------------

#endif
