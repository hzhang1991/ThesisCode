#ifndef __DOG_IO_H__
#define __DOG_IO_H__

#include <stdio.h>

// Is there a higher-level system command to do a file copy?
void copyFile(const char* inFileName, const char* outFileName);

// conversion between strings and arrays
//
class iTensorBase;
class dTensorBase;
void fprint_tensor(FILE* file, const iTensorBase& t,
  char field_sep=',', const char*end="\n");
void fprint_tensor(FILE* file, const dTensorBase& t,
  char field_sep=',', const char*end="\n");
void fprint_array(FILE* file, const int* arr, int first_idx, int last_idx,
  char field_sep=',', const char*end="\n");

bool flagfile_check(const char* filename);
int flagfile_create(const char* filename);
int flagfile_remove(const char* filename);

#endif
