#include "dog_io.h"
#include "debug.h"
#define NDIMS 0
#include "tensors1d.h"
#include <unistd.h> // for unlink
#include <string>

// Is there a higher-level system command to do a file copy?
void copyFile(const char* inFileName, const char* outFileName)
{
  FILE *inFile = fopen(inFileName,"r");
  if(!inFile) eprintf("could not open file for reading: %s",inFileName);
  FILE *outFile = fopen(outFileName,"w");
  if(!outFile) eprintf("could not open file for writing: %s",outFileName);
  int c;
  while(c = fgetc(inFile), c!=EOF)
  {
    fputc(c,outFile);
  }
  fclose(inFile);
  fclose(outFile);
}

// returns true if file exists
//
// need to change this to call opendir so that the directory
// cache will be updated. otherwise this doesn't work over a
// networked file system, if matlab (the requesting process)
// is running on a different machine.  But then we should
// really be using tcpip to transfer data anyway rather than
// writing and reading files.
//
bool flagfile_check(const char* filename)
{
  //DIR *dp = opendir (get_outputdir());
  FILE* file = fopen(filename,"r");
  if(!file) return false;
  fclose(file);
  return true;
}

// returns 0 if successfully created flagfile
int flagfile_create(const char* filename)
{
  FILE* file = fopen(filename,"w");
  if(!file) return 1;
  fclose(file);
  return 0;
}

// return 0 if successfully removed flagfile
int flagfile_remove(const char* filename)
{
  int err = unlink(filename);
  return err;
}

void fprint_array(FILE* file, const int* arr, int first_idx, int last_idx,
  char field_sep, const char*end)
{
  if(last_idx-first_idx>=0) fprintf(file,"%d",arr[first_idx]);
  for(int i=first_idx+1;i<=last_idx;i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%d",arr[i]);
  }
  fprintf(file,"%s",end);
}

void fprint_tensor(FILE* file, const iTensorBase& t,
  char field_sep, const char*end)
{
  if(t.numel()>0) fprintf(file,"%d",t.vget(0));
  for(int i=1;i<t.numel();i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%d",t.vget(i));
  }
  fprintf(file,"%s",end);
}

void fprint_tensor(FILE* file, const dTensorBase& t,
  char field_sep, const char*end)
{
  if(t.numel()>0) fprintf(file,"%24.16e",t.vget(0));
  for(int i=0;i<t.numel();i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%24.16e",t.vget(i));
  }
  fprintf(file,"%s",end);
}
