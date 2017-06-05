#include "DogParamsUnst2.h"

DogParamsUnst2 dogParamsUnst2;

DogParamsUnst2::DogParamsUnst2()
{
  NumElems        = -1;
  NumPhysElems    = -1;
  NumGhostElems   = -1;
  NumNodes        = -1;
  NumPhysNodes    = -1;
  NumBndNodes     = -1;
  NumEdges        = -1;
  NumBndEdges     = -1;
  is_initialized  =  0;
}

DogParamsUnst2::DogParamsUnst2(const DogParamsUnst2& old)
{
  NumElems        = old.NumElems;
  NumPhysElems    = old.NumPhysElems;
  NumGhostElems   = old.NumGhostElems;
  NumNodes        = old.NumNodes;
  NumPhysNodes    = old.NumPhysNodes;
  NumBndNodes     = old.NumBndNodes;
  NumEdges        = old.NumEdges;
  NumBndEdges     = old.NumBndEdges;
  is_initialized  = old.is_initialized;
}

DogParamsUnst2::~DogParamsUnst2()
{
  NumElems        = -1;
  NumPhysElems    = -1;
  NumGhostElems   = -1;
  NumNodes        = -1;
  NumPhysNodes    = -1;
  NumBndNodes     = -1;
  NumEdges        = -1;
  NumBndEdges     = -1;
  is_initialized  =  0;
}

void DogParamsUnst2::init(int NumElems_in,
			  int NumPhysElems_in,
			  int NumGhostElems_in,
			  int NumNodes_in,
			  int NumPhysNodes_in,
			  int NumBndNodes_in,
			  int NumEdges_in,
			  int NumBndEdges_in)
{
  if (is_initialized==0)
    {
      NumElems      = NumElems_in;
      NumPhysElems  = NumPhysElems_in;
      NumGhostElems = NumGhostElems_in;
      NumNodes      = NumNodes_in;
      NumPhysNodes  = NumPhysNodes_in;
      NumBndNodes   = NumBndNodes_in;
      NumEdges      = NumEdges_in;
      NumBndEdges   = NumBndEdges_in;

      printf("        === parameters from Mesh === \n");
      printf("                       NumElems:  %d \n", NumElems);
      printf("                   NumPhysElems:  %d \n", NumPhysElems);
      printf("                  NumGhostElems:  %d \n", NumGhostElems);
      printf("                       NumNodes:  %d \n", NumNodes);
      printf("                   NumPhysNodes:  %d \n", NumPhysNodes);
      printf("                    NumBndNodes:  %d \n", NumBndNodes);
      printf("                       NumEdges:  %d \n", NumEdges);
      printf("                    NumBndEdges:  %d \n", NumBndEdges);
      printf("\n");

      is_initialized = 1;
    }  
}

// data put into qhelp.dat (which is then used by plotting routines)
void DogParamsUnst2::write_qhelp(const char* filename)
{
  FILE* file = fopen(filename,"a");
  if (file!=NULL)
    {
      fprintf(file,"%16d : NumElems\n",NumElems);
      fprintf(file,"%16d : NumPhysElems\n",NumPhysElems);
      fprintf(file,"%16d : NumGhostElems\n",NumGhostElems);
      fprintf(file,"%16d : NumNodes\n",NumNodes);
      fprintf(file,"%16d : NumPhysNodes\n",NumPhysNodes);
      fprintf(file,"%16d : NumBndNodes\n",NumBndNodes);
      fprintf(file,"%16d : NumEdges\n",NumEdges);
      fclose(file);
    }
  else
    {
      printf("\n");
      printf(" Error in DogParamsUnst2::write_qhelp(const char* filename):\n");
      printf("     filename = %s cannot be found\n",filename);
      printf("\n");
      exit(1);
    }
}

const int& DogParamsUnst2::get_NumElems() const
{
  return NumElems;
}

const int& DogParamsUnst2::get_NumPhysElems() const
{
  return NumPhysElems;
}

const int& DogParamsUnst2::get_NumGhostElems() const
{
  return NumGhostElems;
}

const int& DogParamsUnst2::get_NumNodes() const
{
  return NumNodes;
}

const int& DogParamsUnst2::get_NumPhysNodes() const
{
  return NumPhysNodes;
}

const int& DogParamsUnst2::get_NumBndNodes() const
{
  return NumBndNodes;
}

const int& DogParamsUnst2::get_NumEdges() const
{
  return NumEdges;
}

const int& DogParamsUnst2::get_NumBndEdges() const
{
  return NumBndEdges;
}
