#ifndef DOGPARAMSUNST2_H
#define DOGPARAMSUNST2_H

// without these statements, i'm getting a compiler error -DS
#include <stdio.h>
#include <stdlib.h>

#include <string>
using namespace std;

class DogParamsUnst2
{
 private:
  int NumElems;
  int NumPhysElems;
  int NumGhostElems;
  int NumNodes;
  int NumPhysNodes;
  int NumBndNodes;
  int NumEdges;
  int NumBndEdges;
  int is_initialized;
  string outputdir;

 public:
  DogParamsUnst2();

  DogParamsUnst2(const DogParamsUnst2& old);

  ~DogParamsUnst2();

  void init(int NumElems_in,
	    int NumPhysElems_in,
	    int NumGhostElems_in,
	    int NumNodes_in,
	    int NumPhysNodes_in,
	    int NumBndNodes_in,
	    int NumEdges_in,
	    int NumBndEdges_in,
	    string outputdir_in);

  void write_qhelp(const char* filename);

  const int& get_NumElems() const;
  const int& get_NumPhysElems() const;
  const int& get_NumGhostElems() const;
  const int& get_NumNodes() const;
  const int& get_NumPhysNodes() const;
  const int& get_NumBndNodes() const;
  const int& get_NumEdges() const;
  const int& get_NumBndEdges() const;
  const string& get_outputdir() const;
};

extern DogParamsUnst2 dogParamsUnst2;

#endif
