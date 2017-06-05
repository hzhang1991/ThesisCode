#ifndef DOGPARAMSUNST2_H
#define DOGPARAMSUNST2_H

#include <stdio.h>
#include <stdlib.h>

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
	    int NumBndEdges_in);

  void write_qhelp(const char* filename);

  const int& get_NumElems() const;
  const int& get_NumPhysElems() const;
  const int& get_NumGhostElems() const;
  const int& get_NumNodes() const;
  const int& get_NumPhysNodes() const;
  const int& get_NumBndNodes() const;
  const int& get_NumEdges() const;
  const int& get_NumBndEdges() const;
};

extern DogParamsUnst2 dogParamsUnst2;

#endif
