// ===================================================================
// structs.h
// structs header file
//
// These structs are used in mesh.h
//
// ===================================================================

#ifndef STRUCTS_H
#define STRUCTS_H

// ===================================================================
struct point
{
  double x;
  double y;
  double z;
};
// ===================================================================
struct triangle
{
  int n1;
  int n2;
  int n3;
  int ctr;
};
// ===================================================================
struct bar
{
  int n1;
  int n2;
  int bt;
  int bt2;
};
// ===================================================================
struct tetra
{
  int n1;
  int n2;
  int n3;
  int n4;
};
// ===================================================================
struct cart2d
{
  double dx;
  double dy;
  double xlow;
  double ylow;
  int mx;
  int my;
};
// ===================================================================
struct cart3d
{
  double dx;
  double dy;
  double dz;
  double xlow;
  double ylow;
  double zlow;
  int mx;
  int my;
  int mz;
};
// ===================================================================
#endif
