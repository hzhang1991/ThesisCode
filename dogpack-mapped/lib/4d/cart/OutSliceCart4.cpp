#include "OutSliceCart4.h"
#include <cmath>

// This file parses and provides options to access the [outslice] section of the 
// parameters.ini file

OutSliceCart4 outSliceCart4;

// Constructor
OutSliceCart4::OutSliceCart4()
{
  is_initialized = false;
  numslices = 0;
}

// Destructor
OutSliceCart4::~OutSliceCart4()
{
  delete[] slicekind;
  delete[] sliceidx;
}

// initialize (i.e. parse and save) all the parameters from [outslice]
// section of parameters.ini file.  This routine remembers if all these
// parameters have been previously saved.
void OutSliceCart4::init()
{
  if(is_initialized) 
    {
      return;
    }
  
  string section_label = "outslice";
  IniDocument::Section& outslice = ini_doc[section_label];
  
  // options for which we provide non-empty-string defaults
  vector<string> option_names_list;
  option_names_list.push_back("numslices");
  const char* s_numslices = outslice["numslices"].c_str();
  sscanf(s_numslices,"%d",&numslices );  

  // read all slicekind and sliceidx values
  const char* s_slicekind = outslice["slicekind"].c_str();
  const char* s_sliceidx  = outslice["sliceidx" ].c_str();
  if(numslices!=0)
    {
      int numel;
      numel = new_array_from_str(s_slicekind,slicekind,1);
      if(numel<0) eprintf("invalid syntax: slicekind = %s", s_slicekind);
      numel = new_array_from_str(s_sliceidx,sliceidx,1);
      if(numel<0) eprintf("invalid syntax: sliceidx = %s", s_sliceidx);
    }
  output4D = 1;
  const char* s_output4D  = outslice["output4D" ].c_str();
  sscanf(s_output4D,"%d",&output4D );
  
  checkParameters();
  reportParameters();
  
  is_initialized = true;
}

// Check parameters
void OutSliceCart4::checkParameters()
{
  if (numslices<0)
    {
      eprintf("ERROR: numslices=%d must be non-negative.\n",numslices);
    }

  for (int n=1; n<=numslices; n++)
    {
      if (slicekind[n]<1 || slicekind[n]>3)
	{
	  eprintf("ERROR: slicekind[%d]=%d must be 1,2, or 3.\n",n,slicekind[n]);
	}
      int mm = 0;
      switch(slicekind[n])
	{
	case 1:
	  mm = dogParamsCart4.get_mz();
	  break;
	case 2:
	  mm = dogParamsCart4.get_my();
	  break;
	case 3:
	  mm = dogParamsCart4.get_mx();
	  break;
	}
      if (sliceidx[n]<1 || sliceidx[n]>mm)
	{
	  eprintf("ERROR: sliceidx[%d]=%d must be between 1 and %d.\n",n,sliceidx[n],mm);
	}
    }

  if (output4D<0 || output4D>1)
    {
      eprintf("ERROR: output4D=%d must be 0 or 1.\n",output4D);
    }

}

// Prints a list of the parameters to standard output:
void OutSliceCart4::reportParameters()
{
  printf(  "   === parameters from [outslice] ===");
  printf("\n   numslices :  %d", numslices);     
  for (int i=1; i<=numslices; i++)
    {      
      switch(slicekind[i])
	{
	case 1:
	  printf("\n   slice %d   :  z = constant with index k = %d",i,sliceidx[i]);
	  break;
	case 2:
	  printf("\n   slice %d   :  y = constant with index j = %d",i,sliceidx[i]);
	  break;
	case 3:
	  printf("\n   slice %d   :  x = constant with index i = %d",i,sliceidx[i]);
	  break;
	case 4:
	  printf("\n   slice %d   :  w = constant with index l = %d",i,sliceidx[i]);
	  break;
	}
    }
  printf("\n    output4D :  %d", output4D);
  printf("\n\n");
}
  
// Accessors:
bool OutSliceCart4::get_is_initialized() const
{
  return is_initialized;
}

int OutSliceCart4::get_numslices() const
{
  if (is_initialized)
    { return numslices; }
  else
    { eprintf("ERROR: outSliceCart4 has not yet been initialized\n"); }
}

int OutSliceCart4::get_slicekind(int idx) const
{
  if (is_initialized)
    { return slicekind[idx]; }
  else
    { eprintf("ERROR: outSliceCart4 has not yet been initialized\n"); }
}
 
int OutSliceCart4::get_sliceidx(int idx) const
{
  if (is_initialized)
    { return sliceidx[idx]; }
  else
    { eprintf("ERROR: outSliceCart4 has not yet been initialized\n"); }
}

int OutSliceCart4::get_output4D() const
{
  return output4D;
}

//  
// Accessors for plotting routines:
//

void OutSliceCart4::write_plotting_help_files(const char* outputdir)
{
  string qhelp_slice = string(outputdir)+"/qhelp_slice.dat";
  write_qhelp_slice(qhelp_slice.c_str());
}

void OutSliceCart4::write_qhelp_slice(const char* filename)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"%4d : numslices\n",numslices);
  for (int ns=1; ns<=numslices; ns++)
    {
      fprintf(file,"%4d : slicekind[%d]\n",slicekind[ns],ns);
      fprintf(file,"%4d : sliceidx[%d]\n",sliceidx[ns],ns);
    }
  fclose(file);
}

// from DogSolver.h
const char* get_outputdir();

// from DogSolver.cpp
string get_varframe_filename(const char* dir, const char* varname, int nframe);

// Output function
void OutSliceCart4::write_output(const int noutput,
				 const double t,
				 const dTensorBC6& q,
				 const dTensorBC6& aux)
{

  string fileq = get_varframe_filename(get_outputdir(), "q", noutput);
  string filea = get_varframe_filename(get_outputdir(), "a", noutput);
  FILE* file;

/*
  for (int ns=1; ns<=numslices; ns++)
    {      
      assert_ge(ns,0);
      assert_le(ns,9999);

      char nns[5];
      snprintf(nns,5,"%04d",ns);

      string fileq_slice = fileq + "_slice" + nns + ".dat";
      
      dTensor4* qslice;
      slice_q(ns,&q,qslice);

      file = fopen(fileq_slice.c_str(),"w");
      fprintf(file,"%24.16e\n",t);
      for (int n=1; n<=qslice->getsize(4); n++)
	for (int m=1; m<=qslice->getsize(3); m++)
	  for (int j=1; j<=qslice->getsize(2); j++)
	    for (int i=1; i<=qslice->getsize(1); i++)
	      {
		double tmp = qslice->get(i,j,m,n);
		fprintf(file,"%24.16e\n",tmp);
	      }
      fclose(file);    

      if (dogParams.get_maux()>0)
	{
	  string filea_slice = filea + "_slice" + nns + ".dat";

	  dTensor4* aslice;
	  slice_q(ns,&aux,aslice);

	  file = fopen(filea_slice.c_str(),"w");
	  fprintf(file,"%24.16e\n",t);
	  for (int n=1; n<=aslice->getsize(4); n++)
	    for (int m=1; m<=aslice->getsize(3); m++)
	      for (int j=1; j<=aslice->getsize(2); j++)
		for (int i=1; i<=aslice->getsize(1); i++)
		  {
		    double tmp = aslice->get(i,j,m,n);
		    fprintf(file,"%24.16e\n",tmp);
		  }
	  fclose(file);  
	}      
    }

  // Optional extra output
  OutSliceExtra(noutput,t);
*/

}

void OutSliceCart4::slice_q(const int slicenumber,
			    const dTensorBC6* q, dTensor4*& qslice)
{

/*
  const int mx   = q->getsize(1);
  const int my   = q->getsize(2);
  const int mz   = q->getsize(3);
  const int meqn = q->getsize(4);
  const int kmax = q->getsize(5);

  int kmax2d;
  switch(kmax)
    {
    case 1:
      kmax2d = 1;
      break;
    case 4:
      kmax2d = 3;
      break;
    case 10:
      kmax2d = 6;
      break;
    case 20:
      kmax2d = 10;
      break;
    default:
      // should never get here
      eprintf("Error in slice_q, kmax must be 1,4,10,20.  kmax = %d\n",kmax);
    }

  int m1 = 0;
  int m2 = 0;
  switch(slicekind[slicenumber])
    {
    case 1:
      m1 = mx;
      m2 = my;
      break;
    case 2:
      m1 = mx;
      m2 = mz;
      break;
    case 3:
      m1 = my;
      m2 = mz;
      break;
    default:
      // should never get here
      eprintf("Error in SliceQ, slicekdind must be 1,2, or 3.  slicekind = %d\n",slicekind[slicenumber]);
    }

  qslice = new dTensor4(m1,m2,meqn,kmax2d);

  switch(slicekind[slicenumber])
    {
    case 1:
      for (int i=1; i<=mx; i++)
	for (int j=1; j<=my; j++)
	  for (int m=1; m<=meqn; m++)
	    {
	      dTensor1 qin(kmax);
	      dTensor1 qout(kmax2d);
	      
	      for (int n=1; n<=kmax; n++)
		{ qin.set(n, q->get(i,j,sliceidx[slicenumber],m,n) ); }
	      
	      OutSliceCart4Inline::slice_q_local_1(kmax2d,qin,qout);
	      
	      for (int n=1; n<=kmax2d; n++)
		{ qslice->set(i,j,m,n, qout.get(n) ); }
	    }
      break;

    case 2:
      for (int i=1; i<=mx; i++)
	for (int k=1; k<=mz; k++)
	  for (int m=1; m<=meqn; m++)
	    {
	      dTensor1 qin(kmax);
	      dTensor1 qout(kmax2d);
	      
	      for (int n=1; n<=kmax; n++)
		{ qin.set(n, q->get(i,sliceidx[slicenumber],k,m,n) ); }
	      
	      OutSliceCart4Inline::slice_q_local_2(kmax2d,qin,qout);
	      
	      for (int n=1; n<=kmax2d; n++)
		{ qslice->set(i,k,m,n, qout.get(n) ); }
	    }
      break;

    case 3:
      for (int j=1; j<=my; j++)
	for (int k=1; k<=mz; k++)
	  for (int m=1; m<=meqn; m++)
	    {
	      dTensor1 qin(kmax);
	      dTensor1 qout(kmax2d);
	      
	      for (int n=1; n<=kmax; n++)
		{ qin.set(n, q->get(sliceidx[slicenumber],j,k,m,n) ); }
	      
	      OutSliceCart4Inline::slice_q_local_3(kmax2d,qin,qout);
	      
	      for (int n=1; n<=kmax2d; n++)
		{ qslice->set(j,k,m,n, qout.get(n) ); }
	    }
      break;

    default: 
      // should never get here
      eprintf("Error in SliceQ, slicekdind must be 1,2, or 3.  slicekind = %d\n",slicekind[slicenumber]);
    }

*/

}
