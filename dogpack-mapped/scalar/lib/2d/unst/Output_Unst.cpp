#include "dogdefs.h"
#include "DogParams.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"
#include <sstream>
#include <string>

// Output solution on unstructured mesh
void Output_Unst(const mesh& Mesh, const dTensor3& aux,
        const dTensor3& q, double t, int nframe, 
        string outputdir)
{ 
    int i,k,m;
    int NumElems = q.getsize(1);
    int     meqn = q.getsize(2);
    int     kmax = q.getsize(3);
    int     maux = aux.getsize(2);
    int  method1 = int((sqrt(1+8*kmax)-1)/2);
    assert(method1==dogParams.get_space_order());  
    void WriteOutput_Unst(string fname, string varname, 
            const dTensor3& q, double t);

    // Write output to files
    ostringstream basename;
    basename << setfill('0') << setw(4) << nframe;
    string fname1  = outputdir+"/q"+basename.str();
    string fname2  = outputdir+"/a"+basename.str();

    // Output Q values
    WriteOutput_Unst(fname1, "q", q, t);

    // Output aux values
    if (maux>0) 
    {  WriteOutput_Unst(fname2, "a", aux, t);  }
}

void WriteOutput_Unst(string fname, string varname, 
        const dTensor3& q, double t)
{
    int i,m,k;
    int NumElems = q.getsize(1);
    int     meqn = q.getsize(2);
    int     kmax = q.getsize(3);
    fname+=".dat";

    ofstream write_file(fname.c_str(), ios::out);
    write_file << setprecision(16);
    write_file << setw(24) << scientific << t << endl;

    for (k=1; k<=kmax; k++)
    for (m=1; m<=meqn; m++)      
    for (i=1; i<=NumElems; i++)
    {
        double tmp = q.get(i,m,k);
        write_file << setw(24) << scientific << tmp << endl;
    }
    write_file.close();
}
