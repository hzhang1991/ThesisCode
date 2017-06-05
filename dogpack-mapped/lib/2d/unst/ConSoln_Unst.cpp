#include "dogdefs.h"
#include "mesh.h"
#include "DogParams.h"
#include <fstream>

void ConSoln_Unst(const mesh& Mesh, 
        const dTensor3& aux,
        const dTensor3& q, 
        double t, 
        string outputdir)
{
    const int NumElems = q.getsize(1);
    const int     meqn = q.getsize(2);
    const int     kmax = q.getsize(3);
    const int     maux = aux.getsize(2);
    const string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1,write_file2;
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);
    const int NumPhysElems = Mesh.get_NumPhysElems();

    if (t==0) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    if (dogParams.get_mcapa()<1) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m,0.0);

            for (int i=1; i<=NumPhysElems; i++)
            {             
                double dtmp = Mesh.get_area_prim(i);
                double qtmp = q.get(i,m,1);

                qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (fabs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;

    write_file1.close();

}
