#include "dogdefs.h"

// Used for access to ini_doc
#include "IniDocument.h"

// Why include the Cartesian parameters here? (-DS)
#include "DogSolverCart2.h"
#include "DogParamsCart2.h"

// This function performs the following tasks:
//
// 1.) Parse the Parameters file
//
// 2.) Create the output directory
//
// 3.) Calls the "RunDogpack" routine, which executes the code
//
// Currently the only routines that call this function are the 1D
// applications, and a handful of 2D applications, including many of the
// unstructured solvers.  I assume this means that this method will become
// deprecated in the future? (-DS)
//
int main_global(int argc, char* argv[])
{

    // Open parameters.ini file
    ini_doc.initFromFile("parameters.ini");

    // Read the dogParams "section"
    IniDocument::Section& dogParams_section = ini_doc["dogParams"];

    // Find out if running in "Cartesian" or "Unstructured" mode
    string mesh_type = dogParams_section["mesh_type"];
    if (mesh_type != "Cartesian" && mesh_type != "Unstructured")
    {
        printf(
                "\n -------------------------------------------------------- "
                "\n ERROR in parameters.ini "
                "\n -------------------------------------------------------- "
                "\n Mesh type %s is not currently supported."
                "\n Must first set mesh_type in dogpack.data to either "
                "\n     (1) Cartesian; or "
                "\n     (2) Unstructured. "
                "\n No other mesh type is currently supported. "
                "\n -------------------------------------------------------- \n\n",
                mesh_type.c_str());
        exit(1);
    }

    // determine and create the output directory
    DogSolver::parse_arguments(argc,argv);
    void RunStartScript(int ndims);
    RunStartScript(2);

    // Call the ``RunDogpack'' routine, which executes the
    // discontinuous Galerkin code
    int m;
    int RunDogpack_Unst(string outputdir);
    if (mesh_type == "Cartesian")
    {
        // Cartesian mesh
        DogSolverCart2 dogSolverCart2;
        m = dogSolverCart2.RunDogpack();
    }
    else
    {
        // Unstructured mesh
        m = RunDogpack_Unst(DogSolver::get_outputdir());
    }

    return m;

}
