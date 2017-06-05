#include "dogdefs.h"
#include "ext_time.h"
#include "DogParams.h"
#include "IniDocument.h"

/* This version copied from lib/2d/main_global.cpp */

//
// All hybrid applications call this code.
//
// The purpose of this code is to parse the dogParams section of the ini file,
// and then call the correct RunDogpack* routine buried deeper in the library
// section.  The current options are RunDogpack_Unst or TODO (for the Hybrid
// code).
//
int main_global_hybrid(int argc, char* argv[])
{

    // Open parameters.ini file
    ini_doc.initFromFile("parameters.ini");
    IniDocument::Section& ini_sec = ini_doc["dogParams"];

    // Find out if running in "Cartesian" or "Unstructured" mode
    string mesh_type = ini_sec["mesh_type"];
    void check_mesh_type( string mesh_type );
    check_mesh_type( mesh_type );

    // Determine the output directory, create it if it doesn't exist and parse
    // the dogParams section of the parameters.ini file.
    dogParams.parse_arguments(argc,argv);
    void RunStartScript(const char* outputdir, int ndims);
    RunStartScript(dogParams.get_outputdir(), 2);

    // Call the ``RunDogpack'' routine, which executes the
    // discontinuous Galerkin code
    int m = -1;
    if( mesh_type == "Hybrid" )
    {
        int RunDogpackHybrid( string outputdir );
        m = RunDogpackHybrid( dogParams.get_outputdir() );
    }
//  else
//  {
//      printf("Mesh needs to be hybrid to use this portion of the code\n");
//      int RunDogpack_Unst(string outputdir);
//      m = RunDogpack_Unst( dogParams.get_outputdir() );
//  }
    return m;

}

// quick little error check (repeated from main_global)
void check_mesh_type( string mesh_type )
{

    if (mesh_type != "Hybrid" && mesh_type != "Unstructured")
    {
        printf(
                "\n -------------------------------------------------------- "
                "\n ERROR in parameters.ini "
                "\n -------------------------------------------------------- "
                "\n Mesh type %s is not currently supported."
                "\n Must first set mesh_type in dogpack.data to either "
                "\n     (1) Hybrid ; or "
                "\n     (2) Unstructured. "
                "\n No other mesh type is currently supported. "
                "\n -------------------------------------------------------- \n\n",
                mesh_type.c_str());
        exit(1);
    }
    return;
}
