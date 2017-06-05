#include <stdio.h>
#include <stdlib.h>

int RunDogpack_Unst(const char* outputdir)
{

    // Output title information
    printf( "\n"
            "\n  Error: mesh_type = Unstructured in parameters.ini,"
            "\n  but COMPILE_WITH_UNST was not set when linking.  Put the line"
            "\n"
            "\n    COMPILE_WITH_UNST = 1"
            "\n"
            "\n  in the Makefile prior to including Makefile.defs"
            "\n  from the 2d library, rm dog.exe, and recompile."
            "\n\n");
    exit(1);

}
