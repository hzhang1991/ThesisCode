#include <stdlib.h>
#include "assert.h"

void RunMeshCopyScript(const char* outputdir)
{

    char command_str[1024];

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "${DOGPACK}/scripts/meshcopy_script %s\n",
            outputdir);
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str); 

}
