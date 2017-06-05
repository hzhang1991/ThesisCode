#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

void RunStartScript(const char* outputdir, int ndims)
{
    char command_str[1024];

    int numchars;
    int exit_status;

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
//  numchars = snprintf(command_str,1024,
//          "if test -f startscript && test -x startscript;\n"
//          "then ./startscript %s %d\n"
//          "else ${DOGPACK}/scripts/startscript %s %d\n"
//          "fi", get_outputdir(),ndims,get_outputdir(),ndims);
    numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %d\n"
            "else ${DOGPACK}/scripts/startscript %s %d\n"
            "fi", outputdir,ndims,outputdir,ndims);
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    exit_status = system(command_str);
}
