#include <stdlib.h>
#include <sstream>
using namespace std;

// TODO - merge this with what's written in the main branch (i.e. use C-strings
// in place of C++-strings
void RunMeshCopyScript(const string& outputdir)
{
  ostringstream command;
  command
    << "${DOGPACK}/scripts/meshcopy_script " << outputdir;
  
  int exit_status = system(command.str().c_str());  
}

