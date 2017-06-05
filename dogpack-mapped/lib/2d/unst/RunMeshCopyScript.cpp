#include <stdlib.h>
#include <sstream>
using namespace std;

void RunMeshCopyScript(const string& outputdir)
{
  ostringstream command;
  command
    << "${DOGPACK}/scripts/meshcopy_script " << outputdir;
  
  int exit_status = system(command.str().c_str());  
}

