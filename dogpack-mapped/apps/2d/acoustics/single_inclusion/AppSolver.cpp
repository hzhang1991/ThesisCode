#include "IniDocument.h"
#include "AcousticParams.h"
#include "AppSolver.h"

AcousticParams acousticParams;

void AppSolverCart2::initParams()
{

    DogSolverCart2::initParams();
    acousticParams.init(ini_doc);

}
