#include "IniDocument.h"
#include "EulerParams.h"
#include "AppSolver.h"

EulerParams eulerParams;

void AppSolverCart2::initParams()
{
    eulerParams.init(ini_doc);
}
