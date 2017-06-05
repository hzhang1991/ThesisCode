#include "IniDocument.h"
#include "EulerParams.h"

EulerParams eulerParams;

void InitApp(IniDocument& ini_doc)
{
  eulerParams.init(ini_doc);
}
