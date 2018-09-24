
#ifndef ORBIT_TOOLS_DEFINES_H_
#define ORBIT_TOOLS_DEFINES_H_

#include <qd/globals.h>
#include <cstdint>


#undef double
typedef double double_type;

#ifndef ORBIT_TOOLS_NAMESPACE_BEGIN
#define ORBIT_TOOLS_NAMESPACE_BEGIN namespace Zeptomoby { namespace OrbitTools {
#endif

#ifndef ORBIT_TOOLS_NAMESPACE_END
#define ORBIT_TOOLS_NAMESPACE_END }}
#endif

#include "orbitTools/SysDefine.h"

#endif
