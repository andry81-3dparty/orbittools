
#ifndef ORBIT_TOOLS_DEFINES_H_
#define ORBIT_TOOLS_DEFINES_H_

#ifndef QD_INTEGRATION_ENABLED
#ifdef ORBIT_TOOLS_ENABLE_QD_QD_INTEGRATION
#include <qd/qd_real.h>

static inline qd_real from_double(double d) { return d; }

#define QD_INTEGRATION_ENABLED 1

#elif defined(ORBIT_TOOLS_ENABLE_QD_DD_INTEGRATION)
#include <qd/dd_real.h>

static inline dd_real from_double(double d) { return d; }

#define QD_INTEGRATION_ENABLED 1

#else
static inline double to_double(double d) { return d; }
static inline double from_double(double d) { return 0; }

#define QD_INTEGRATION_ENABLED 0

#endif
#endif

#ifndef ORBIT_TOOLS_NAMESPACE_BEGIN
#define ORBIT_TOOLS_NAMESPACE_BEGIN namespace Zeptomoby { namespace OrbitTools {
#endif

#ifndef ORBIT_TOOLS_NAMESPACE_END
#define ORBIT_TOOLS_NAMESPACE_END }}
#endif

#include "orbitTools/SysDefine.h"

#endif
