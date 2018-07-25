
#ifndef ORBIT_TOOLS_DEFINES_H_
#define ORBIT_TOOLS_DEFINES_H_

#include <qd/globals.h>
#include <cstdint>

#ifndef QD_INTEGRATION_ENABLED
#ifdef ORBIT_TOOLS_ENABLE_QD_QD_INTEGRATION
namespace qd
{
    static inline qd_real to_real(double d) { return d; }
}

#define QD_INTEGRATION_ENABLED 1

#elif defined(ORBIT_TOOLS_ENABLE_QD_DD_INTEGRATION)
namespace qd
{
    static inline dd_real to_real(double d) { return d; }
}

#define QD_INTEGRATION_ENABLED 1

#else

static inline int to_double(int i) { return i; }
static inline int64_t to_double(int64_t i) { return i; }
static inline double to_double(double d) { return d; }

namespace qd
{
    static inline double to_real(double d) { return d; }
}

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
