* README_EN.txt
* 2018.07.16
* orbittools

1. DESCRIPTION
2. LICENSE
3. REPOSITORIES
4. INSTALLATION
5. AUTHOR EMAIL

-------------------------------------------------------------------------------
1. DESCRIPTION
-------------------------------------------------------------------------------
OrbitTools patched sources fork from: http://www.zeptomoby.com/satellites/

The OrbitTools Libraries
NORAD SGP4/SDP4 Implementations in C++ and C#
by
Michael F. Henry

The OrbitTools Libraries are implementations of NORAD algorithms for determining
satellite location and velocity in Earth orbit. The algorithms come from the
December, 1980 NORAD document "Space Track Report No. 3". The orbital algorithms
implemented in OrbitTools are: SGP4, for "near-Earth" objects, and SDP4 for
"deep space" objects. These algorithms are widely used in the satellite tracking
community and produce very accurate results when provided with current NORAD
two-line element data.

The original library patched to fix these issue:

1. All double's replaced by dd_real from QD library to decrease precision
   fluctuation to the minimum in certain cases.

2. time_t in cJulian class replaced by floating point value to avoid truncation
   to seconds.

All patches improved precision from ~400 meters per 100km altitude along
velocity vector in certain routines up to 10^-7 meters per 100km altitude along
velocity vector.

Cmake scripts uses the cmake modules from the tacklelib library:

https://svn.code.sf.net/p/tacklelib/cmake/trunk

Third-party QD library hosts here:

https://svn.code.sf.net/p/orbittools/qd_/trunk

-------------------------------------------------------------------------------
2. LICENSE
-------------------------------------------------------------------------------
The MIT license (see included text file "license.txt" or
https://en.wikipedia.org/wiki/MIT_License)

-------------------------------------------------------------------------------
3. REPOSITORIES
-------------------------------------------------------------------------------
Primary:
  * https://svn.code.sf.net/p/orbittools/orbittools/trunk
First mirror:
  * https://github.com/andry81/orbittools--orbittools.git
Second mirror:
  * https://bitbucket.org/andry81/orbittools-orbittools.git

-------------------------------------------------------------------------------
4. INSTALLATION
-------------------------------------------------------------------------------
N/A

-------------------------------------------------------------------------------
5. AUTHOR EMAIL
-------------------------------------------------------------------------------
Andrey Dibrov (andry at inbox dot ru)
