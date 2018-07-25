AC_INIT(polyclip,[1.9-1])
dnl Determine version of R
CXX11FUN=`"${R_HOME}/bin/Rscript" -e 'if(getRversion() < "3.4") cat("CXX1X", fill=TRUE) else cat("CXX11", fill=TRUE)'`
dnl Check for presence of C++11 
CXX1X=`"${R_HOME}/bin/R" CMD config "${CXX11FUN}"`
if test "${CXX1X}" != ""; then
   echo "compiling under C++11"
   CXX1XSTD=`"${R_HOME}/bin/R" CMD config "${CXX11FUN}STD"`
   CXX="${CXX1X} ${CXX1XSTD}"
   POLYCLIP_CXX_DECLAR="CXX_STD = CXX11"
   CXXFLAGS=`"${R_HOME}/bin/R" CMD config "${CXX11FUN}FLAGS"`
else
   echo "C++11 is not supported; compiling under vanilla C++"
   CXX=`"${R_HOME}/bin/R" CMD config CXX`
   POLYCLIP_CXX_DECLAR=""
   CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
fi
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CXX
AC_LANG([C++])
dnl Check for pkgconf or pkg-config
if test -z "$PKG_CONFIG"; then
   if command -v pkgconf 2>&1 > /dev/null; then
      PKG_CONFIG=pkgconf
   else
      if command -v pkg-config 2>&1 > /dev/null; then
         PKG_CONFIG=pkg-config
      fi
   fi
fi
echo "Using PKG_CONFIG: $PKG_CONFIG"
dnl Check if polyclipping library is installed on system.
if test -n "$PKG_CONFIG" && $PKG_CONFIG --exists polyclipping; then
   echo Compiling against system copy of polyclipping library.
   POLYCLIP_LIBS=`$PKG_CONFIG --libs polyclipping`
   POLYCLIP_CPPFLAGS=`$PKG_CONFIG --cflags polyclipping`
   dnl The polyclip type declarations are needed to compile the R interface
   LLTYPES="-DPOLYCLIP_LONG64=\"signed long long\" -DPOLYCLIP_ULONG64=\"unsigned long long\""
   POLYCLIP_CPPFLAGS="${POLYCLIP_CPPFLAGS} ${LLTYPES}"
else
   echo Compiling against bundled copy of clipper library.
   dnl Check for availability of 64-bit integer types in C++
   dnl Signed 64-bit integer
   long64=""
   name64="signed 64-bit integers (cInt)"
   if test "${CXX1X}" != ""; then
      long64="signed long long"
      POLYCLIP_LONG64="${long64}"
   else
      AC_CHECK_TYPES([int64_t],[long64="int64_t"],[],[#include <stdint.h>])
      if test "${long64}" != ""; then
         POLYCLIP_LONG64="${long64}"
      else
         AC_CHECK_TYPES([int_fast64_t],[long64="int_fast64_t"],[],[#include <stdint.h>])
         if test "${long64}" != ""; then
            POLYCLIP_LONG64="${long64}"
         else
            AC_CHECK_TYPES([int_least64_t],[long64="int_least64_t"],[],[#include <stdint.h>])
            if test "${long64}" != ""; then
               POLYCLIP_LONG64="${long64}"
            else
               AC_CHECK_TYPES([long long],[long64="long long"])
               if test "${long64}" != ""; then
                  POLYCLIP_LONG64="${long64}"
               else
                  echo "Error: unable to find a C++ data type for ${name64}"
                  exit 1
               fi
            fi
         fi
      fi
   fi
   echo "     In the clipper library, ${name64}"
   echo "     will be declared as '${long64}'"
   dnl Unsigned 64-bit integer
   ulong64=""
   uname64="unsigned 64-bit integers (cUInt)"
   if test "${CXX1X}" != ""; then
      ulong64="unsigned long long"
      POLYCLIP_ULONG64="${ulong64}"
   else
      AC_CHECK_TYPES([uint64_t],[ulong64="uint64_t"],[],[#include <stdint.h>])
      if test "${ulong64}" != ""; then
         POLYCLIP_ULONG64="${ulong64}"
      else
         AC_CHECK_TYPES([uint_fast64_t],[ulong64="uint_fast64_t"],[],[#include <stdint.h>])
         if test "${ulong64}" != ""; then
            POLYCLIP_ULONG64="${ulong64}"
         else
            AC_CHECK_TYPES([uint_least64_t],[ulong64="uint_least64_t"],[],[#include <stdint.h>])
            if test "${ulong64}" != ""; then
               POLYCLIP_ULONG64="${ulong64}"
            else
               AC_CHECK_TYPES([unsigned long long],[ulong64="unsigned long long"])
               if test "${ulong64}" != ""; then
                  POLYCLIP_ULONG64="${ulong64}"
               else
                  echo "Error: unable to find a C++ data type for ${uname64}"
                  exit 1
               fi
            fi
         fi
      fi
   fi
   echo "     In the clipper library, ${uname64}"
   echo "     will be declared as '${ulong64}'"
   dnl Put results in C++ preprocessor flags
   POLYCLIP_CPPFLAGS="-DPOLYCLIP_LONG64=\"${POLYCLIP_LONG64}\" -DPOLYCLIP_ULONG64=\"${POLYCLIP_ULONG64}\""
   POLYCLIP_OBJECTS="clipper.o"
fi
POLYCLIP_CXXFLAGS="${CXXFLAGS}"
AC_SUBST(POLYCLIP_CPPFLAGS)
AC_SUBST(POLYCLIP_CXXFLAGS)
AC_SUBST(POLYCLIP_LIBS)
AC_SUBST(POLYCLIP_OBJECTS)
AC_SUBST(POLYCLIP_CXX_DECLAR)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT


