AC_INIT(polyclip,[1.3-2])
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
dnl Check for availability of 64-bit integer types in C++
AC_PROG_CXX
AC_LANG([C++])
dnl Check for presence of C++11 
dnl ................................
dnl I can't get the following to work, even after reinstalling autoconf/m4
dnl AX_CXX_COMPILE_STDCXX_11([noext],[optional])
dnl if test ${HAVE_CXX11}; then
dnl    echo "C++11 present"
dnl else
dnl    echo "C++11 absent"
dnl fi
dnl ................................
dnl Signed 64-bit integer 
long64=""
name64="signed 64-bit integers (cInt)"
dnl if test "${HAVE_CXX11}"; then
dnl    POLYCLIP_LONG64="long long"
dnl else 
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
dnl fi
echo "     In src/clipper.h, ${name64}"
echo "     will be declared as '${long64}'"
dnl Unsigned 64-bit integer 
ulong64=""
uname64="unsigned 64-bit integers (cUInt)"
dnl if test ${HAVE_CXX11}; then
dnl   POLYCLIP_ULONG64="unsigned long long"
dnl else 
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
dnl fi
echo "     In src/clipper.h, ${uname64}"
echo "     will be declared as '${ulong64}'"
dnl Put results in C++ preprocessor flags 
POLYCLIP_CPPFLAGS="-DPOLYCLIP_LONG64=\"${POLYCLIP_LONG64}\" -DPOLYCLIP_ULONG64=\"${POLYCLIP_ULONG64}\""
AC_SUBST(POLYCLIP_CPPFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT


