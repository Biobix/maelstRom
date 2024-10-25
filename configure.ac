# autoconf this file --> configure
AC_INIT
AC_CONFIG_SRCDIR([MAGE])
AC_CONFIG_SRCDIR([DESCRIPTION])

dnl Select an optional include path, from a configure option
dnl or from an environment variable  >> mpfr_INCLUDE <<
AC_ARG_WITH([mpfr-include],
            AS_HELP_STRING([--with-mpfr-include=INCLUDE_PATH],[the location of MPFR header files]),
            [mpfr_include_path=$withval])

if test -n "$mpfr_include_path" ; then
   MPFR_CPPFLAGS="-I${mpfr_include_path}"
elif test -n "${mpfr_INCLUDE}" ; then
   MPFR_CPPFLAGS="-I${mpfr_INCLUDE}"
fi

dnl Ditto for a library path  environment variable  >> mpfr_LDFLAGS <<
AC_ARG_WITH([mpfr-lib],
            AS_HELP_STRING([--with-mpfr-lib=LIB_PATH],[the location of MPFR libraries]),
            [mpfr_lib_path=$withval])

if test -n "$mpfr_lib_path" ; then
   MPFR_LDFLAGS="-L$mpfr_lib_path ${LDFLAGS}"
elif test -n "${mpfr_LDFLAGS}" ; then
   MPFR_LDFLAGS="-L${mpfr_LDFLAGS} ${LDFLAGS}"
fi

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`


dnl Setup the compilers
AC_PROG_CC
AC_PROG_CPP

dnl use the include and lib
CPPFLAGS="${CPPFLAGS} ${MPFR_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${MPFR_LDFLAGS}"

PKG_CPPFLAGS=$CPPFLAGS
PKG_CFLAGS=$CFLAGS
PKG_LDFLAGS=$LDFLAGS

AC_CHECK_HEADER([mpfr.h], ,
    [AC_MSG_ERROR([Header file mpfr.h not found; maybe use --with-mpfr-include=INCLUDE_PATH])])

AC_CHECK_HEADER([gmp.h], ,
    [AC_MSG_ERROR([Header file gmp.h not found; maybe use --with-mpfr-include=INCLUDE_PATH])])

AC_CHECK_LIB(gmp, __gmpz_init, ,
    [AC_MSG_ERROR([GNU MP not found, see README])])

AC_CHECK_LIB(mpfr, mpfr_init, ,
    [AC_MSG_ERROR([MPFR Library not found, see README])])

AC_CHECK_LIB(mpfr, mpfr_digamma, ,
    [AC_MSG_ERROR([MPFR Library must be at least version 3.0.0, see README])])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LDFLAGS)



## Use gsl-config to find arguments for
## compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, lets use it
if test "${GSL_CONFIG}" != ""; then
  # Use gsl-config for header and linker args
  # (without BLAS which we get from R)
  GSL_CFLAGS=`${GSL_CONFIG} --cflags`
  GSL_LIBS=`${GSL_CONFIG} --libs-without-cblas`
else
  AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)



AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT
