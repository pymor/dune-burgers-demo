# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/dune_burgers.cc])
AC_CONFIG_HEADERS([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-pdelab
# this implies checking for [dune-common], [dune-geometry], [dune-grid], [dune-istl], [dune-localfunctions], [dune-pdelab]
DUNE_CHECK_ALL

AX_BOOST_PYTHON()

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/burgers/Makefile
  m4/Makefile
  dune-burgers.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
