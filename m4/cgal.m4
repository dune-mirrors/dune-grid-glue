dnl This test searches for CGAL in the standard places.
AC_DEFUN([ACX_CHECK_CGAL],
[

AC_ARG_WITH(cgal,
    AC_HELP_STRING([--with-cgal=PATH],[directory where cgal is installed]))
    
acx_cgal_found=no

AC_LANG_PUSH([C++])

dnl if test "$acx_cgal_found" == no; then
	AC_CHECK_HEADER(CGAL/Exact_predicates_inexact_constructions_kernel.h, cgal_have_header=yes, cgal_have_header=no)
	if test "$cgal_have_header" == yes; then
		AC_CHECK_LIB(CGAL, main, cgal_have_lib=yes, cgal_have_lib=no)
		if test "$cgal_have_lib" == no; then
			save_LIBS="$LIBS"; LIBS="$LIBS -lgmp -lmpfr -lm"
        		AC_CHECK_LIB(CGAL, main, [CGAL_LIBS="-lCGAL -lgmp -lmpfr"
						  cgal_have_lib=yes], cgal_have_lib=no)
        		LIBS="$save_LIBS"
		else
			CGAL_LIBS="-lCGAL"
			AC_CHECK_LIB(mpfr, main, [CGAL_LIBS="$CGAL_LIBS -lmpfr"])
			AC_CHECK_LIB(gmp, main, [CGAL_LIBS="$CGAL_LIBS -lgmp"])
			AC_CHECK_LIB(gmpxx, main, [CGAL_LIBS="$CGAL_LIBS -lgmpxx"])
		fi

		if test "$cgal_have_lib" == yes; then 
			acx_cgal_found=yes
		fi
	fi 
	if test "$acx_cgal_found" == yes; then 
		AC_CHECK_LIB(Core, main, [CGAL_LIBS="$CGAL_LIBS -lCore"])
	fi
dnl fi

AC_LANG_POP([C++])

# Store LIBS
AC_SUBST(CGAL_LIBS)

# Define the preprocessor macro HAVE_CGAL in config.h if cgal has been found
if test "$acx_cgal_found" == yes; then
        AC_DEFINE(HAVE_CGAL, 1, [Define to 1 if the cgal library is found])
fi

# Mention in the module summary that cgal has been found
DUNE_ADD_SUMMARY_ENTRY([cgal],[$acx_cgal_found])

AC_MSG_CHECKING(CGAL)
if test "$acx_cgal_found" == yes; then 
	AC_MSG_RESULT(yes);
	$1
else
	AC_MSG_RESULT(no)
	$2
fi])


dnl CHECK CGAL END
