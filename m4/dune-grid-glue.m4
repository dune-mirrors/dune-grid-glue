# Additional checks needed to build the module
AC_DEFUN([DUNE_GRID_GLUE_CHECKS],[
    AC_REQUIRE([ACX_CHECK_CGAL])
    AC_REQUIRE([DUNE_PATH_PSURFACE])
])


# Additional checks needed to find the module
AC_DEFUN([DUNE_GRID_GLUE_CHECK_MODULE],
[
    DUNE_CHECK_MODULES([dune-grid-glue], [grid-glue/merging/cgalmerge.hh], [[
      &CGALMerge<1,double>::build;]])
])
