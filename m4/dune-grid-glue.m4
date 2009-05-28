# Additional checks needed to build the module
AC_DEFUN([DUNE_GRID_GLUE_CHECKS])
# Additional checks needed to find the module
AC_DEFUN([DUNE_GRID_GLUE_CHECK_MODULE],
[
    DUNE_CHECK_MODULES([dune-grid-glue], [glue/misc/geometry.hh])    
])
