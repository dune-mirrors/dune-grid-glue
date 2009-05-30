// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    conceptchecking.hh
 *  Version:     1.0
 *  Created on:  Jan 25, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: for preprocessor based concept checking
 *  subversion:  $Id$
 *
 */
/**
 * @file conceptchecking.hh
 * @brief for preprocessor based concept checking
 */

/**
 * A concept check class' name following a uniform naming scheme might be "<concept_name>Concept".
 * Such a class is supposed to contain a method "void constraints()" in which all functionality is tested, i.e. calls to
 * all members of the concpet's interface are supposed to be made. Variables needed for this appear as data members of
 * the concept class.
 *
 * This macro may be used in the class body of a concept. It defines a function pointer that points to the concept check class'
 * "constraints()" method thus ensuring that all required functionality is provided by the concept without any runtime overhead.
 *
 * @param type_var the type to be checked, if this is a template
 *                 it is better to introduce another type via typedef
 * @param concept the concept class of which type_var maybe is a model
 */
#define CLASS_REQUIRE(type_var, concept)                            \
  typedef void (concept<type_var>::* func ## type_var ## concept)();  \
                                                                    \
  template<func ## type_var ## concept FuncPtr>                       \
  struct dummy_struct_ ## type_var ## concept {};                     \
                                                                    \
  typedef dummy_struct_ ## type_var ## concept<&concept<type_var>::constraints> dummy_typdef_ ## type_var ## concept;
