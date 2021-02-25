/// @file suffix_decl_hand.h

/* @source suffix_handler.c
** beta 01
** May 26th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Description
** Description
**
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @return something
*******************************************************************************/
#ifndef SUF_DECL_HAND
#define SUF_DECL_HAND

#include "asl.h"
/**
 * Declares the array of Suffixes
 * @param suf_ptr the array of SufDecl type from ASL.
 * @param ptr_name Suffix names.
 * @param rhs_name RHS names
 * @param n_reg_suf number of suffixes for RHS.
 * @param n_rhs number of rhs
 */
void suffix_decl_hand(SufDecl *suf_ptr, char **ptr_name, char **rhs_name, unsigned n_reg_suf, int n_rhs);

#endif /* SUF_DECL_HAND */
