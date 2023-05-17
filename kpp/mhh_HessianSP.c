/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Hessian Data Structures File                              */
/*                                                                  */
/* Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : mhh_HessianSP.c                           */
/* Time                 : Tue Oct 18 15:45:13 2022                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : mhh.kpp                                   */
/* Output root filename : mhh                                       */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mhh_Parameters.h"
#include "mhh_Global.h"
#include "mhh_Sparse.h"



/* Hessian Sparse Data                                              */
/*                                                                  */
 /* Index i of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_I[] = {
       0,  0,  1,  2,  2,  2,  3,  3,  3,  3,  4,  4,
       5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,
       7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
       9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
      10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12,
      12, 12, 12, 12, 12, 12, 13, 13, 13, 13 }; 

 /* Index j of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_J[] = {
       0, 10,  7,  2,  5,  6,  3,  6,  7, 10,  4, 10,
       5,  5,  5,  4,  5,  5,  6,  6, 10, 11, 12, 12,
       7,  7,  7,  8, 10, 11, 11, 12,  5,  7,  8,  8,
       8,  0,  2,  3,  4,  5,  5,  6,  7,  8,  8,  9,
      10,  0,  2,  5,  6,  8,  8,  9, 10, 10, 10, 10,
      11, 12, 12,  3,  5,  6,  7,  7, 10, 11, 11,  4,
       5,  5, 10, 11, 12, 12,  8, 10, 11, 12 }; 

 /* Index k of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_K[] = {
       9, 10, 11,  9,  8,  9,  9, 11,  9, 11,  9, 12,
       8,  9, 11,  9,  8,  9,  9, 11, 12, 12, 12, 13,
       8,  9, 11, 13, 13, 12, 13, 13,  8,  8,  9, 10,
      13,  9,  9,  9,  9,  8,  9,  9,  9,  9, 10, 10,
      13,  9,  9,  8,  9,  9, 10, 10, 10, 11, 12, 13,
      12, 12, 13,  9, 11, 11,  8, 11, 11, 12, 13,  9,
       8,  9, 12, 12, 12, 13, 13, 13, 13, 13 }; 

