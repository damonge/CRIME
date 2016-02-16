///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CRIME.                                       //
//                                                                   //
// CRIME is free software: you can redistribute it and/or modify it  //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CRIME is distributed in the hope that it will be useful, but      //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CRIME.  If not, see <http://www.gnu.org/licenses/>.    //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef _JOINT_COMMON_
#define _JOINT_COMMON_

#include "common.h"

typedef struct {
  char fname_input[256];
  char prefix_cosmo[256];
  char prefix_synchrotron[256];
  char prefix_galactic_freefree[256];
  char prefix_extragalactic_freefree[256];
  char prefix_point_sources[256];
  char prefix_custom[256];
  char prefix_out[256];
  char fname_nutable[256];
  char fname_mask[256];
  unsigned int seed;
  int n_nu;
  double **nutable;
  int nside;
  int do_polarization;
  double polarization_leakage;
  double dish_diameter;
  int n_dish;
  double tsys;
  double time_tot;
  int *mask;
  double fsky;
} ParamsJoinT;


//////
// Functions defined in common_jt.c
ParamsJoinT *set_default_params_JoinT(void);

void free_params_JoinT(ParamsJoinT *pars);

//////
// Functions defined in io_jt.c
void write_maps(char *prefix_out,int nside,int n_nu,flouble **maps);

ParamsJoinT *read_input_params_JoinT(char *fname);

double **read_nutable(char *fname,int *n_nu);

#endif //_JOINT_COMMON_
