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
#ifndef _FORGET_COMMON_
#define _FORGET_COMMON_

#include "common.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define NSIDE_TEMPLATE 512
#define NSIDE_POL 128
#define NU_HASLAM 408.0
#define SPECIN_OFFSET 0.0

typedef struct {
  char fname_input[256];
  char fname_nutable[256];
  char fname_haslam[256];
  char fname_specin[256];
  char prefix_out[256];
  char cl_model[256];
  int lmin;
  int lmax;
  int nside;
  unsigned int seed;
  int do_galaxy;
  double amp;
  double beta;
  double alpha;
  double xi;
  double nu_ref;
  int lref;
  int n_nu;
  double **nutable;
  int do_polarization;
  double xi_polarization;
  double beta_polarization;
  char fname_faraday[256];
} ParamsForGet;


//////
// Functions defined in common_fg.c
ParamsForGet *set_default_params_ForGet(void);

void free_params_ForGet(ParamsForGet *pars);


//////
// Functions defined in sck_maps.c
flouble **get_sck_maps(ParamsForGet *pars);


//////
// Functions defined in gsync_i.c
void add_haslam(ParamsForGet *pars,flouble **maps);

void constrain_alms(ParamsForGet *pars,gsl_matrix *eigenvec,
		    gsl_vector *eigenval,fcomplex **alms);

void init_galaxy_template(ParamsForGet *pars);

void end_galaxy_template(void);


//////
// Functions defined in gsync_qu.c
void do_polarization(ParamsForGet *pars,flouble ***maps_nuspace_r,
		     flouble ***maps_nuspace_i);


//////
// Functions defined in io_fg.c
void write_maps(char *prefix_out,int nside,int n_nu,
		flouble **maps_i,flouble **maps_q,
		flouble **maps_u);

ParamsForGet *read_input_params_ForGet(char *fname);

double **read_nutable(char *fname,int *n_nu);


#endif //_FORGET_COMMON_
