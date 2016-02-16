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
#ifndef _GH_COMMON_
#define _GH_COMMON_

#include "common.h"
#ifdef _HAVE_MPI
#include <mpi.h>
#include <fftw3-mpi.h>
#endif //_HAVE_MPI

#define TWOPIPIINVLOGTEN  0.1166503235296796 //ln(10)/(2*pi^2)
#define TWOPIPIINV  0.05066059182116889 //1/(2*pi^2)
#define DZ 0.001
#define NZ 5001
#define NZ_PSOURCES 256

#ifdef _HAVE_MPI
#ifdef _SPREC
#define FLOUBLE_MPI MPI_FLOAT
#else //_SPREC
#define FLOUBLE_MPI MPI_DOUBLE
#endif //_SPREC
#endif //_HAVE_MPI

//Defined in common.c
extern int NodeThis;
extern int NodeLeft;
extern int NodeRight;
extern int NNodes;
extern int IThread0;
extern int MPIThreadsOK;

void mpi_init(int* p_argc,char*** p_argv);
void print_info(char *fmt,...);
void report_error(int level,char *fmt,...);


typedef struct {
  char fnamePk[256];
  double OmegaM;
  double OmegaL;
  double OmegaB;
  double hhub;
  double weos;
  double n_scal;
  double sig8;
  //Derived parameters
  double fgrowth_0;
  double hubble_0;
  double z_max;
  double z_min;
  double r_max;
  double r_min;
  double r2_smooth;
  int do_smoothing;

  //Only used in common.c
  int numk;
  double logkmax;
  double logkmin;
  double idlogk;
  double *logkarr;
  double *pkarr;
  double z_arr_z2r[NZ];
  double r_arr_z2r[NZ];
  double z_arr_r2z[NZ];
  double r_arr_r2z[NZ];
  double growth_d_arr[NZ];
  double growth_v_arr[NZ];
  double glob_idr;

  unsigned int seed_rng;

#ifdef _IRREGULAR_NUTABLE
  char fnameNuTable[256];
  double *nu0_arr;
  double *nuf_arr;
#endif //_IRREGULAR_NUTABLE
  long n_side;
  double nu_max;
  double nu_min;
  int n_nu;

  int n_grid;
  double l_box;
  int nz_here;
  int iz0_here;

  char prefixOut[256];
  double pos_obs[3];

  int do_psources;
  double nz_psources_arr[NZ_PSOURCES];
  double max_Lpdf_arr[NZ_PSOURCES];
  double glob_inv_dz;
  double glob_dz;

  dftw_complex *grid_dens_f;
  flouble *grid_dens;
  dftw_complex *grid_vpot_f;
  flouble *grid_vpot;
  flouble *slice_left;
  flouble *slice_right;
  flouble *grid_rvel;
  double sigma2_gauss;
  int *nsources;
  flouble *maps_HI;
  flouble *maps_PS;
} ParamGetHI;


//////
// Functions defined in cosmo.c
double pk_linear0(ParamGetHI *par,double lgk);
void cosmo_set(ParamGetHI *par);
double r_of_z(ParamGetHI *par,double z);
double z_of_r(ParamGetHI *par,double r);
double dgrowth_of_r(ParamGetHI *par,double r);
double vgrowth_of_r(ParamGetHI *par,double r);


//////
// Functions defined in io_gh.c
ParamGetHI *read_run_params(char *fname);
void write_maps(ParamGetHI *par);
void param_gethi_free(ParamGetHI *par);


//////
// Functions defined in user_defined.c
double fraction_HI(double z);
double bias_HI(double z);


//////
// Functions defined in fourier.c
void init_fftw(ParamGetHI *par);
void create_d_and_vr_fields(ParamGetHI *par);
void end_fftw(void);


//////
// Functions defined in grid_tools.c
void get_HI(ParamGetHI *par);
void get_point_sources(ParamGetHI *par);


//////
// Functions defined in pixelize.c
void mk_T_maps(ParamGetHI *par);
void mk_psources_maps(ParamGetHI *par);


//////
// Functions defined in psources.c
void setup_psources(ParamGetHI *par);
double n_of_z_psources(ParamGetHI *par,double z);
double bias_psources(double z);
double draw_luminosity(ParamGetHI *par,double redshift,gsl_rng *rng);
double temp_of_l(ParamGetHI *par,double L0,double nu_obs,double z,double r,double dOmega);


#endif //_GH_COMMON_
