#ifndef _GH_COMMON_
#define _GH_COMMON_

#include "common.h"

#define TWOPIPIINVLOGTEN  0.1166503235296796 //ln(10)/(2*pi^2)
#define N_SUBPART_HI 10

typedef struct {
  char fnamePk[256]; //
  double OmegaM; //
  double OmegaL; //
  double OmegaB; //
  double OmegaK; //
  double hhub; //
  double w0; //
  double wa; //
  double n_scal; //
  double sig8; //
  double growth_0; //
  double fgrowth_0; //
  double hubble_0; //
  double rmin;
  double rmax;

  SplPar *spl_rz; //
  SplPar *spl_dgf; //
  SplPar *spl_vgf; //
  SplPar *spl_pk; //
  double *r_arr;
  double *z_arr;
  double idr;

  char fname_nulist[256]; //
  ParamMaps *pmap; //

  int n_grid; //
  int nz_here; //
  int iz0_here; //
  double l_box; //

  char prefix_out[256]; //
  unsigned int seed_rng; //

  double pos_obs[3]; //

  int do_psources; //

  dftw_complex *grid_dens_f;
  flouble *grid_dens;
  dftw_complex *grid_vpot_f;
  flouble *grid_vpot;
  flouble *slice_left;
  flouble *slice_right;
  flouble *grid_rvel;
  flouble *maps_HI;
  flouble *maps_PS;
} ParamGetHI;

//Defined in io_gh.c
ParamGetHI *read_run_params(char *fname);
void param_gethi_free(ParamGetHI *par);
void write_output(ParamGetHI *par);

//Defined in cosmo.c
double pk_linear0(ParamGetHI *par,double lgk);
double z_of_r(ParamGetHI *par,double r);
double r_of_z(ParamGetHI *par,double z);
double dgrowth_of_z(ParamGetHI *par,double z);
double vgrowth_of_z(ParamGetHI *par,double r);
void cosmo_set(ParamGetHI *par);

//Defined in fourier.c
void init_fftw(ParamGetHI *par);
void create_density_and_radvel(ParamGetHI *par);
void end_fftw(void);

//Defined in grid2maps.c
void mk_HI_maps(ParamGetHI *par);

//Defined in user_defined.c
double fraction_HI(double z);
double bias_HI(double z);


#endif //_GH_COMMON_
