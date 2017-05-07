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
#include "common_jt.h"
#include <chealpix.h>

#define PREFAC_TNU 3.6E9 // (1 h)*(1 MHz)
#define PREFAC_FWHM 1030610.55 //c/(1 m * 1 MHz) * (radians to arcmin)

static void write_arrays(ParamsJoinT *pars)
{
  FILE *fil;
  char fname[256];
  int ii;

  sprintf(fname,"%s_instrument.dat",pars->prefix_out);
  fil=fopen(fname,"w");
  if(fil==NULL) error_open_file(fname);
  for(ii=0;ii<pars->n_nu;ii++) {
    double nu=pars->nutable[0][ii];
    double dnu=pars->nutable[2][ii]-pars->nutable[1][ii];
    double fwhm_rad;
    double n_pointings,sigma2_noise_pointing;
    double n_pixels,sigma2_noise_pixel;
    n_pixels=pars->fsky*nside2npix(pars->nside);
    if(pars->dish_diameter>0) {
      double fwhm=PREFAC_FWHM/(nu*pars->dish_diameter);
      fwhm_rad=fwhm*M_PI/(180*60);
      n_pointings=(pars->fsky)*4*M_PI/(1.1330900354567985*fwhm_rad*fwhm_rad);
    }
    else {
      fwhm_rad=0;
      n_pointings=pars->fsky*nside2npix(pars->nside);
    }

    if(pars->tsys>0) {
      sigma2_noise_pointing=n_pointings*(pars->tsys*pars->tsys)/
	(PREFAC_TNU*dnu*pars->time_tot*pars->n_dish);
      sigma2_noise_pixel=n_pixels*(pars->tsys*pars->tsys)/
	(PREFAC_TNU*dnu*pars->time_tot*pars->n_dish);
    }
    else {
      sigma2_noise_pointing=0;
      sigma2_noise_pixel=0;
    }

    fprintf(fil,"%lE %lE %lE %lE\n",nu,180*fwhm_rad/M_PI,
	    sqrt(sigma2_noise_pointing),
	    sqrt(sigma2_noise_pixel));
  }
  fclose(fil);
}

static void assert_nside(long nside1,long nside2)
{
  if(nside1!=nside2) {
    fprintf(stderr,"CRIME: read nside is not the same as provided : %ld != %ld\n",
	    nside1,nside2);
    exit(1);
  }
}

static void add_maps(long npix,flouble *map_to_sum,flouble *map_result)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(npix,map_to_sum,map_result)
#endif //_HAVE_OMP
  {
    long ii;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<npix;ii++)
      map_result[ii]+=map_to_sum[ii];
  }
}

void merge_maps(ParamsJoinT *pars,flouble **map_out)
{
  int inu;
  long nside_in,npix_in;
  char fname[256];
  int add_cosmo=1;
  int add_egfree=1;
  int add_gfree=1;
  int add_psources=1;
  int add_gsync=1;
  int add_custom=1;
 
  if(!strcmp(pars->prefix_cosmo,"none")) add_cosmo=0;
  if(!strcmp(pars->prefix_extragalactic_freefree,"none")) add_egfree=0;
  if(!strcmp(pars->prefix_galactic_freefree,"none")) add_gfree=0;
  if(!strcmp(pars->prefix_point_sources,"none")) add_psources=0;
  if(!strcmp(pars->prefix_synchrotron,"none")) add_gsync=0;
  if(!strcmp(pars->prefix_custom,"none")) add_custom=0;
  
  if(add_cosmo)
    sprintf(fname,"%s_001.fits",pars->prefix_cosmo);
  else if(add_egfree)
    sprintf(fname,"%s_001.fits",pars->prefix_extragalactic_freefree);
  else if(add_gfree)
    sprintf(fname,"%s_001.fits",pars->prefix_galactic_freefree);
  else if(add_psources)
    sprintf(fname,"%s_001.fits",pars->prefix_point_sources);
  else if(add_gsync)
    sprintf(fname,"%s_001.fits",pars->prefix_synchrotron);
  else if(add_custom)
    sprintf(fname,"%s_001.fits",pars->prefix_custom);
  else return;

  flouble *dum=he_read_healpix_map(fname,&nside_in,0);
  free(dum);
  npix_in=nside2npix(nside_in);

  for(inu=0;inu<pars->n_nu;inu++) {
    long nside_test;
    flouble *map_read;
    flouble *map_in=(flouble *)calloc(npix_in,sizeof(flouble));

#ifdef _VERBOSE    
    printf("%d\n",inu);
#endif //_VERBOSE    

    //Cosmological signal
    if(add_cosmo) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_cosmo,inu+1);
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    //Extragalactic free-free
    if(add_egfree) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_extragalactic_freefree,inu+1);
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    //Galactic free-free
    if(add_gfree) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_galactic_freefree,inu+1);
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    //Point sources
    if(add_psources) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_point_sources,inu+1);
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    //Synchrotron
    if(add_gsync) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_synchrotron,inu+1);
      if((pars->do_polarization)&&(pars->polarization_leakage>0)) {
	long ii;
	
	map_read=he_read_healpix_map(fname,&nside_test,1);
	assert_nside(nside_in,nside_test);
	for(ii=0;ii<npix_in;ii++) map_read[ii]*=pars->polarization_leakage;
	add_maps(npix_in,map_read,map_in);
	free(map_read);
      }
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    //Custom
    if(add_custom) {
      sprintf(fname,"%s_%03d.fits",pars->prefix_custom,inu+1);
      map_read=he_read_healpix_map(fname,&nside_test,0);
      assert_nside(nside_in,nside_test);
      add_maps(npix_in,map_read,map_in);
      free(map_read);
    }

    he_udgrade(map_in,nside_in,map_out[inu],pars->nside,0);

    free(map_in);
  }
}

void apply_beam(ParamsJoinT *pars,flouble **maps)
{
  int inu;
  int lmax=3*(pars->nside)-1;
  fcomplex **alms;

  alms=(fcomplex **)my_malloc(pars->n_nu*sizeof(fcomplex *));
  for(inu=0;inu<pars->n_nu;inu++)
    alms[inu]=(fcomplex *)my_malloc(he_nalms(lmax)*sizeof(fcomplex));

  he_map2alm(pars->nside,lmax,pars->n_nu,maps,alms);

  //alter alms
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(pars,lmax,alms)
#endif //_HAVE_OMP
  {
    int ii;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<pars->n_nu;ii++) {
      double nu=pars->nutable[0][ii];
      double fwhm=PREFAC_FWHM/(nu*pars->dish_diameter);

#ifdef _VERBOSE
      printf("%d\n",ii);
#endif //_VERBOSE
      he_alter_alm(lmax,fwhm,alms[ii],NULL);
    }// end omp for
  } // end omp parallel

  he_alm2map(pars->nside,lmax,pars->n_nu,maps,alms);

  for(inu=0;inu<pars->n_nu;inu++)
    free(alms[inu]);
  free(alms);
}

void apply_mask(ParamsJoinT *pars,flouble **maps)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(pars,maps)
#endif //_HAVE_OMP
  {
    int inu;
    long npix=nside2npix(pars->nside);
    
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(inu=0;inu<pars->n_nu;inu++) {
      long ii;

      for(ii=0;ii<npix;ii++) {
	if(pars->mask[ii]==0) maps[inu][ii]=0;
      }
    }// end omp for
  }// end omp parallel
}

void add_noise(ParamsJoinT *pars,flouble **maps)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(pars,maps)
#endif //_HAVE_OMP
  {
    int inu;
    long npix=nside2npix(pars->nside);
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=pars->seed+ithr;
    gsl_rng *rng=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(inu=0;inu<pars->n_nu;inu++) {
      long ii;
      double dnu=pars->nutable[2][inu]-pars->nutable[1][inu];
      double n_elements=pars->fsky*npix;
#ifdef _ADD_ATM
      double nu=pars->nutable[0][inu];
      double t_sky=6E4*pow(nu/300.,-2.5);
      double t_sys=pars->tsys+t_sky;
#else //_ADD_ATM
      double t_sys=pars->tsys;
#endif //_ADD_ATM     
      double sigma2_noise=n_elements*t_sys*t_sys/(PREFAC_TNU*dnu*pars->time_tot*pars->n_dish);

      for(ii=0;ii<npix;ii++) {
	double mod,phase,noise;
	rng_delta_gauss(&mod,&phase,rng,sigma2_noise);
	noise=mod*cos(phase);

	maps[inu][ii]+=noise;
      }
    }// end omp for

    end_rng(rng);
  }// end omp parallel
}

int main(int argc,char **argv)
{
  char fname_input[256];
  ParamsJoinT *pars;
  flouble **maps_out;
  int ii;

  if(argc!=2) {
    fprintf(stderr,"Usage: ./JoinT file_name\n");
    exit(1);
  }
  sprintf(fname_input,"%s",argv[1]);

  setbuf(stdout,NULL);
  printf("\n");
  printf("|-------------------------------------------------|\n");
  printf("|                      JoinT                      |\n");
  printf("|-------------------------------------------------|\n\n");
  
  timer(4);

  pars=read_input_params_JoinT(fname_input);

  maps_out=(flouble **)my_malloc(pars->n_nu*sizeof(flouble *));
  for(ii=0;ii<pars->n_nu;ii++)
    maps_out[ii]=(flouble *)my_calloc(nside2npix(pars->nside),sizeof(flouble));

  printf("Merging maps\n");
  timer(0);
  merge_maps(pars,maps_out);
  timer(1);
  
  if(pars->dish_diameter>0) {
    printf("Applying beam for a dish of size %.2lf m\n",
	   pars->dish_diameter);
    timer(0);
    apply_beam(pars,maps_out);
    timer(1);
  }

  if(pars->tsys>0) {
    printf("Adding Gaussian noise\n");
    timer(0);
    add_noise(pars,maps_out);
    timer(1);
  }

  printf("Applying mask\n");
  timer(0);
  apply_mask(pars,maps_out);
  timer(1);
  
  write_maps(pars->prefix_out,pars->nside,pars->n_nu,maps_out);
  write_arrays(pars);

  for(ii=0;ii<pars->n_nu;ii++)
    free(maps_out[ii]);
  free(maps_out);
  
  free_params_JoinT(pars);
  timer(5);

  printf("\n");
  printf("|-------------------------------------------------|\n\n");

  return 0;
}
