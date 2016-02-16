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
#include "common_fg.h"
#include <gsl/gsl_blas.h>
#include <chealpix.h>

#define NSIDE_REGIONS 16
#define LMAX_TR_HASLAM 1024
#define LMAX_HASLAM 200
#define NU_HASLAM 408.0
static flouble *galaxy_amplitude;
static flouble *galaxy_template;
static flouble *specin_template;
static int i_haslam;

void add_haslam(ParamsForGet *pars,flouble **maps)
{
#pragma omp parallel default(none) \
  shared(pars,maps,galaxy_template,specin_template)
  {
    int inu;

#pragma omp for
    for(inu=0;inu<pars->n_nu;inu++) {
      long ii;
      double nu=pars->nutable[0][inu];
      for(ii=0;ii<nside2npix(pars->nside);ii++) {
	long ipix_haslam;
	double theta,phi;
	double t_haslam,specin,amplitude;
	double map0=maps[inu][ii];
	
	pix2ang_ring(pars->nside,ii,&theta,&phi);
	ang2pix_ring(NSIDE_TEMPLATE,theta,phi,&ipix_haslam);
	t_haslam=galaxy_template[ipix_haslam];
	specin=specin_template[ipix_haslam];
	//      amplitude=galaxy_amplitude[ipix_haslam];
	amplitude=1;
	
	maps[inu][ii]=amplitude*map0+t_haslam*pow(nu/NU_HASLAM,specin);
      }
    }
  }
}

void constrain_alms(ParamsForGet *pars,gsl_matrix *eigenvec,
		    gsl_vector *eigenval,fcomplex **alms)
{
  int ii,k0=-1;
  double emax=-10000;
  double ev_h_k0;
  gsl_vector *ev_haslam=gsl_vector_alloc(pars->n_nu);

  //Get B(i_haslam,k) and B(i_haslam,k0)
  for(ii=0;ii<pars->n_nu;ii++) {
    double eval=gsl_vector_get(eigenval,ii);
    if(eval>=emax) {
      emax=eval;
      k0=ii;
    }
    gsl_vector_set(ev_haslam,ii,gsl_matrix_get(eigenvec,i_haslam,ii));
  }
  if(k0==-1) {
    fprintf(stderr,"CRIME: can't be\n");
    exit(1);
  }
  ev_h_k0=gsl_matrix_get(eigenvec,i_haslam,k0);

#pragma omp parallel default(none)		\
  shared(pars,ev_haslam,ev_h_k0,alms,k0)
  {
    int l;
    gsl_vector *aux_re=gsl_vector_alloc(pars->n_nu);
    gsl_vector *aux_im=gsl_vector_alloc(pars->n_nu);

#pragma omp for
    for(l=0;l<LMAX_HASLAM;l++) {
      int m;
      for(m=0;m<=l;m++) {
	int k;
	double re,im;
	long index=he_indexlm(l,m,pars->lmax);
	fcomplex a0=alms[k0][index];
	for(k=0;k<pars->n_nu;k++) {
	  gsl_vector_set(aux_re,k,creal(alms[k][index]));
	  gsl_vector_set(aux_im,k,cimag(alms[k][index]));
	}
	gsl_blas_ddot(ev_haslam,aux_re,&re);
	gsl_blas_ddot(ev_haslam,aux_im,&im);
	alms[k0][index]=a0-(re+I*im)/ev_h_k0;
      }
    }// end omp for

    gsl_vector_free(aux_re);
    gsl_vector_free(aux_im);
  }// end omp parallel
  
  gsl_vector_free(ev_haslam);
}

#ifdef _DEBUG
static void test_galaxy_maps(void)
{
  write_maps("haslam_read",NSIDE_TEMPLATE,1,&galaxy_template,NULL,NULL);
  write_maps("specin_read",NSIDE_TEMPLATE,1,&specin_template,NULL,NULL);
  write_maps("amplitude_read.fits",NSIDE_TEMPLATE,1,&galaxy_amplitude,NULL,NULL);
}
#endif //_DEBUG

void init_galaxy_template(ParamsForGet *pars)
{
  flouble *map_aux,*map_mean,*map_sigma;
  int *n_pixels;
  long nside_aux;
  long ii;

  printf("  * Setting galaxy templates\n");

  i_haslam=-1;
  for(ii=0;ii<pars->n_nu;ii++) {
    double nu0=pars->nutable[1][ii];
    double nuf=pars->nutable[2][ii];
    if((NU_HASLAM<=nuf)&&(NU_HASLAM>nu0)) i_haslam=ii;
  }
  if(i_haslam==-1) {
    fprintf(stderr,"CRIME: the Haslam map is outside the frequency range\n");
    exit(1);
  }

  galaxy_amplitude=(flouble *)my_malloc(nside2npix(NSIDE_TEMPLATE)*sizeof(flouble));
  galaxy_template=(flouble *)my_malloc(nside2npix(NSIDE_TEMPLATE)*sizeof(flouble));
  specin_template=(flouble *)my_malloc(nside2npix(NSIDE_TEMPLATE)*sizeof(flouble));
  
  printf("    Reading Haslam map\n");
  map_aux=he_read_healpix_map(pars->fname_haslam,&nside_aux,0);
  if(nside_aux!=NSIDE_TEMPLATE) {
    fprintf(stderr,"CRIME: can't be!\n");
    exit(1);
  }
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    galaxy_template[ii]=1000*map_aux[ii];
  }
  free(map_aux);

  printf("    Reading spectral index map\n");
  map_aux=he_read_healpix_map(pars->fname_specin,&nside_aux,0);
  if(nside_aux!=NSIDE_TEMPLATE) {
    fprintf(stderr,"CRIME: can't be!\n");
    exit(1);
  }
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    specin_template[ii]=map_aux[ii]-2.0+SPECIN_OFFSET;
  }
  free(map_aux);

  printf("    Calculating rms in regions\n");
  map_mean=(flouble *)my_calloc(nside2npix(NSIDE_REGIONS),sizeof(flouble));
  map_sigma=(flouble *)my_calloc(nside2npix(NSIDE_REGIONS),sizeof(flouble));
  n_pixels=(int *)my_calloc(nside2npix(NSIDE_REGIONS),sizeof(int));
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    double theta,phi,tval;
    long ipix_hi,ipix_lo;
    ipix_hi=ii;
    pix2ang_ring(NSIDE_TEMPLATE,ipix_hi,&theta,&phi);
    ang2pix_ring(NSIDE_REGIONS,theta,phi,&ipix_lo);
    tval=galaxy_template[ipix_hi];
    n_pixels[ipix_lo]++;
    map_mean[ipix_lo]+=tval;
    map_sigma[ipix_lo]+=tval*tval;
  }
  double max_sigma=0;
  for(ii=0;ii<nside2npix(NSIDE_REGIONS);ii++) {
    map_mean[ii]/=n_pixels[ii];
    map_sigma[ii]=map_sigma[ii]/n_pixels[ii]-
      map_mean[ii]*map_mean[ii];
    if(map_sigma[ii]>max_sigma) max_sigma=map_sigma[ii];
  }
  for(ii=0;ii<nside2npix(NSIDE_REGIONS);ii++) {
    map_sigma[ii]=sqrt(map_sigma[ii]/max_sigma);
  }
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    double theta,phi;
    long ipix_hi,ipix_lo;
    ipix_hi=ii;
    pix2ang_ring(NSIDE_TEMPLATE,ipix_hi,&theta,&phi);
    ang2pix_ring(NSIDE_REGIONS,theta,phi,&ipix_lo);
    galaxy_amplitude[ipix_hi]=map_sigma[ipix_lo];
  }
  free(map_mean);
  free(map_sigma);
  free(n_pixels);

  printf("    Removing redundant multipoles\n");
  fcomplex *alm_haslam=(fcomplex *)my_malloc(he_nalms(LMAX_TR_HASLAM)*
					     sizeof(fcomplex));
  he_map2alm(NSIDE_TEMPLATE,LMAX_TR_HASLAM,1,&galaxy_template,&alm_haslam);
  for(ii=LMAX_HASLAM;ii<=LMAX_TR_HASLAM;ii++) {
    int mm;
    for(mm=0;mm<=ii;mm++)
      alm_haslam[he_indexlm(ii,mm,LMAX_TR_HASLAM)]=0;
  }
  he_alm2map(NSIDE_TEMPLATE,LMAX_TR_HASLAM,1,&galaxy_template,&alm_haslam);

#ifdef _DEBUG
  test_galaxy_maps();
#endif //_DEBUG
}

void end_galaxy_template(void)
{
  free(galaxy_template);
  free(specin_template);
  free(galaxy_amplitude);
}
