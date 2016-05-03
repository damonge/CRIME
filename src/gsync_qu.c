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
#include <chealpix.h>

#define NL_POL 383 //3*NSIDE-1
#define RES_X_FAC 0.1
#define RANGE_X_FAC 5
#define CLIGHT 299.792458 // c/(1 MHz * 1 m)
#define NU_HI 23000.0 //23 GHz
#define P_HILAT 0.24 //High latitude polarized fraction
#define HILAT_ANGLE 20.0 // Angle defining the NGC

static int glb_nx;
static double glb_range_x,glb_x_min,glb_sigmax;
static double xi_pol; //Correlation length in Faraday space
static double t2mean_hilat_hifreq; //Mean temperature of the high-frequency map
static flouble *map_sigma;
static flouble *map_specin;
static flouble cls_pol[NL_POL+1];

static void end_polarization(void)
{
  free(map_sigma);
  free(map_specin);
}

static void init_polarization(ParamsForGet *pars)
{
  long ii,nside_aux;
  flouble *map_T_dum,*map_specin_dum;
  int *npix_aux;

  xi_pol=pars->xi_polarization;

  printf("  * Setting templates\n");
  printf("    Reading Faraday widths from %s\n",pars->fname_faraday);
  map_T_dum=he_read_healpix_map(pars->fname_faraday,&nside_aux,0);
  if(nside_aux!=NSIDE_POL) {
    fprintf(stderr,"CRIME: can't be!\n");
    exit(1);
  }
  map_sigma=(flouble *)my_malloc(nside2npix(NSIDE_POL)*sizeof(flouble));
  for(ii=0;ii<nside2npix(NSIDE_POL);ii++) {
    map_sigma[ii]=map_T_dum[ii];
  }
  free(map_T_dum);

  double sigmax=10000,sigmin=-10000,xmeanmax=-10000,xmeanmin=10000;
  for(ii=0;ii<nside2npix(NSIDE_POL);ii++) {
    if(map_sigma[ii]>=sigmin)
      sigmin=map_sigma[ii];
    if(map_sigma[ii]<=sigmax)
      sigmax=map_sigma[ii];
  }
  sigmin=1./sigmin; sigmax=1./sigmax; glb_sigmax=sigmax;

  for(ii=0;ii<pars->n_nu;ii++) {
    double nu=pars->nutable[0][ii];
    double x_nu=2*(CLIGHT*CLIGHT)/(nu*nu);
    if(x_nu>=xmeanmax)
      xmeanmax=x_nu;
    if(x_nu<=xmeanmin)
      xmeanmin=x_nu;
  }

  double xmin,xmax,dx;
  xmin=xmeanmin-RANGE_X_FAC*sigmax;
  xmax=xmeanmax+RANGE_X_FAC*sigmax;
  dx=RES_X_FAC*sigmin;
  glb_x_min=xmin;
  glb_range_x=xmax-xmin;
  glb_nx=(int)(glb_range_x/dx+0.5);
#ifdef _DEBUG
  printf("    X range: (%lE , %lE), %lE %lE %lE %d\n",xmin,xmax,xmeanmin,xmeanmax,dx,glb_nx);
#endif //_DEBUG

  printf("    Reading spectral indices from %s\n",pars->fname_specin);
  map_specin_dum=he_read_healpix_map(pars->fname_specin,&nside_aux,0);
  if(nside_aux!=NSIDE_TEMPLATE) {
    fprintf(stderr,"CRIME: can't be!\n");
    exit(1);
  }
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    map_specin_dum[ii]+=SPECIN_OFFSET-2;
  }

  printf("    Reading Haslam map from %s\n",pars->fname_haslam);
  map_T_dum=he_read_healpix_map(pars->fname_haslam,&nside_aux,0);
  if(nside_aux!=NSIDE_TEMPLATE) {
    fprintf(stderr,"CRIME: can't be!\n");
    exit(1);
  }
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    map_T_dum[ii]*=1000;
  }

  printf("    Computing high-latitude mean temperature at 23 GHz: ");
  long npix_strip;
  long *listpix=he_query_strip(NSIDE_TEMPLATE,0.0,
			       DTOR*HILAT_ANGLE,&npix_strip);
  t2mean_hilat_hifreq=0;
  for(ii=0;ii<npix_strip;ii++) {
    long ipix=listpix[ii];
    double tmp=map_T_dum[ipix];
    tmp*=pow(NU_HI/NU_HASLAM,map_specin_dum[ipix]);
    t2mean_hilat_hifreq+=tmp*tmp;
  }
  t2mean_hilat_hifreq/=npix_strip;
  printf("%lE mK\n",sqrt(t2mean_hilat_hifreq));
  free(listpix);
  free(map_T_dum);

  npix_aux=(int *)my_calloc(nside2npix(NSIDE_POL),sizeof(int));
  map_specin=(flouble *)my_calloc(nside2npix(NSIDE_POL),sizeof(flouble));
  for(ii=0;ii<nside2npix(NSIDE_TEMPLATE);ii++) {
    long iph,ipl;
    double theta,phi;

    iph=ii;
    pix2ang_ring(NSIDE_TEMPLATE,iph,&theta,&phi);
    ang2pix_ring(NSIDE_POL,theta,phi,&ipl);
    npix_aux[ipl]++;
    map_specin[ipl]+=map_specin_dum[iph];
  }
  for(ii=0;ii<nside2npix(NSIDE_POL);ii++)
    map_specin[ii]/=npix_aux[ii];
  free(npix_aux);
  free(map_specin_dum);

  for(ii=0;ii<=NL_POL;ii++) {
    if(ii<pars->lmin) {
      cls_pol[ii]=pow((double)(pars->lref)/(pars->lmin),
		      pars->beta_polarization);
    }
    else {
      cls_pol[ii]=pow((double)(pars->lref)/ii,
		      pars->beta_polarization);
    }
  }
}

static void get_uncorrelated_maps(ParamsForGet *pars,flouble ***maps_out_q,
				  flouble ***maps_out_u)
{
  int ii;

  fcomplex **alms_q=(fcomplex **)my_malloc(glb_nx*sizeof(fcomplex *));
  fcomplex **alms_u=(fcomplex **)my_malloc(glb_nx*sizeof(fcomplex *));
  flouble **maps_q=(flouble **)my_malloc(glb_nx*sizeof(flouble *));
  flouble **maps_u=(flouble **)my_malloc(glb_nx*sizeof(flouble *));

  for(ii=0;ii<glb_nx;ii++) {
    alms_q[ii]=(fcomplex *)my_malloc(he_nalms(NL_POL)*sizeof(fcomplex));
    alms_u[ii]=(fcomplex *)my_malloc(he_nalms(NL_POL)*sizeof(fcomplex));
    maps_q[ii]=(flouble *)my_malloc(nside2npix(NSIDE_POL)*sizeof(flouble));
    maps_u[ii]=(flouble *)my_malloc(nside2npix(NSIDE_POL)*sizeof(flouble));
  }

#ifdef _HAVE_OMP
#pragma omp parallel default(none)				\
  shared(pars,alms_q,alms_u,cls_pol,glb_nx,maps_q,maps_u)
#endif //_HAVE_OMP
  {
    int ix;
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
    for(ix=0;ix<glb_nx;ix++) {
      int l;
      for(l=0;l<=NL_POL;l++) {
	int m;
	double mod,phase;
	double sigma2=cls_pol[l];
	long index=he_indexlm(l,0,NL_POL);
	
	rng_delta_gauss(&mod,&phase,rng,sigma2);
	alms_q[ix][index]=mod*cos(phase);
	alms_u[ix][index]=mod*sin(phase);
	
	for(m=1;m<=l;m++) {
	  fcomplex zp,zm;
	  index=he_indexlm(l,m,NL_POL);
	  
	  rng_delta_gauss(&mod,&phase,rng,sigma2);
	  zp=mod*(cos(phase)+I*sin(phase));
	  rng_delta_gauss(&mod,&phase,rng,sigma2);
	  zm=mod*(cos(phase)+I*sin(phase));
	  
	  alms_q[ix][index]=0.5*(zp+zm);
	  alms_u[ix][index]=-0.5*I*(zp-zm);
	}
      }
#ifndef _OLD_ORDER
      he_alm2map(NSIDE_POL,NL_POL,1,&(maps_q[ix]),&(alms_q[ix]));
      he_alm2map(NSIDE_POL,NL_POL,1,&(maps_u[ix]),&(alms_u[ix]));
#endif //_OLD_ORDER
    }// end omp for
  }// end omp parallel

#ifdef _OLD_ORDER
  printf("Transforming to real space\n");
  he_alm2map(NSIDE_POL,NL_POL,glb_nx,maps_q,alms_q);
  he_alm2map(NSIDE_POL,NL_POL,glb_nx,maps_u,alms_u);
#endif //_OLD_ORDER

  for(ii=0;ii<glb_nx;ii++) {
    free(alms_q[ii]);
    free(alms_u[ii]);
  }
  free(alms_q);
  free(alms_u);
 
  *maps_out_q=maps_q;
  *maps_out_u=maps_u;
}

static void integrate_LOS_maps(double nu,flouble **maps_in_q,flouble **maps_in_u,
			       flouble *map_out_q,flouble *map_out_u)
{
  int ix;
  long ipix;
  double x_nu=2*(CLIGHT*CLIGHT)/(nu*nu);
  double dx=glb_range_x/glb_nx;
#ifdef _OLD_ORDER
  double x_min=x_nu-RANGE_X_FAC*glb_sigmax;
  double x_max=x_nu+RANGE_X_FAC*glb_sigmax;

  for(ix=0;ix<glb_nx;ix++) {
    double x=glb_x_min+(ix+0.5)*dx;
    if((x<=x_max) && (x>=x_min)) {
      for(ipix=0;ipix<nside2npix(NSIDE_POL);ipix++) {
	double sigma=map_sigma[ipix];
	double prefac=sqrt(sigma)*dx;
	double kernel=prefac*exp(-0.5*(x-x_nu)*(x-x_nu)*sigma*sigma-
				 0.25*x*x*xi_pol*xi_pol);
	map_out_q[ipix]+=kernel*maps_in_q[ix][ipix];
	map_out_u[ipix]+=kernel*maps_in_u[ix][ipix];
      }
    }
  }
#else //_OLD_ORDER

  for(ipix=0;ipix<nside2npix(NSIDE_POL);ipix++) {
    double sigma=map_sigma[ipix];
    double x_min=x_nu-RANGE_X_FAC/sigma;
    double x_max=x_nu+RANGE_X_FAC/sigma;
    double prefac=sqrt(sigma)*dx;
    for(ix=0;ix<glb_nx;ix++) {
      double x=glb_x_min+(ix+0.5)*dx;
      if((x<=x_max) && (x>=x_min)) {
	double kernel=prefac*exp(-0.5*(x-x_nu)*(x-x_nu)*sigma*sigma-
				 0.25*x*x*xi_pol*xi_pol);
	map_out_q[ipix]+=kernel*maps_in_q[ix][ipix];
	map_out_u[ipix]+=kernel*maps_in_u[ix][ipix];
      }
    }
  }
#endif //_OLD_ORDER
}

void do_polarization(ParamsForGet *pars,flouble ***maps_nuspace_q,
		     flouble ***maps_nuspace_u)
{
  long ii;
  flouble **maps_x_q,**maps_x_u;
  flouble **maps_nu_q,**maps_nu_u;
  flouble *map_aux_q,*map_aux_u;

  printf("*** Doing polarized synchrotron\n");
  init_polarization(pars);

  printf("  Generating uncorrelated alms\n");
  timer(0);
  get_uncorrelated_maps(pars,&maps_x_q,&maps_x_u);
  timer(1);
  
  printf("  Normalizing to polarized fraction at high latitudes\n");
  timer(0);
  map_aux_q=(flouble *)my_calloc(nside2npix(NSIDE_POL),sizeof(flouble));
  map_aux_u=(flouble *)my_calloc(nside2npix(NSIDE_POL),sizeof(flouble));
  integrate_LOS_maps(NU_HI,maps_x_q,maps_x_u,map_aux_q,map_aux_u);

  long npix_strip;
  long *listpix=he_query_strip(NSIDE_POL,0.0,DTOR*HILAT_ANGLE,&npix_strip);
  double normalization=0;
  for(ii=0;ii<npix_strip;ii++) {
    long index=listpix[ii];
    double q=map_aux_q[index];
    double u=map_aux_u[index];
    normalization+=q*q+u*u;
  }
  normalization/=npix_strip;
  normalization=sqrt(P_HILAT*P_HILAT/(1-P_HILAT*P_HILAT)*(t2mean_hilat_hifreq/normalization));
  free(listpix);
  free(map_aux_q);
  free(map_aux_u);
  timer(1);

  printf("  Integrating LOSs and applying spectral indices\n");
  timer(0);
  maps_nu_q=(flouble **)my_malloc(pars->n_nu*sizeof(flouble *));
  maps_nu_u=(flouble **)my_malloc(pars->n_nu*sizeof(flouble *));
  for(ii=0;ii<pars->n_nu;ii++) {
    maps_nu_q[ii]=(flouble *)my_calloc(nside2npix(NSIDE_POL),sizeof(flouble));
    maps_nu_u[ii]=(flouble *)my_calloc(nside2npix(NSIDE_POL),sizeof(flouble));
  }

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(maps_nu_q,maps_nu_u,maps_x_q,maps_x_u) \
  shared(pars,map_specin,normalization)
#endif //_HAVE_OMP
  {
    long i;

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)    
#endif //_HAVE_OMP
    for(i=0;i<pars->n_nu;i++) {
      long ip;
      double nufrac=pars->nutable[0][i]/NU_HI;

      integrate_LOS_maps(pars->nutable[0][i],maps_x_q,maps_x_u,
      			 maps_nu_q[i],maps_nu_u[i]);
      for(ip=0;ip<nside2npix(NSIDE_POL);ip++) {
	double powerlaw=pow(nufrac,map_specin[ip]);
	maps_nu_q[i][ip]*=powerlaw*normalization;
	maps_nu_u[i][ip]*=powerlaw*normalization;
      }
    }// end omp for
  }// end omp parallel
  timer(1);
  
  for(ii=0;ii<glb_nx;ii++) {
    free(maps_x_q[ii]);
    free(maps_x_u[ii]);
  }
  free(maps_x_q);
  free(maps_x_u);

  *maps_nuspace_q=maps_nu_q;
  *maps_nuspace_u=maps_nu_u;

  end_polarization();
}
