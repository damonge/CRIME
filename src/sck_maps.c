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
#include <gsl/gsl_eigen.h>
#include <chealpix.h>

static double *cls_arr;

static void set_cls(ParamsForGet *pars)
{
  int l;
  
  cls_arr=(double *)my_malloc((pars->lmax+1)*sizeof(double));
  for(l=0;l<=pars->lmax;l++) {
    if(l<pars->lmin) {
      cls_arr[l]=(pars->amp)*
	pow(((double)(pars->lref))/(pars->lmin),pars->beta);
    }
    else {
      cls_arr[l]=(pars->amp)*
	pow(((double)(pars->lref))/l,pars->beta);
    }
  }
}

static void end_cls(void)
{
  free(cls_arr);
}

static fcomplex *get_alms(ParamsForGet *pars,double scale,gsl_rng *rng)
{
  fcomplex *alms;

  if(scale<=0) {
    alms=(fcomplex *)my_calloc(he_nalms(pars->lmax),sizeof(fcomplex));
  }
  else {
    int ll;

    alms=(fcomplex *)my_malloc(he_nalms(pars->lmax)*sizeof(fcomplex));
    for(ll=0;ll<=pars->lmax;ll++) {
      int mm;
      double sigma2=cls_arr[ll]*scale;
      for(mm=0;mm<=ll;mm++) {
	double mod,phase;
	rng_delta_gauss(&mod,&phase,rng,sigma2);
	alms[he_indexlm(ll,mm,pars->lmax)]=mod*(cos(phase)+I*sin(phase));
      }
    }
  }

  return alms;
}

static void eigen_freq(ParamsForGet *pars,gsl_vector *eigenval,gsl_matrix *eigenvec)
{
  int ii;
  gsl_matrix *corrmat=gsl_matrix_alloc(pars->n_nu,pars->n_nu);
  gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(pars->n_nu);

  for(ii=0;ii<pars->n_nu;ii++) {
    int jj;
    double nu1=pars->nutable[0][ii];
    for(jj=0;jj<pars->n_nu;jj++) {
      double nu2=pars->nutable[0][jj];
      double dum=log(nu1/nu2)/pars->xi;
      double m_element=pow((nu1*nu2)/(pars->nu_ref*pars->nu_ref),
			   -pars->alpha)*exp(-0.5*dum*dum);
      gsl_matrix_set(corrmat,ii,jj,m_element);
    }
  }
  gsl_eigen_symmv(corrmat,eigenval,eigenvec,w);

  for(ii=0;ii<pars->n_nu;ii++) {
    double ev=gsl_vector_get(eigenval,ii);
    if(ev<=0) gsl_vector_set(eigenval,ii,0);
  }

  gsl_eigen_symmv_free(w);
  gsl_matrix_free(corrmat);
}

flouble **get_sck_maps(ParamsForGet *pars)
{
  //CBI generate alms and transform only for non-zero eigenvalues
  int ii;
  flouble **maps_corr;
  fcomplex **alms_uncorr=(fcomplex **)my_malloc(pars->n_nu*sizeof(flouble *));
  fcomplex **alms_corr=(fcomplex **)my_malloc(pars->n_nu*sizeof(flouble *));
  gsl_vector *eigenval=gsl_vector_alloc(pars->n_nu);
  gsl_matrix *eigenvec=gsl_matrix_alloc(pars->n_nu,pars->n_nu);

  printf("*** Generating SCK maps\n");
  printf("  Diagonalizing covariance\n");
  timer(0);
  eigen_freq(pars,eigenval,eigenvec);
  timer(1);

  printf("  Generating uncorrelated alms\n");
  timer(0);
  set_cls(pars);
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(pars,alms_uncorr,eigenval)
#endif //_HAVE_OMP
  {
    int i;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=pars->seed+ithr;
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(i=0;i<pars->n_nu;i++) {
      double scale=gsl_vector_get(eigenval,i);
      alms_uncorr[i]=get_alms(pars,scale,rng_thr);
    }// end omp for
  }//end omp parallel
  end_cls();
  timer(1);

  if(pars->do_galaxy) {
    printf("  Constraining alms in Haslam scales\n");
    timer(0);
    init_galaxy_template(pars);
    constrain_alms(pars,eigenvec,eigenval,alms_uncorr);
    timer(1);
  }

  printf("  Correlating alms\n");
  timer(0);
  for(ii=0;ii<pars->n_nu;ii++)
    alms_corr[ii]=(fcomplex *)my_malloc(he_nalms(pars->lmax)*sizeof(fcomplex));

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(pars,alms_uncorr,eigenvec,alms_corr)
#endif //_HAVE_OMP
  {
    int i;
    gsl_vector *aux1_re=gsl_vector_alloc(pars->n_nu);
    gsl_vector *aux1_im=gsl_vector_alloc(pars->n_nu);
    gsl_vector *aux2_re=gsl_vector_alloc(pars->n_nu);
    gsl_vector *aux2_im=gsl_vector_alloc(pars->n_nu);

#ifdef _HAVE_OMP    
#pragma omp for
#endif //_HAVE_OMP
    for(i=0;i<he_nalms(pars->lmax);i++) {
      int j;
      for(j=0;j<pars->n_nu;j++) {
	gsl_vector_set(aux1_re,j,creal(alms_uncorr[j][i]));
	gsl_vector_set(aux1_im,j,cimag(alms_uncorr[j][i]));
      }
      gsl_blas_dgemv(CblasNoTrans,1,eigenvec,aux1_re,0,aux2_re);
      gsl_blas_dgemv(CblasNoTrans,1,eigenvec,aux1_im,0,aux2_im);
      for(j=0;j<pars->n_nu;j++) {
	double re=gsl_vector_get(aux2_re,j);
	double im=gsl_vector_get(aux2_im,j);
	alms_corr[j][i]=re+I*im;
      }
    }//end omp for

    gsl_vector_free(aux1_re);
    gsl_vector_free(aux1_im);
    gsl_vector_free(aux2_re);
    gsl_vector_free(aux2_im);
  }//end omp parallel
  timer(1);

  maps_corr=(flouble **)my_malloc(pars->n_nu*sizeof(flouble *));
  for(ii=0;ii<pars->n_nu;ii++) {
    free(alms_uncorr[ii]);
    maps_corr[ii]=(flouble *)my_malloc(nside2npix(pars->nside)*sizeof(flouble));
  }
  free(alms_uncorr);

  printf("  Harmonic-transforming\n");
  timer(0);
  he_alm2map(pars->nside,pars->lmax,pars->n_nu,maps_corr,alms_corr);
  timer(1);

  if(pars->do_galaxy) {
    printf("  Adding extrapolated Haslam map\n");
    timer(0);
    add_haslam(pars,maps_corr);
    timer(1);
    end_galaxy_template();
  }

  for(ii=0;ii<pars->n_nu;ii++)
    free(alms_corr[ii]);
  free(alms_corr);

  printf("\n");
  return maps_corr;
}
