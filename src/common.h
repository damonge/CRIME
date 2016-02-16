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
#ifndef _COMMON_
#define _COMMON_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#ifdef _HAVE_OMP
#include <omp.h>
#endif //_HAVE_OMP
#include <gsl/gsl_rng.h>
#include "fftw3.h"

#define NU_21 1420.40575177
#define DYNAMIC_SIZE 1
#define RTOD 57.2957795
#define DTOR 0.01745329251
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
typedef fftwf_complex dftw_complex;
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
typedef fftw_complex dftw_complex;
#endif //_SPREC

void *my_malloc(size_t size);

void *my_calloc(size_t nmemb,size_t size);

void error_open_file(char *fname);

void error_read_line(char *fname,int nlin);

int linecount(FILE *f);

void timer(int i);

gsl_rng *init_rng(unsigned int seed);

double rng_01(gsl_rng *rng);

int rng_poisson(double lambda,gsl_rng *rng);

void rng_delta_gauss(double *module,double *phase,
		     gsl_rng *rng,double sigma2);

void end_rng(gsl_rng *rng);


//////
// Functions defined in healpix_extra.c
long he_nalms(int lmax);

long he_indexlm(int l,int m,int lmax);

void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);

void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);

void he_write_healpix_map(float **tmap,int nfields,long nside,char *fname);

flouble *he_read_healpix_map(char *fname,long *nside,int nfield);

int he_ring_num(long nside,double z);

long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip);

void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest);

double *he_generate_beam_window(int lmax,double fwhm_amin);

void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window);

#endif //_COMMON_
