#ifndef _COMMON_CRIME_
#define _COMMON_CRIME_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdarg.h>
#ifdef _HAVE_OMP
#include <omp.h>
#endif //_HAVE_OMP
#include <fftw3.h>
#ifdef _HAVE_MPI
#include <mpi.h>
#include <fftw3-mpi.h>
#endif //_HAVE_MPI
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <fitsio.h>
#include <chealpix.h>

#ifdef _VERBOSE
#define VERBOSITY 1
#else //_VERBOSE
#define VERBOSITY 0
#endif //_VERBOSE

#define NU_21 1420.40575177
#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
typedef fftwf_complex dftw_complex;
#ifdef _HAVE_MPI
#define FLOUBLE_MPI MPI_FLOAT
#endif //_HAVE_MPI
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
typedef fftw_complex dftw_complex;
#ifdef _HAVE_MPI
#define FLOUBLE_MPI MPI_DOUBLE
#endif //_HAVE_MPI
#endif //_SPREC

extern int NodeThis;
extern int NodeLeft;
extern int NodeRight;
extern int NNodes;
extern int IThread0;
extern int MPIThreadsOK;

//Defined in common.c
typedef struct {
#ifdef _MY_SPLINE
  int n;
  double *x;
  double *y;
  double idx;
#else //_MY_SPLINE
  gsl_interp_accel *intacc;
  gsl_spline *spline;
#endif //_MY_SPLINE
  double x0,xf;
  double y0,yf;
} SplPar;
typedef struct {
  int n_nu;
  double *nu_0;
  double *nu_f;
  long n_side;
} ParamMaps;
void mpi_init(int* p_argc,char*** p_argv);
void report_error(int level,char *fmt,...);
void print_info(int level,char *fmt,...);
void *my_malloc(size_t size);
size_t my_fread(void *ptr,size_t size,size_t nmemb,FILE *stream);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
int linecount(FILE *f);
SplPar *spline_init(int n,double *x,double *y,double y0,double yf);
double spline_eval(double x,SplPar *spl);
void spline_free(SplPar *spl);
void read_nulist(char *fname,ParamMaps *pmap);
int get_inu(ParamMaps *pmap,double nu,int inu_start);
void param_maps_free(ParamMaps *pmap);

//Defined in rng.c
#define RNG_NRAN 624
#define RNG_MRAN 397
#define RNG_MATRIX_A 0x9908b0df
#define RNG_UPPER_MASK 0x80000000UL
#define RNG_LOWER_MASK 0x7fffffffUL
typedef struct {
  unsigned long mt[RNG_NRAN];
  int mti;
  int calc_gauss;
  double u;
  double phi;
} Rng;
Rng *init_rng(unsigned long seed);
void end_rng(Rng *rng);
unsigned long rand_ulong(Rng *rng);
double rand_real01(Rng *rng);
double rand_gauss(Rng *rng);

//Defined in healpix_extra.c
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,
		int pol_1,int pol_2,
		flouble **cls,int nside,int lmax);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest);
double *he_generate_beam_window(int lmax,double fwhm_amin);
void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window);


#endif //_COMMON_CRIME_
