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
#include "common_gh.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "cosmo_mad.h"

static double z_of_r_provisional(ParamGetHI *par,double r)
{
  if(r<=0) return 0;
  else if(r>=par->r_arr_z2r[NZ-1]) return par->z_arr_z2r[NZ-1];
  else {
    int iz=0;
    while(r>=par->r_arr_z2r[iz])
      iz++;
    return par->z_arr_z2r[iz-1]+(par->z_arr_z2r[iz]-par->z_arr_z2r[iz-1])*
      (r-par->r_arr_z2r[iz-1])/(par->r_arr_z2r[iz]-par->r_arr_z2r[iz-1]);
  }
}

double r_of_z(ParamGetHI *par,double z)
{
  if(z<=0) return 0;
  else if(z>=par->z_arr_z2r[NZ-1]) return par->r_arr_z2r[NZ-1];
  else {
    int iz=(int)(z/DZ);
    double r=par->r_arr_z2r[iz]+(par->r_arr_z2r[iz+1]-par->r_arr_z2r[iz])*
      (z-par->z_arr_z2r[iz])/DZ;
    return r;
  }
}

double z_of_r(ParamGetHI *par,double r)
{
  if(r<=0) return 0;
  else if(r>=par->r_arr_r2z[NZ-1]) return par->z_arr_r2z[NZ-1];
  else {
    int ir=(int)(r*par->glob_idr);
    double z=par->z_arr_r2z[ir]+(par->z_arr_r2z[ir+1]-par->z_arr_r2z[ir])*
      (r-par->r_arr_r2z[ir])*par->glob_idr;
    return z;
  }
}

double dgrowth_of_r(ParamGetHI *par,double r)
{
  if(r<=0) return 1;
  else if(r>=par->r_arr_r2z[NZ-1]) return par->growth_d_arr[NZ-1];
  else {
    int ir=(int)(r*par->glob_idr);
    double gd=par->growth_d_arr[ir]+(par->growth_d_arr[ir+1]-par->growth_d_arr[ir])*
      (r-par->r_arr_r2z[ir])*par->glob_idr;
    return gd;
  }
}

double vgrowth_of_r(ParamGetHI *par,double r)
{
  if(r<=0) return 1;
  else if(r>=par->r_arr_r2z[NZ-1]) return par->growth_v_arr[NZ-1];
  else {
    int ir=(int)(r*par->glob_idr);
    double gv=par->growth_v_arr[ir]+(par->growth_v_arr[ir+1]-par->growth_v_arr[ir])*
      (r-par->r_arr_r2z[ir])*par->glob_idr;
    return gv;
  }
}

static void int_error_handle(int status,double result,
                             double error)
{
  //////
  // Error handler for gsl
  if(isnan(result)) {
    fprintf(stderr,"CRIME: NAN found \n");
    exit(1);
  }
  else{
    if(status==GSL_EROUND)
      fprintf(stderr,"CRIME: Roundoff error: %lE %lE \n",result,error);
    else if(status==GSL_EMAXITER)
      fprintf(stderr,"CRIME: Ran out of iterations: %lE %lE \n",result,error);
    else if(status==GSL_ESING)
      fprintf(stderr,"CRIME: Singularity found: %lE %lE \n",result,error);
    else if(status==GSL_EDIVERGE) {
      fprintf(stderr,"CRIME: Integral seems to diverge: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ETOL) {
      fprintf(stderr,"CRIME: Can't reach tolerance: %lE %lE : %lE %%\n",
              result,error,100*error/result);
    }
    else if(status==GSL_EUNDRFLW)
      fprintf(stderr,"CRIME: Underflow: %lE %lE \n",result,error);
    else if(status==GSL_EDOM) {
      fprintf(stderr,"CRIME: Outside interpolation range!! %lE %lE\n",
	      result,error);
      exit(1);
    }
    else if(status) {
      fprintf(stderr,"CRIME: Unknown error code %d %lf %lf \n",
	      status,result,error);
      exit(1);
    }
  }
}

static double wind(double x,int setwf)
{
  //////
  // Window function:
  //  setwf=0 -> top hat
  //  setwf=1 -> gaussian
  //  setwf=2 -> sharp-k
  if(setwf==1) { //Gaussian
    return exp(-0.5*x*x);
  }
  else if(setwf==0) { //TopHat
    if(x<0.1) {
      return 1.-0.1*x*x+0.003571429*x*x*x*x-6.61376E-5*x*x*x*x*x*x
        +7.51563E-7*x*x*x*x*x*x*x*x;
    }
    else
      return 3*(sin(x)-x*cos(x))/(x*x*x);
  }
  else if(setwf==2) { //Sharp-k
    if(x<1) return 1;
    else return 0;
  }
  else
    return -1;
}

double pk_linear0(ParamGetHI *par,double lgk)
{
  //////
  // Linear power spectrum at redshift 0.
  // Extrapolated to ~k^ns for small k and
  // to k^{-3} for large k
  double pk;
  int ik=(int)((lgk-par->logkmin)*par->idlogk);
  
  if(ik<0)
    pk=par->pkarr[0]*pow(10,par->n_scal*(lgk-par->logkmin));
  else if(ik<par->numk)
    pk=par->pkarr[ik]+(lgk-par->logkarr[ik])*(par->pkarr[ik+1]-par->pkarr[ik])*par->idlogk;
  else
    pk=par->pkarr[par->numk-1]*pow(10,-3*(lgk-par->logkmax));

  return pk;
}

static double j_bessel_0(double x)
{
  //////
  // Bessel's j0 function
  if(x>0.1)
    return sin(x)/x;
  else
    return 1.-0.166667*x*x+0.00833333*x*x*x*x-
      0.000198413*x*x*x*x*x*x+2.75573E-6*x*x*x*x*x*x*x*x;
}

typedef struct { //Param struct for integrals
  double r;
  double R1;
  double R2;
  int wf1;
  int wf2;
  ParamGetHI *par;
} xiparam;

static double integxiL_O(double kk,void *params)
{
  //////
  // integrand for xi_L (for oscillatory integration)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  double lgk=log10(kk);
  par=(xiparam *)params;

  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);

  dum=TWOPIPIINV*pk_linear0(par->par,lgk)*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)/xr;

  return dum;
}

static double integxiL_NO(double logk,void *params)
{
  //////
  // integrand for xi_L (including oscillatory term in j0)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  par=(xiparam *)params;

  double kk=pow(10,logk);
  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);

  dum=TWOPIPIINVLOGTEN*pk_linear0(par->par,logk)*kk*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)*j_bessel_0(xr);

  return dum;
}

static double xi2p_L(ParamGetHI *par,double r,double R1,double R2,
		     char *wf1,char *wf2,double errfac)
{
  //////
  // Correlation function between the linear density contrast smoothed
  // with window function (wf1,R1) and with window function (wf2,R2)
  // at two points separated by a distance r:
  //              <delta_(R1,wf1)(x)*delta_(R2,wf2)(x+r)>
  gsl_function integrand;
  double relerrt=1E-4;
  double integral,errintegral;
  xiparam xpar;
  double lim=MIN(R1,R2);
  lim/=r;

  xpar.r=r;
  xpar.R1=R1;
  xpar.R2=R2;
  xpar.par=par;
  if(!strcmp(wf1,"Gauss"))
    xpar.wf1=1;
  else if(!strcmp(wf1,"TopHat"))
    xpar.wf1=0;
  else {
    fprintf(stderr,"CRIME: Unknown window function %s \n",wf1);
    exit(1);
  }
  if(!strcmp(wf2,"Gauss"))
    xpar.wf2=1;
  else if(!strcmp(wf2,"TopHat"))
    xpar.wf2=0;
  else {
    fprintf(stderr,"CRIME: Unknown window function %s \n",wf2);
    exit(1);
  }

  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(1000);
  integrand.params=&xpar;
  if(lim>=1) {
    integrand.function=&integxiL_NO;
    int stat=gsl_integration_qagil(&integrand,par->logkmax,0,relerrt,1000,w,
                                   &integral,&errintegral);
    int_error_handle(stat,integral,errintegral);
  }
  else {
    lim*=errfac;
    gsl_integration_workspace *cw
      =gsl_integration_workspace_alloc(1000);
    gsl_integration_qawo_table *wf
      =gsl_integration_qawo_table_alloc(r,0.1,GSL_INTEG_SINE,100);

    integrand.function=&integxiL_O;
    int stat=gsl_integration_qawf(&integrand,0,relerrt*lim,1000,
                           w,cw,wf,&integral,&errintegral);
    int_error_handle(stat,integral,errintegral);

    gsl_integration_qawo_table_free(wf);
    gsl_integration_workspace_free(cw);
  }
  gsl_integration_workspace_free(w);

  return integral;
}

static double sigL2(ParamGetHI *par,double R1,double R2,char *wf1,char *wf2)
{
  //////
  // Covariance between the linear density contrast smoothed with
  // window function (wf1,R1) and with window function (wf2,R2) at
  // the same point:  <delta_(R1,wf1)(x)*delta_(R2,wf2)(x)>
  return xi2p_L(par,0,R1,R2,wf1,wf2,1);
}

static void pk_linear_set(ParamGetHI *par)
{
  //////
  // Reads linear power spectrum. CAMB format expected.
  int ii;
  double kk,ppk;
  FILE *fpk;

  print_info("Reading P_k from file: %s\n",par->fnamePk);
  fpk=fopen(par->fnamePk,"r");
  if(fpk==NULL) error_open_file(par->fnamePk);

  par->numk=linecount(fpk);
  par->logkarr=(double *)my_malloc(par->numk*sizeof(double));
  par->pkarr=(double *)my_malloc(par->numk*sizeof(double));
  rewind(fpk);
  for(ii=0;ii<par->numk;ii++) {
    int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
    if(stat!=2) error_read_line(par->fnamePk,ii+1);
    par->pkarr[ii]=ppk;
    par->logkarr[ii]=log10(kk); //log(k) in h Mpc^-1
  }
  par->logkmin=par->logkarr[0];
  par->logkmax=par->logkarr[par->numk-1];
  par->idlogk=(par->numk-1)/(par->logkarr[par->numk-1]-par->logkarr[0]);
  fclose(fpk);

  // normalize
  double norm_pk=par->sig8*par->sig8/sigL2(par,8,8,"TopHat","TopHat");
  print_info("  Original sigma8=%lf\n",
	     sqrt(sigL2(par,8,8,"TopHat","TopHat")));
  for(ii=0;ii<par->numk;ii++)
    par->pkarr[ii]*=norm_pk;
}

void cosmo_set(ParamGetHI *par)
{
  //////
  // This initializes the cosmological model
  // at redshift z_s

  int ii;

  Csm_params *pars=csm_params_new();
  csm_unset_gsl_eh();
  csm_background_set(pars,par->OmegaM,par->OmegaL,par->OmegaB,par->weos,0,par->hhub,2.275);

  par->fgrowth_0=csm_f_growth(pars,1);
  par->hubble_0=csm_hubble(pars,1);

  par->z_min=NU_21/par->nu_max-1;
  par->z_max=NU_21/par->nu_min-1;
  par->r_min=csm_radial_comoving_distance(pars,1/(1+par->z_min));
  par->r_max=csm_radial_comoving_distance(pars,1/(1+par->z_max));

  par->l_box=2*par->r_max*(1+2./par->n_grid);
  par->pos_obs[0]=0.5*par->l_box;
  par->pos_obs[1]=0.5*par->l_box;
  par->pos_obs[2]=0.5*par->l_box;

  //Set z-dependent functions
  for(ii=0;ii<NZ;ii++) {
    double z=((double)ii)*DZ;
    double a=1/(1+z);
    double rz=csm_radial_comoving_distance(pars,a);
    par->z_arr_z2r[ii]=z;
    par->r_arr_z2r[ii]=rz;
  }
  if((par->z_arr_z2r[NZ-1]<=par->z_max)||(par->r_arr_z2r[NZ-1]<=par->r_max))
    report_error(1,"OMG!\n");
  
  double growth0=csm_growth_factor(pars,1);
  par->glob_idr=(NZ-1)/(par->r_arr_z2r[NZ-1]-par->r_arr_z2r[0]);
  for(ii=0;ii<NZ;ii++) {
    double r=ii/par->glob_idr;
    double z=z_of_r_provisional(par,r);
    double a=1/(1+z);
    par->z_arr_r2z[ii]=z;
    par->r_arr_r2z[ii]=r;

    double gz=csm_growth_factor(pars,a)/growth0;
    double fz=csm_f_growth(pars,a);
    double hhz=csm_hubble(pars,a);
    par->growth_d_arr[ii]=gz;
    //This is for the comoving velocity
    par->growth_v_arr[ii]=(gz*hhz*fz)/(par->fgrowth_0*par->hubble_0);
  }
  if((par->z_arr_r2z[NZ-1]<=par->z_max)||(par->r_arr_r2z[NZ-1]<=par->r_max))
    report_error(1,"OMG!\n");

#ifdef _DEBUG
  int nz=500;
  FILE *fil=fopen("test_cosmo.dat","w");
  for(ii=0;ii<nz;ii++) {
    double z1=5.*(ii+0.5)/nz;
    double r=csm_radial_comoving_distance(pars,1/(1+z1));
    double z2=z_of_r(par,r);
    double gz1=dgrowth_of_r(par,r);
    double gz2=csm_growth_factor(pars,1/(1+z1))/growth0;

    fprintf(fil,"%lE %lE %lE %d %lE %lE\n",z1,r,z2,nz,gz1,gz2);
  }
  fclose(fil);
#endif //_DEBUG
  csm_params_free(pars);

  pk_linear_set(par);
}
