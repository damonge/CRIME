#include "common.h"
#include "common_gh.h"

#define AINIT_GROWTH 1E-6
#define RELERRT_DIST 1E-6
#define HMPC 2997.92458 //H0^-1 in Mpc/h
#define DZ 0.001
#define NZ 5001

static void int_error_handle(int status,double result,
                             double error)
{
  if(isnan(result)) {
    report_error(1,"NAN found \n");
  }
  else{
    if(status==GSL_EROUND) {
      report_error(0,"Roundoff error: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EMAXITER) {
      report_error(0,"Ran out of iterations: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ESING) {
      report_error(0,"Singularity found: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EDIVERGE) {
      report_error(0,"Integral seems to diverge: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ETOL) {
      report_error(0,"Can't reach tolerance: %lE %lE\n",
	      result,error);
    }
    else if(status==GSL_EUNDRFLW)
      report_error(0,"Underflow: %lE %lE \n",result,error);
    else if(status==GSL_EDOM) {
      report_error(1,"Outside interpolation range!! %lE %lE\n",
	      result,error);
    }
    else if(status) {
      report_error(1,"Unknown error code %d %lf %lf \n",
	      status,result,error);
    }
  }
}

static double hubble(ParamGetHI *par,double a)
{
  return sqrt((par->OmegaM+par->OmegaL*pow(a,-3*(par->w0+par->wa))*
	       exp(3*(a-1)*par->wa)+par->OmegaK*a)/(a*a*a));
}

static double omega_matter(ParamGetHI *par,double a)
{
  return par->OmegaM/(par->OmegaM+par->OmegaL*pow(a,-3*(par->w0+par->wa))*
		      exp(3*(a-1)*par->wa)+par->OmegaK*a);
}

static double integ_chi(double z,void *params)
{
  ParamGetHI *par=(ParamGetHI *)params;
  return 1./hubble(par,1./(1+z));
}

static double chi(ParamGetHI *par,double z)
{
  if(z==0)
    return 0;
  else {
    double result,error;
    size_t sdum;
    gsl_function integrand;
    int stat;
    
    integrand.function=&integ_chi;
    integrand.params=par;
    stat=gsl_integration_qng(&integrand,0,z,0,RELERRT_DIST,
			     &result,&error,&sdum);
    int_error_handle(stat,result,error);
    
    return result*HMPC;
  }
}

static int growth_ode_system(double a,const double y[],double dydt[],void *params)
{
  ParamGetHI *par=(ParamGetHI *)params;
  double hnorm=hubble(par,a);
  
  dydt[0]=y[1]/(a*a*a*hnorm);
  dydt[1]=1.5*hnorm*a*omega_matter(par,a)*y[0];

  return GSL_SUCCESS;
}

static void growth_factor_and_growth_rate(ParamGetHI *par,double z,double *gf,double *fg)
{
  double a=1./(1+z);

  if(a<AINIT_GROWTH) {
    *gf=a;
    *fg=1;
  }
  else {
    double y[2];
    double ainit=AINIT_GROWTH;
    gsl_odeiv2_system sys={growth_ode_system,NULL,2,par};
    gsl_odeiv2_driver *d=
      gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkck,0.1*AINIT_GROWTH,0,1E-10);

    y[0]=AINIT_GROWTH;
    y[1]=AINIT_GROWTH*AINIT_GROWTH*AINIT_GROWTH*hubble(par,AINIT_GROWTH);

    int status=gsl_odeiv2_driver_apply(d,&ainit,a,y);
    if(status!=GSL_SUCCESS) {
      fprintf(stderr,"CosmoMad: ODE didn't converge\n");
      exit(1);
    }
    
    *gf=y[0];
    *fg=y[1]/(a*a*hubble(par,a)*y[0]);
  }
}      


static double wind(double x,int setwf)
{
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
  double pk;

  if(lgk<par->spl_pk->x0)
    pk=par->spl_pk->y0*pow(10,par->n_scal*(lgk-par->spl_pk->x0));
  else if(lgk>=par->spl_pk->xf)
    pk=par->spl_pk->yf*pow(10,-3*(lgk-par->spl_pk->xf));
  else
    pk=spline_eval(lgk,par->spl_pk);

  return pk;
}

typedef struct {
  double R;
  ParamGetHI *par;
} SigmaPar;

static double integ_sigma(double logk,void *params)
{
  SigmaPar *sp=(SigmaPar *)params;
  double dum;
  double x,w;
  double kk=pow(10,logk);
  x=kk*sp->R;
  w=wind(x,0);

  dum=TWOPIPIINVLOGTEN*pk_linear0(sp->par,logk)*kk*kk*kk*w*w;

  return dum;
}

static double sigma2(ParamGetHI *par,double R)
{
  gsl_function integrand;
  double relerrt=1E-4;
  double result,error;
  SigmaPar sp;
  sp.par=par;
  sp.R=R;

  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(1000);
  integrand.params=&sp;
  integrand.function=&integ_sigma;
  int stat=gsl_integration_qagil(&integrand,3,0,relerrt,1000,w,
				 &result,&error);
  int_error_handle(stat,result,error);
  gsl_integration_workspace_free(w);

  return result;
}

double z_of_r(ParamGetHI *par,double r) {
  if(r<par->r_arr[0])
    return par->z_arr[0];
  else if(r>=par->r_arr[NZ-1])
    return par->z_arr[NZ-1];
  else {
    int ir=(int)((r-par->r_arr[0])*par->idr);
    double u=(r-par->r_arr[ir])*par->idr;
    if(ir>=NZ-1)
      return par->z_arr[NZ-1];
    else
      return par->z_arr[ir]*(1-u)+par->z_arr[ir+1]*u;
  }
}

double r_of_z(ParamGetHI *par,double z) {
  return spline_eval(z,par->spl_rz);
}

double dgrowth_of_z(ParamGetHI *par,double z) {
  return spline_eval(z,par->spl_dgf);
}

double vgrowth_of_z(ParamGetHI *par,double z) {
  return spline_eval(z,par->spl_vgf);
}

void cosmo_set(ParamGetHI *par)
{
  int ii;
  double z_max,z_min;
  double *x_arr,*y_arr,*yb_arr;

  par->OmegaK=1-par->OmegaM-par->OmegaL;
  growth_factor_and_growth_rate(par,0.,&(par->growth_0),&(par->fgrowth_0));
  par->hubble_0=hubble(par,1.)/HMPC;

  z_min=NU_21/par->pmap->nu_f[par->pmap->n_nu-1]-1;
  z_max=NU_21/par->pmap->nu_0[0]-1;
  par->rmin=chi(par,z_min);
  par->rmax=chi(par,z_max);

  par->l_box=2*(1+1./par->n_grid)*par->rmax;
  par->pos_obs[0]=0.5*par->l_box;
  par->pos_obs[1]=0.5*par->l_box;
  par->pos_obs[2]=0.5*par->l_box;

  x_arr=my_malloc(NZ*sizeof(double));
  y_arr=my_malloc(NZ*sizeof(double));
  yb_arr=my_malloc(NZ*sizeof(double));
  par->r_arr=my_malloc(NZ*sizeof(double));
  par->z_arr=my_malloc(NZ*sizeof(double));

  //Distance-redshift relation
  for(ii=0;ii<NZ;ii++) {
    double z=((double)ii)*DZ;
    x_arr[ii]=z;
    y_arr[ii]=chi(par,z);
  }
  par->spl_rz=spline_init(NZ,y_arr,x_arr,0,x_arr[NZ-1]);
  for(ii=0;ii<NZ;ii++) {
    par->r_arr[ii]=y_arr[0]+(y_arr[NZ-1]-y_arr[0])*(ii+0.5)/NZ;
    par->z_arr[ii]=spline_eval(par->r_arr[ii],par->spl_rz);
  }
  par->idr=NZ/(y_arr[NZ-1]-y_arr[0]);
  spline_free(par->spl_rz);
  par->spl_rz=spline_init(NZ,x_arr,y_arr,0,x_arr[NZ-1]);

  //Growth
  for(ii=0;ii<NZ;ii++) {
    double z=((double)ii)*DZ;
    double hh,gf,fg;
    x_arr[ii]=z;
    growth_factor_and_growth_rate(par,z,&gf,&fg);
    hh=hubble(par,1./(1+z))/HMPC;
    y_arr[ii]=gf/par->growth_0;
    yb_arr[ii]=gf*hh*fg/(par->growth_0*par->fgrowth_0*par->hubble_0);
  }
  par->spl_dgf=spline_init(NZ,x_arr,y_arr,0,y_arr[NZ-1]);
  par->spl_vgf=spline_init(NZ,x_arr,yb_arr,0,yb_arr[NZ-1]);

  free(x_arr);
  free(y_arr);
  free(yb_arr);

  //Pk
  FILE *fpk=my_fopen(par->fnamePk,"r");
  int nk=linecount(fpk); rewind(fpk);
  x_arr=my_malloc(nk*sizeof(double));
  y_arr=my_malloc(nk*sizeof(double));
  for(ii=0;ii<nk;ii++) {
    double kk,pk;
    int stat=fscanf(fpk,"%lf %lf",&kk,&pk);
    if(stat!=2) report_error(1,"Error reading file %s, line %d\n",par->fnamePk,ii+1);
    x_arr[ii]=log10(kk);
    y_arr[ii]=pk;
  }
  fclose(fpk);
  par->spl_pk=spline_init(nk,x_arr,y_arr,0,0);

  double norm_pk=par->sig8*par->sig8/sigma2(par,8.);
  for(ii=0;ii<nk;ii++)
    y_arr[ii]*=norm_pk;
  spline_free(par->spl_pk);
  par->spl_pk=spline_init(nk,x_arr,y_arr,y_arr[0],y_arr[nk-1]);
  
  free(x_arr);
  free(y_arr);
}
