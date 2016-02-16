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
#include <chealpix.h>
#include <fitsio.h>

static float *T_map_single;

#ifdef _IRREGULAR_NUTABLE
static void read_nutable(ParamGetHI *par)
{
  int inu;
  FILE *fi;
  fi=fopen(par->fnameNuTable,"r");
  if(fi==NULL) error_open_file(par->fnameNuTable);
  par->n_nu=linecount(fi)-1;
  rewind(fi);
  par->nu0_arr=my_malloc(par->n_nu*sizeof(double));
  par->nuf_arr=my_malloc(par->n_nu*sizeof(double));
  for(inu=0;inu<=par->n_nu;inu++) {
    int stat;
    double nu;
    stat=fscanf(fi,"%lf ",&nu);
    if(stat!=1)
      report_error(1,"Error reading file %s, line %d\n",par->fnameNuTable,inu+1);
    if(inu!=par->n_nu)
      par->nu0_arr[inu]=nu;
    if(inu!=0)
      par->nuf_arr[inu-1]=nu;
  }
  for(inu=0;inu<par->n_nu;inu++) {
    if(par->nuf_arr[inu]<=par->nu0_arr[inu])
      report_error(1,"Frequency bins don't make sense\n");
  }

  par->nu_max=par->nuf_arr[par->n_nu-1];
  par->nu_min=par->nu0_arr[0];
}
#endif //_IRREGULAR_NUTABLE

static void allocate_maps(ParamGetHI *par)
{
  lint n_pix_total=par->n_nu*((lint)(12*par->n_side*par->n_side));

  par->maps_HI=my_calloc(n_pix_total,sizeof(flouble));
  if(par->do_psources)
    par->maps_PS=my_calloc(n_pix_total,sizeof(flouble));
}

static void write_single_fit(char *prefix,int i_nu,flouble *T_maps,long n_side)
{
  char fname[256];
  long n_pix_ang=nside2npix(n_side);
  long npix_start=i_nu*n_pix_ang;

  sprintf(fname,"%s_%03d.fits",prefix,i_nu+1);

  long ii;
  for(ii=0;ii<n_pix_ang;ii++)
    T_map_single[ii]=(float)(T_maps[npix_start+ii]);

  he_write_healpix_map(&T_map_single,1,n_side,fname);
}

static void write_nu_table(ParamGetHI *par,char *prefix)
{
  int ii;
  char fname_nutable[256];
  FILE *ftab;

  sprintf(fname_nutable,"%s_nuTable.dat",prefix);
  ftab=fopen(fname_nutable,"w");
  if(ftab==NULL) error_open_file(fname_nutable);
  for(ii=0;ii<par->n_nu;ii++) {
#ifdef _IRREGULAR_NUTABLE
    double nu0=par->nu0_arr[ii];
    double nuf=par->nuf_arr[ii];
#else //_IRREGULAR_NUTABLE
    double nu0=par->nu_min+(par->nu_max-par->nu_min)*(ii+0.0)/par->n_nu;
    double nuf=par->nu_min+(par->nu_max-par->nu_min)*(ii+1.0)/par->n_nu;
#endif //_IRREGULAR_NUTABLE
    double z0=NU_21/nuf-1;
    double zf=NU_21/nu0-1;
    
    fprintf(ftab,"%d %lf %lf %lf %lf\n",ii+1,nu0,nuf,z0,zf);
  }
  fclose(ftab);
}

void write_maps(ParamGetHI *par)
{
  int i_nu;
#ifdef _DEBUG
  write_nu_table(par,par->prefixOut);
#endif //_DEBUG
  
  T_map_single=(float *)my_malloc(nside2npix(par->n_side)*sizeof(float));
  
  print_info("*** Writing files %s_###.fits\n",par->prefixOut);
  for(i_nu=0;i_nu<par->n_nu;i_nu++)
    write_single_fit(par->prefixOut,i_nu,par->maps_HI,par->n_side);
  
  if(par->do_psources) {
    char prefix_psources[256];
    sprintf(prefix_psources,"%s_ps",par->prefixOut);
    print_info("*** Writing single slice files %s_###.fits\n",prefix_psources);
    for(i_nu=0;i_nu<par->n_nu;i_nu++)
      write_single_fit(prefix_psources,i_nu,par->maps_PS,par->n_side);
  }
  
  free(T_map_single);
}

static ParamGetHI *param_gethi_new(void)
{
  ParamGetHI *par=my_malloc(sizeof(ParamGetHI));

  sprintf(par->fnamePk,"default");
  par->OmegaM=0.3;
  par->OmegaL=0.7;
  par->OmegaB=0.05;
  par->hhub=0.7;
  par->weos=-1.;
  par->n_scal=0.96;
  par->sig8=0.83;
  par->fgrowth_0=-1;
  par->hubble_0=-1;
  par->z_max=1.5;
  par->z_min=0.5;
  par->r_max=-1;
  par->r_min=-1;
  par->r2_smooth=2.0;
  par->do_smoothing=1;
  par->numk=0;
  par->logkmax=1;
  par->logkmin=-3;
  par->idlogk=100;
  par->logkarr=NULL;
  par->pkarr=NULL;
  par->glob_idr=-1;
  par->seed_rng=1234;
  par->n_side=128;
  par->nu_max=1050.;
  par->nu_min=350.;
  par->n_nu=150;
#ifdef _IRREGULAR_NUTABLE
  sprintf(par->fnameNuTable,"default");
  par->nu0_arr=NULL;
  par->nuf_arr=NULL;
#endif //_IRREGULAR_NUTABLE
  par->n_grid=512;
  par->l_box=-1;
  par->nz_here=512;
  par->iz0_here=0;
  sprintf(par->prefixOut,"default");
  par->grid_dens_f=NULL;
  par->grid_dens=NULL;
  par->grid_vpot_f=NULL;
  par->grid_vpot=NULL;
  par->grid_rvel=NULL;
  par->sigma2_gauss=-1;
  par->nsources=NULL;
  par->maps_HI=NULL;
  par->maps_PS=NULL;

  return par;
}

ParamGetHI *read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  ParamGetHI *par=param_gethi_new();
  
  //Read parameters from file
  print_info("*** Reading run parameters \n");
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"prefix_out="))
      sprintf(par->prefixOut,"%s",s2);
    else if(!strcmp(s1,"pk_filename="))
      sprintf(par->fnamePk,"%s",s2);
    else if(!strcmp(s1,"omega_M="))
      par->OmegaM=atof(s2);
    else if(!strcmp(s1,"omega_L="))
      par->OmegaL=atof(s2);
    else if(!strcmp(s1,"omega_B="))
      par->OmegaB=atof(s2);
    else if(!strcmp(s1,"h="))
      par->hhub=atof(s2);
    else if(!strcmp(s1,"w="))
      par->weos=atof(s2);
    else if(!strcmp(s1,"ns="))
      par->n_scal=atof(s2);
    else if(!strcmp(s1,"sigma_8="))
      par->sig8=atof(s2);
    else if(!strcmp(s1,"r_smooth="))
      par->r2_smooth=atof(s2);
#ifdef _IRREGULAR_NUTABLE
    else if(!strcmp(s1,"frequencies_filename="))
      sprintf(par->fnameNuTable,"%s",s2);
#else //_IRREGULAR_NUTABLE
    else if(!strcmp(s1,"nu_min="))
      par->nu_min=atof(s2);
    else if(!strcmp(s1,"nu_max="))
      par->nu_max=atof(s2);
    else if(!strcmp(s1,"n_nu="))
      par->n_nu=atoi(s2);
#endif //_IRREGULAR_NUTABLE
    else if(!strcmp(s1,"n_grid="))
      par->n_grid=atoi(s2);
    else if(!strcmp(s1,"n_side="))
      par->n_side=atoi(s2);
    else if(!strcmp(s1,"seed="))
      par->seed_rng=atoi(s2);
    else if(!strcmp(s1,"do_psources="))
      par->do_psources=atoi(s2);
    else
      fprintf(stderr,"CRIME: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  if(par->r2_smooth>0) {
    par->r2_smooth=pow(par->r2_smooth,2);
    par->do_smoothing=1;
  }
  else
    par->do_smoothing=0;
#ifdef _IRREGULAR_NUTABLE
  read_nutable(par);
#endif //_IRREGULAR_NUTABLE
  cosmo_set(par);
  init_fftw(par);
  allocate_maps(par);

  double dk=2*M_PI/par->l_box;
  double dtheta=sqrt(41253./(12*par->n_side*par->n_side));
  double drTmax=par->r_max*dtheta/RTOD;
  double drTmin=par->r_min*dtheta/RTOD;
  double dnu=(par->nu_max-par->nu_min)/par->n_nu;
  double drL=(par->r_max-par->r_min)/par->n_nu;
  double output_size=(flouble)(12*par->n_side*par->n_side*par->n_nu*sizeof(flouble))/(1024*1024*1024);
  print_info("Run parameters: \n");
  print_info("  %.3lf < nu/MHz < %.3lf\n",par->nu_min,par->nu_max);
  print_info("  %.3lf < z < %.3lf\n",par->z_min,par->z_max);
  print_info("  %.3lf < r/(Mpc/h) < %.3lf\n",par->r_min,par->r_max);
  print_info("  L_box = %.3lf Mpc/h, N_grid = %d \n",par->l_box,par->n_grid);
  print_info("  Scales resolved: %.3lE < k < %.3lE h/Mpc\n",dk,0.5*(par->n_grid-1)*dk);
  print_info("  Fourier-space resolution: dk = %.3lE h/Mpc\n",dk);
  print_info("  Real-space resolution: dx = %.3lE Mpc/h\n",par->l_box/par->n_grid);
  if(par->do_smoothing)
    print_info("  Density field pre-smoothed on scales: x_s = %.3lE Mpc/h\n",sqrt(par->r2_smooth));
  else
    print_info("  No extra smoothing\n");
  print_info("  n_nu = %d, d_nu= %.3lf MHz, drL ~ %.3lf Mpc/h\n",
	     par->n_nu,dnu,drL);
  print_info("  n_side = %ld, dtheta = %.3lf deg, %.3lf < drT/(Mpc/h) < %.3lf\n",
	     par->n_side,dtheta,drTmin,drTmax);
  print_info("  Estimated output size ~ %.1lf GB \n",output_size);

  print_info("\n");
  
  return par;
}

void param_gethi_free(ParamGetHI *par)
{
  free(par->logkarr);
  free(par->pkarr);
#ifdef _SPREC
  fftwf_free(par->grid_dens_f);
  fftwf_free(par->grid_vpot_f);
#else //_SPREC
  fftw_free(par->grid_dens_f);
  fftw_free(par->grid_vpot_f);
#endif //_SPREC
#ifdef _HAVE_MPI
  free(par->slice_left);
  free(par->slice_right);
#endif //_HAVE_MPI
  free(par->grid_rvel);
  free(par->maps_HI);
#ifdef _IRREGULAR_NUTABLE
  free(par->nu0_arr);
  free(par->nuf_arr);
#endif //_IRREGULAR_NUTABLE
  if(par->do_psources) {
    free(par->nsources);
    free(par->maps_PS);
  }
  end_fftw();
}
