#include "common_gh.h"

void write_output(ParamGetHI *par)
{
  int inu;
  flouble *map;
  char fname[256];
  long n_pix_angular=nside2npix(par->pmap->n_side);

  print_info(0,"*** Writing output maps into %s_XXX.fits\n",par->prefix_out);
  for(inu=0;inu<par->pmap->n_nu;inu++) {
    sprintf(fname,"!%s_%03d.fits",par->prefix_out,inu+1);
    map=&(par->maps_HI[inu*n_pix_angular]);
    he_write_healpix_map(&map,1,par->pmap->n_side,fname);
  }
}

static ParamGetHI *param_get_hi_new(void)
{
  ParamGetHI *par=my_malloc(sizeof(ParamGetHI));
  sprintf(par->fnamePk,"default");
  sprintf(par->fname_nulist,"default");
  sprintf(par->prefix_out,"default");
  par->OmegaM=0.3;
  par->OmegaL=0.7;
  par->OmegaB=0.05;
  par->OmegaK=0;
  par->hhub=0.67;
  par->w0=-1.;
  par->wa=0.;
  par->n_scal=0.96;
  par->sig8=0.83;
  par->n_grid=512;
  par->nz_here=par->n_grid;
  par->iz0_here=0;
  par->seed_rng=1234;
  par->do_psources=0;
  par->pmap=my_malloc(sizeof(ParamMaps));
  
  par->growth_0=-1;
  par->fgrowth_0=-1;
  par->hubble_0=-1;
  par->l_box=-1;
  par->r_arr=NULL;
  par->z_arr=NULL;
  par->spl_rz=NULL;
  par->spl_dgf=NULL;
  par->spl_vgf=NULL;
  par->spl_pk=NULL;
  par->grid_dens_f=NULL;
  par->grid_vpot_f=NULL;
  par->grid_dens=NULL;
  par->grid_vpot=NULL;
  par->slice_left=NULL;
  par->slice_right=NULL;
  par->grid_rvel=NULL;
  par->maps_HI=NULL;
  par->maps_PS=NULL;

  return par;
}

static void allocate_maps(ParamGetHI *par)
{
  lint n_pix_total=par->pmap->n_nu*((lint)(12*par->pmap->n_side*par->pmap->n_side));

  par->maps_HI=my_calloc(n_pix_total,sizeof(flouble));
  if(par->do_psources)
    par->maps_PS=my_calloc(n_pix_total,sizeof(flouble));
}

ParamGetHI *read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  ParamGetHI *par=param_get_hi_new();

  //Read parameters from file
  print_info(0,"*** Reading run parameters \n");
  fi=my_fopen(fname,"r");
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);

    if(!strcmp(s1,"prefix_out="))
      sprintf(par->prefix_out,"%s",s2);
    else if(!strcmp(s1,"nu_list="))
      sprintf(par->fname_nulist,"%s",s2);
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
    else if(!strcmp(s1,"w0="))
      par->w0=atof(s2);
    else if(!strcmp(s1,"wa="))
      par->wa=atof(s2);
    else if(!strcmp(s1,"ns="))
      par->n_scal=atof(s2);
    else if(!strcmp(s1,"sigma_8="))
      par->sig8=atof(s2);
    else if(!strcmp(s1,"n_grid="))
      par->n_grid=atoi(s2);
    else if(!strcmp(s1,"n_side="))
      par->pmap->n_side=atoi(s2);
    else if(!strcmp(s1,"seed="))
      par->seed_rng=atoi(s2);
    else if(!strcmp(s1,"do_psources="))
      par->do_psources=atoi(s2);
    else
      report_error(0,"Unknown parameter %s\n",s1);
  }
  fclose(fi);

  print_info(1,"  Reading frequency list\n");
  read_nulist(par->fname_nulist,par->pmap);
  print_info(1,"  Setting cosmology\n");
  cosmo_set(par);
  print_info(1,"  Initializing FFTW\n");
  init_fftw(par);
  print_info(1,"  Allocating big memory\n");
  allocate_maps(par);

  return par;
}

void param_gethi_free(ParamGetHI *par)
{
  free(par->r_arr);
  free(par->z_arr);
  spline_free(par->spl_rz);
  spline_free(par->spl_dgf);
  spline_free(par->spl_vgf);
  spline_free(par->spl_pk);
  param_maps_free(par->pmap);
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
  if(par->do_psources)
    free(par->maps_PS);
  free(par);
  end_fftw();
}
