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

static float **map_aux;

static void write_map(char *fname,int n_side,flouble *map_i,
		      flouble *map_q,flouble *map_u)
{
  int ii;
  int with_pol,nfields;

  if(map_q==NULL) {
    with_pol=0;
    nfields=1;
  }
  else {
    with_pol=1;
    nfields=3;
  }

  for(ii=0;ii<nside2npix(n_side);ii++) {
    map_aux[0][ii]=(float)(map_i[ii]);
    if(with_pol) {
      map_aux[1][ii]=(float)(map_q[ii]);
      map_aux[2][ii]=(float)(map_u[ii]);
    }
  }

  he_write_healpix_map(map_aux,nfields,n_side,fname);
}

void write_maps(char *prefix_out,int nside,int n_nu,
		flouble **maps_i,flouble **maps_q,
		flouble **maps_u)
{
  long ii;
  int with_pol,nfields;
  flouble *map_aux_q,*map_aux_u;

  printf("*** Writing maps %s_###.fits\n",prefix_out);

  if(maps_q==NULL) {
    with_pol=0;
    nfields=1;
    map_aux_q=NULL;
    map_aux_u=NULL;
  }
  else {
    with_pol=1;
    nfields=3;
    map_aux_q=(flouble *)my_malloc(nside2npix(nside)*sizeof(flouble));
    map_aux_u=(flouble *)my_malloc(nside2npix(nside)*sizeof(flouble));
  }
      
  map_aux=(float **)my_malloc(nfields*sizeof(float *));
  for(ii=0;ii<nfields;ii++)
    map_aux[ii]=(float *)my_malloc(nside2npix(nside)*sizeof(float));

  timer(0);
  for(ii=0;ii<n_nu;ii++) {
    char fname[256];

    if(with_pol) {
      he_udgrade(maps_q[ii],NSIDE_POL,map_aux_q,nside,0);
      he_udgrade(maps_u[ii],NSIDE_POL,map_aux_u,nside,0);
    }
    sprintf(fname,"%s_%03ld.fits",prefix_out,ii+1);
    write_map(fname,nside,maps_i[ii],map_aux_q,map_aux_u);
  }

  for(ii=0;ii<nfields;ii++)
    free(map_aux[ii]);
  free(map_aux);

  if(with_pol) {
    free(map_aux_q);
    free(map_aux_u);
  }

  timer(1);
  printf("\n");
}

static void set_cl_model(ParamsForGet *pars)
{
  if(!strcmp(pars->cl_model,"galactic_synchrotron")) {
    pars->do_galaxy=1;
    //    pars->amp=700.0; //This is in SCK
    pars->amp=1600.0; //I think this fits better
    pars->alpha=2.8;
    //    pars->beta=2.4; //This is in SCK
    pars->beta=3.3; //I think this fits better
    pars->xi=4.0;
    pars->nu_ref=130.0;
    pars->lref=1000;
  }
  else if(!strcmp(pars->cl_model,"galactic_freefree")) {
    pars->do_galaxy=0;
    pars->amp=0.088;
    pars->alpha=2.15;
    pars->beta=3.0;
    pars->xi=35.0;
    pars->nu_ref=130.0;
    pars->lref=1000;
  }
  else if(!strcmp(pars->cl_model,"extragalactic_freefree")) {
    pars->do_galaxy=0;
    pars->amp=0.014;
    pars->alpha=2.10;
    pars->beta=1.0;
    pars->xi=35.0;
    pars->nu_ref=130.0;
    pars->lref=1000;
  }
  else if(!strcmp(pars->cl_model,"point_sources")) {
    pars->do_galaxy=0;
    pars->amp=57.0;
    pars->alpha=2.07;
    pars->beta=1.1;
    pars->xi=1.0;
    pars->nu_ref=130.0;
    pars->lref=1000;
  }
  else if(!strcmp(pars->cl_model,"custom")) {
    pars->do_galaxy=0;
  }
  else {
    fprintf(stderr,"CRIME: unknown Cl model %s\n",pars->cl_model);
    exit(1);
  }
  
  if(pars->do_galaxy==0)
    pars->do_polarization=0;
  if(pars->do_polarization!=1) pars->do_polarization=0;

  printf("  cl_model = %s\n",pars->cl_model);
  printf("  amp = %.3lE\n",pars->amp);
  printf("  beta = %.3lE\n",pars->beta);
  printf("  alpha = %.3lE\n",pars->alpha);
  printf("  xi = %.3lE\n",pars->xi);
  printf("  lref = %d\n",pars->lref);
  printf("  nu_ref = %.3lE\n",pars->nu_ref);
  if(pars->do_polarization) {
    printf("  xi_polarization = %.3lE\n",pars->xi_polarization);
    printf("  beta_polarization = %.3lE\n",pars->beta_polarization);
  }
}

ParamsForGet *read_input_params_ForGet(char *fname)
{
  FILE *fi;
  int nlin,ii;
  int seed_signed=-1;
  ParamsForGet *pars=set_default_params_ForGet();
  
  printf("*** Reading run parameters\n");
  sprintf(pars->fname_input,"%s",fname);

  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  nlin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<nlin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"fname_nutable=")) {
      sprintf(pars->fname_nutable,"%s",s2);
      printf("  fname_nutable = %s\n",pars->fname_nutable);
    }
    else if(!strcmp(s1,"fname_haslam=")) {
      sprintf(pars->fname_haslam,"%s",s2);
      printf("  fname_haslam = %s\n",pars->fname_haslam);
    }
    else if(!strcmp(s1,"fname_specin=")) {
      sprintf(pars->fname_specin,"%s",s2);
      printf("  fname_specin = %s\n",pars->fname_specin);
    }
    else if(!strcmp(s1,"prefix_out=")) {
      sprintf(pars->prefix_out,"%s",s2);
      printf("  prefix_out = %s\n",pars->prefix_out);
    }
    else if(!strcmp(s1,"lmin=")) {
      pars->lmin=atoi(s2);
      printf("  lmin = %d\n",pars->lmin);
    }
    else if(!strcmp(s1,"lmax=")) {
      pars->lmax=atoi(s2);
      printf("  lmax = %d\n",pars->lmax);
    }
    else if(!strcmp(s1,"nside=")) {
      pars->nside=atoi(s2);
      printf("  nside = %d\n",pars->nside);
    }
    else if(!strcmp(s1,"seed="))
      seed_signed=atoi(s2);
    else if(!strcmp(s1,"cl_model="))
      sprintf(pars->cl_model,"%s",s2);
    else if(!strcmp(s1,"amp="))
      pars->amp=atof(s2);
    else if(!strcmp(s1,"beta="))
      pars->beta=atof(s2);
    else if(!strcmp(s1,"alpha="))
      pars->alpha=atof(s2);
    else if(!strcmp(s1,"xi="))
      pars->xi=atof(s2);
    else if(!strcmp(s1,"nu_ref="))
      pars->nu_ref=atof(s2);
    else if(!strcmp(s1,"lref="))
      pars->lref=atoi(s2);
    else if(!strcmp(s1,"do_polarization="))
      pars->do_polarization=atoi(s2);
    else if(!strcmp(s1,"xi_polarization="))
      pars->xi_polarization=atof(s2);
    else if(!strcmp(s1,"beta_polarization="))
      pars->beta_polarization=atof(s2);
    else if(!strcmp(s1,"fname_faraday=")) {
      sprintf(pars->fname_faraday,"%s",s2);
      printf("  fname_faraday = %s\n",pars->fname_faraday);
    }
    else {
      fprintf(stderr,"CRIME: unknown param %s in line %d\n",
	      s1,ii+1);
    }
  }

  set_cl_model(pars);

  if(seed_signed<=0) pars->seed=time(NULL);
  else pars->seed=(unsigned int)seed_signed;
  printf("  seed = %u\n",pars->seed);
  
  read_nutable(pars->fname_nutable, pars);
  if(pars->nside<=0) {
    fprintf(stderr,"CRIME: Wrong nside = %d\n",pars->nside);
    exit(1);
  }

  printf("\n");
  return pars;
}

void read_nutable(char *fname, ParamsForGet *pars)
{
  int ii,nnu;
  pars->nutable=(double **)my_malloc(3*sizeof(double *));
  FILE *fin=fopen(fname,"r");
  if(fin==NULL) error_open_file(fname);
  nnu=linecount(fin);
  rewind(fin);

  pars->n_nu=nnu;
  /* We'll just allocate one extra slot in case we need if for
     Haslam freq */
  for(ii=0;ii<3;ii++)
    pars->nutable[ii]=(double *)my_malloc((nnu+1)*sizeof(double)); 

  // only add haslam freq if doing galaxy
  pars->haslam_nu_added=pars->do_galaxy;
  for(ii=0;ii<nnu;ii++) {
    int inu;
    double nu0,nuf,dum;
    int stat=fscanf(fin,"%d %lf %lf %lf %lf\n",
		    &inu,&nu0,&nuf,&dum,&dum);
    if(stat!=5) error_read_line(fname,ii+1);
    pars->nutable[0][ii]=0.5*(nu0+nuf);
    pars->nutable[1][ii]=nu0;
    pars->nutable[2][ii]=nuf;
    if((NU_HASLAM<=nuf)&&(NU_HASLAM>nu0)) pars->haslam_nu_added=0;
  }
  if (pars->haslam_nu_added) {
    pars->n_nu++;
    pars->nutable[0][ii]=NU_HASLAM;
    pars->nutable[1][ii]=NU_HASLAM-0.5; /* dummy freq differences */
    pars->nutable[2][ii]=NU_HASLAM+0.5;
  }
  
  fclose(fin);
}
