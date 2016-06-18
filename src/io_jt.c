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
#include "common_jt.h"
#include <chealpix.h>

void write_maps(char *prefix_out,int nside,int n_nu,
		flouble **maps)
{
  int ii;
  float *map_aux;

  printf("*** Writing maps %s_###.fits\n",prefix_out);

  map_aux=(float *)my_malloc(nside2npix(nside)*sizeof(float));

  timer(0);
  for(ii=0;ii<n_nu;ii++) {
    long jj;
    char fname[256];
    sprintf(fname,"%s_%03d.fits",prefix_out,ii+1);
    for(jj=0;jj<nside2npix(nside);jj++)
      map_aux[jj]=(float)(maps[ii][jj]);
    he_write_healpix_map(&map_aux,1,nside,fname);
  }
  free(map_aux);

  timer(1);
  printf("\n");
}

static void read_mask(ParamsJoinT *pars)
{
  long ii;

  pars->mask=(int *)my_malloc(nside2npix(pars->nside)*sizeof(int));

  if((!strcmp(pars->fname_mask,"default"))||
     (!strcmp(pars->fname_mask,"none"))) {
    for(ii=0;ii<nside2npix(pars->nside);ii++)
      pars->mask[ii]=1;
    pars->fsky=1;
  }
  else {
    long nside_test,n_theta;
    flouble *mask_f=he_read_healpix_map(pars->fname_mask,&nside_test,0);
    if(nside_test!=pars->nside) {
      fprintf(stderr,"CRIME: read n_side is not the same as provided: %ld != %ld "
	      "in file %s\n",nside_test,(long)(pars->nside),pars->fname_mask);
      exit(1);
    }
    
    n_theta=0;
    for(ii=0;ii<nside2npix(pars->nside);ii++) {
      if((mask_f[ii]<-0.01)||
	 ((mask_f[ii]>0.01)&&(mask_f[ii]<0.99))||
	 (mask_f[ii]>1.01)) {
	fprintf(stderr,"CRIME: Mask must be 0 or 1 : %E\n",mask_f[ii]);
	exit(1);
      }
      if(mask_f[ii]>0.99) {
	pars->mask[ii]=1;
	n_theta++;
      }
      else pars->mask[ii]=0;
    }
    
    pars->fsky=(double)n_theta/nside2npix(pars->nside);
  }
}

ParamsJoinT *read_input_params_JoinT(char *fname)
{
  FILE *fi;
  int nlin,ii;
  int seed_signed=-1;
  ParamsJoinT *pars=set_default_params_JoinT();
  
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
    else if(!strcmp(s1,"fname_mask=")) {
      sprintf(pars->fname_mask,"%s",s2);
      printf("  fname_mask = %s\n",pars->fname_mask);
    }
    else if(!strcmp(s1,"prefix_cosmo=")) {
      sprintf(pars->prefix_cosmo,"%s",s2);
      printf("  prefix_cosmo = %s\n",pars->prefix_cosmo);
    }
    else if(!strcmp(s1,"prefix_synchrotron=")) {
      sprintf(pars->prefix_synchrotron,"%s",s2);
      printf("  prefix_synchrotron = %s\n",pars->prefix_synchrotron);
    }
    else if(!strcmp(s1,"prefix_galactic_freefree=")) {
      sprintf(pars->prefix_galactic_freefree,"%s",s2);
      printf("  prefix_galactic_freefree = %s\n",pars->prefix_galactic_freefree);
    }
    else if(!strcmp(s1,"prefix_extragalactic_freefree=")) {
      sprintf(pars->prefix_extragalactic_freefree,"%s",s2);
      printf("  prefix_extragalactic_freefree = %s\n",
	     pars->prefix_extragalactic_freefree);
    }
    else if(!strcmp(s1,"prefix_point_sources=")) {
      sprintf(pars->prefix_point_sources,"%s",s2);
      printf("  prefix_point_sources = %s\n",pars->prefix_point_sources);
    }
    else if(!strcmp(s1,"prefix_custom=")) {
      sprintf(pars->prefix_custom,"%s",s2);
      printf("  prefix_custom = %s\n",pars->prefix_custom);
    }
    else if(!strcmp(s1,"prefix_out=")) {
      sprintf(pars->prefix_out,"%s",s2);
      printf("  prefix_out = %s\n",pars->prefix_out);
    }
    else if(!strcmp(s1,"nside=")) {
      pars->nside=atoi(s2);
      printf("  nside = %d\n",pars->nside);
    }
    else if(!strcmp(s1,"seed="))
      seed_signed=atoi(s2);
    else if(!strcmp(s1,"do_polarization="))
      pars->do_polarization=atoi(s2);
    else if(!strcmp(s1,"polarization_leakage="))
      pars->polarization_leakage=atof(s2);
    else if(!strcmp(s1,"dish_diameter="))
      pars->dish_diameter=atof(s2);
    else if(!strcmp(s1,"n_dishes="))
      pars->n_dish=atoi(s2);
    else if(!strcmp(s1,"t_system="))
      pars->tsys=atof(s2);
    else if(!strcmp(s1,"time_total="))
      pars->time_tot=atof(s2);
    else {
      fprintf(stderr,"CRIME: unknown param %s in line %d\n",
	      s1,ii+1);
    }
  }

  if(pars->do_polarization!=1) pars->do_polarization=0;
  if(pars->do_polarization) {
    if((pars->polarization_leakage<0)||(pars->polarization_leakage>1)) {
      fprintf(stderr,"CRIME: wrong polarization leakage %.3lE\n",
	      pars->polarization_leakage);
      exit(1);
    }
    else {
      printf("  polarization leakage (%%) : %.3lE\n",100*(pars->polarization_leakage));
    }
  }

  if(seed_signed<=0) pars->seed=time(NULL);
  else pars->seed=(unsigned int)seed_signed;
  printf("  seed : %u\n",pars->seed);
  printf("  Dish diameter (m) : %.3lE\n",pars->dish_diameter);
  printf("  Number of dishes : %d\n",pars->n_dish);
  printf("  System temperature (K) : %.3lE\n",0.001*pars->tsys);
  printf("  Integration time (h) : %.3lE\n",pars->time_tot);

  //SHOULD WE VERIFY THE PRESENCE OF ALL THE FILES?
  pars->nutable=read_nutable(pars->fname_nutable,&(pars->n_nu));

  read_mask(pars);
  printf("  Sky fraction : %.2lE%%\n",100*pars->fsky);

  if(pars->n_nu<=0) {
    fprintf(stderr,"CRIME: Wrong n_nu = %d\n",pars->n_nu);
    exit(1);
  }

  if(pars->nside<=0) {
    fprintf(stderr,"CRIME: Wrong nside = %d\n",pars->nside);
    exit(1);
  }

  printf("\n");
  return pars;
}

double **read_nutable(char *fname,int *n_nu)
{
  int ii,nnu;
  double **nutable=(double **)my_malloc(3*sizeof(double *));
  FILE *fin=fopen(fname,"r");
  if(fin==NULL) error_open_file(fname);
  nnu=linecount(fin);
  rewind(fin);

  *n_nu=nnu;
  for(ii=0;ii<3;ii++) {
    nutable[ii]=(double *)my_malloc(nnu*sizeof(double));
  }

  for(ii=0;ii<nnu;ii++) {
    int inu;
    double nu0,nuf,dum;
    int stat=fscanf(fin,"%d %lf %lf %lf %lf\n",
		    &inu,&nu0,&nuf,&dum,&dum);
    if(stat!=5) error_read_line(fname,ii+1);
    nutable[0][ii]=0.5*(nu0+nuf);
    nutable[1][ii]=nu0;
    nutable[2][ii]=nuf;
  }
  fclose(fin);

  return nutable;
}
