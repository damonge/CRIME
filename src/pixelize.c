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

#define N_SUBPART 10

#ifdef _IRREGULAR_NUTABLE
static int get_inu(ParamGetHI *par,double nu,int inu_start)
{
  int gotit=0;
  int inu0;
  if(inu_start<0)
    inu0=0;
  else if(inu_start>=par->n_nu)
    inu0=par->n_nu-1;
  else
    inu0=inu_start;
  
  while(!gotit) {
    if((inu0==-1) || (inu0==par->n_nu))
      gotit=1;
    else {
      if(nu<par->nu0_arr[inu0])
	inu0--;
      else {
	if(nu>=par->nuf_arr[inu0])
	  inu0++;
	else
	  gotit=1;
      }
    }
  }

  return inu0;
}
#endif //_IRREGULAR_NUTABLE

void mk_psources_maps(ParamGetHI *par)
{
  int ii;
  double *nu_arr=(double *)my_malloc(par->n_nu*sizeof(double));
  
  for(ii=0;ii<par->n_nu;ii++) {
#ifdef _IRREGULAR_NUTABLE
    nu_arr[ii]=(par->nu0_arr[ii]+par->nuf_arr[ii])*0.5;
#else //_IRREGULAR_NUTABLE
    nu_arr[ii]=par->nu_min+(par->nu_max-par->nu_min)*(ii+0.5)/par->n_nu;
#endif //_IRREGULAR_NUTABLE
  }

  print_info("*** Making source maps\n");

  if(NodeThis==0) timer(0);

  //dbg
#ifdef _HAVE_OMP
#pragma omp parallel default(none)			\
  shared(par,nu_arr,IThread0)
#endif //_HAVE_OMP
  {
    int iz;
    double dx=par->l_box/par->n_grid;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+IThread0+ithr;
    gsl_rng *rng_thr=init_rng(seed_thr);
    long n_pix_ang=nside2npix(par->n_side);
    double dOmega=4*M_PI/n_pix_ang;
    
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      double z0=dx*(iz+par->iz0_here)-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=(lint)(iy*2*(par->n_grid/2+1));
	double y0=dx*iy-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=dx*ix-par->pos_obs[0];
	  int np=par->nsources[index];
	  if(np>0) {
	    int ip;
	    double dz_rsd=(double)(par->grid_rvel[index]);
	    for(ip=0;ip<np;ip++) {
	      double pos[3];
	      double r,red_cosmo,red_true;
	      double l0;
	      int inu;
	      long ipix;
	      pos[0]=x0+dx*rng_01(rng_thr);
	      pos[1]=y0+dx*rng_01(rng_thr);
	      pos[2]=z0+dx*rng_01(rng_thr);
	      r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	      red_cosmo=z_of_r(par,r);
	      red_true=red_cosmo+dz_rsd;
	      vec2pix_ring(par->n_side,pos,&ipix);
	      //Luminosity in units of 10^22 W/Hz
	      l0=draw_luminosity(par,red_cosmo,rng_thr);
	      
	      for(inu=0;inu<par->n_nu;inu++) {
		double nu=nu_arr[inu];
		double temperature=temp_of_l(par,l0,nu,red_true,r,dOmega);
		if(temperature>0) {
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		  par->maps_PS[ipix+n_pix_ang*inu]+=temperature;
		  //		  par->maps_PS[inu+par->n_nu*ipix]+=temperature;
		}
	      }
	    }
	  }
	}
      }
    }//end omp for
  }//end omp parallel
  
  free(nu_arr);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void mk_T_maps(ParamGetHI *par)
{
  int ii;
  double x_sub[N_SUBPART],y_sub[N_SUBPART],z_sub[N_SUBPART];
  //Use this if the mass grid is in units of 2.77459457E11*OmegaB/hhub M_sun
  double M2T_factor=90.057156*par->OmegaB*par->hhub*nside2npix(par->n_side)/(4*M_PI);  //HEY!!!!

  gsl_rng *rng=init_rng(par->seed_rng);
  double lcell=par->l_box/par->n_grid;
  for(ii=0;ii<N_SUBPART;ii++) {
    x_sub[ii]=lcell*(rng_01(rng)-0.5);
    y_sub[ii]=lcell*(rng_01(rng)-0.5);
    z_sub[ii]=lcell*(rng_01(rng)-0.5);
  }
  end_rng(rng);

  print_info("*** Making maps\n");

  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(M2T_factor,x_sub,y_sub,z_sub,par)
#endif //_HAVE_OMP
  {
    long n_pix_ang=nside2npix(par->n_side);
    double dx=par->l_box/par->n_grid;
#ifndef _IRREGULAR_NUTABLE
    double inv_dnu=par->n_nu/(par->nu_max-par->nu_min);
#endif //_IRREGULAR_NUTABLE

    int iz;
#ifdef _HAVE_OMP
#pragma omp single 
#endif //_HAVE_OMP
    {
      print_info(" Collecting masses\n");
    }
#ifdef _HAVE_OMP
#pragma omp for nowait
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      int inu=0;
      double z0=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      lint indexz=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	double y0=dx*(iy+0.5)-par->pos_obs[1];
	lint indexy=(lint)(iy*2*(par->n_grid/2+1));
	for(ix=0;ix<par->n_grid;ix++) {
	  int isub;
	  lint index=ix+indexy+indexz;
	  double x0=dx*(ix+0.5)-par->pos_obs[0];
	  double mass_sub=par->grid_dens[index]/N_SUBPART;
	  //	  double dz_rsd=0;//w.o. RSD
	  double dz_rsd=(double)(par->grid_rvel[index]);//w. RSD
	  for(isub=0;isub<N_SUBPART;isub++) {
	    double x=x0+x_sub[isub];
	    double y=y0+y_sub[isub];
	    double z=z0+z_sub[isub];
	    double r=sqrt(x*x+y*y+z*z);
	    double redshift=z_of_r(par,r)+dz_rsd;
	    double nu=NU_21/(1+redshift);
#ifdef _IRREGULAR_NUTABLE
	    inu=get_inu(par,nu,inu);
#else //_IRREGULAR_NUTABLE
	    inu=(int)(inv_dnu*(nu-par->nu_min));
#endif //_IRREGULAR_NUTABLE
	    if((inu>=0)&&(inu<par->n_nu)) {
	      long ipix;
	      double pos[3]={x,y,z};
	    
	      vec2pix_ring(par->n_side,pos,&ipix);
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
	      par->maps_HI[ipix+n_pix_ang*inu]+=mass_sub;
	    }
	  }
	}
      }
    } //end pragma omp for
#ifdef _HAVE_OMP
#pragma omp barrier
#endif //_HAVE_OMP

    int inu;
#ifdef _HAVE_OMP
#pragma omp single 
#endif //_HAVE_OMP
    {
      print_info(" Normalizing to temperature\n");
    }
#ifdef _HAVE_OMP
#pragma omp for nowait
#endif //_HAVE_OMP
    for(inu=0;inu<par->n_nu;inu++) {
      long ipix;
#ifdef _IRREGULAR_NUTABLE
      double dnu=par->nuf_arr[inu]-par->nu0_arr[inu];
      double nu=(par->nuf_arr[inu]+par->nu0_arr[inu])*0.5;
#else //_IRREGULAR_NUTABLE
      double dnu=(par->nu_max-par->nu_min)/par->n_nu;
      double nu=par->nu_min+(inu+0.5)*dnu;
#endif //_IRREGULAR_NUTABLE
      double r=r_of_z(par,NU_21/nu-1);
      double prefac=M2T_factor/(r*r*dnu);
      for(ipix=0;ipix<n_pix_ang;ipix++) {
	long index=ipix+n_pix_ang*inu;
	par->maps_HI[index]*=prefac;
      }
    } //end pragma omp for
  } //end pragma omp parallel
  if(NodeThis==0) timer(2);
  print_info("\n");

#ifdef _HAVE_MPI
#define REDUCE_BATCH 1073741824
  int remainder;
  long ipix0_here=0;
  long npix_all=par->n_nu*nside2npix(par->n_side);
  while(ipix0_here+REDUCE_BATCH<=npix_all) {
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,&(par->maps_HI[ipix0_here]),REDUCE_BATCH,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(&(par->maps_HI[ipix0_here]),NULL,REDUCE_BATCH,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    ipix0_here+=REDUCE_BATCH;
  }
  remainder=(int)(npix_all-ipix0_here);
  if(remainder) {
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,&(par->maps_HI[ipix0_here]),remainder,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(&(par->maps_HI[ipix0_here]),NULL,remainder,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
#endif //_HAVE_MPI
}
