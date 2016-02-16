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

void get_point_sources(ParamGetHI *par)
{
  //////
  // Uses the gaussian matter density field to obtain a
  // poisson sampling of point sources (returns an integer array
  // with the number of sources in each cell).
  lint np_tot=0;

  print_info("*** Getting point sources\n");
  if(NodeThis==0) timer(0);
  print_info("Poisson-sampling\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(np_tot,par,IThread0)
#endif //_HAVE_OMP
  {
    lint iz;
    double dx=par->l_box/par->n_grid;
    double cell_vol=dx*dx*dx;
    lint np_tot_thr=0;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+IThread0+omp_get_thread_num();
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=(ix+0.5)*dx-par->pos_obs[1];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,r);
	  double ndens=n_of_z_psources(par,redshift);
	  int npp;
	  if(ndens>0) {
	    double bias=bias_psources(redshift);
	    double gfb=dgrowth_of_r(par,r)*bias;
	    double lambda=ndens*cell_vol*
	      exp(gfb*(par->grid_dens[index]-0.5*gfb*par->sigma2_gauss));
	    npp=rng_poisson(lambda,rng_thr);
	  }
	  else npp=0;

	  par->nsources[index]=npp;
	  np_tot_thr+=npp;
	}
      }
    }//end omp for

#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      np_tot+=np_tot_thr;
    }//end omp critical
    end_rng(rng_thr);
  }//end omp parallel
  if(NodeThis==0) timer(2);

#ifdef _LONGIDS
  print_info("  There will be %ld particles in total \n",np_tot);
#else //LONGIDS
  print_info("  There will be %d particles in total \n",np_tot);
#endif //_LONGIDS
}

void get_HI(ParamGetHI *par)
{
  //////
  // 1 - Applies lognormal transformation to density grid
  // 2 - Transforms radial velocity into redshift distortion
  // 3 - Substitutes the overdensity grid for the corresponding
  //     HI mass in each cell.
  print_info("*** Gettin' HI\n");

  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par)
#endif //_HAVE_OMP
  {
    lint iz;
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    //Using this the mass grid will be in units of 2.77459457E11*OmegaB/hhub M_sun
    double mass_prefac=dx*dx*dx;

#ifdef _HAVE_OMP
#pragma omp for nowait
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      double z=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      lint iz0=iz*((lint)(ngx*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	double y=dx*(iy+0.5)-par->pos_obs[1];
	int iy0=iy*ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=iz0+iy0+ix;
	  double x=dx*(ix+0.5)-par->pos_obs[0];
	  double r=sqrt(x*x+y*y+z*z);
	  double redshift=z_of_r(par,r);
	  double gfd=dgrowth_of_r(par,r)*bias_HI(redshift);
	  double delta_gauss=par->grid_dens[index];
	  double dz_rsd=par->grid_rvel[index]*vgrowth_of_r(par,r);
	  double dens_LN=exp(gfd*(delta_gauss-0.5*gfd*par->sigma2_gauss));
	  double mass_HI=mass_prefac*fraction_HI(redshift)*dens_LN;
	  par->grid_dens[index]=mass_HI;
	  par->grid_rvel[index]=dz_rsd;
	}
      }
    } //end pragma omp for
  } //end pragma omp parallel
  if(NodeThis==0) timer(2);
  print_info("\n");
}
