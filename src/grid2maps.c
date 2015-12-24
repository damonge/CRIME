#include "common_gh.h"

double compute_sigma_dens(ParamGetHI *par)
{
  double sigma2_gauss=0;
  double mean_gauss=0;

  //Compute Gaussian variance
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,sigma2_gauss,mean_gauss)
#endif //_HAVE_OMP
  {
    int iz;
    double sigma2_thr=0;
    double mean_thr=0;
    lint ng_tot=par->n_grid*((lint)(par->n_grid*par->n_grid));

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz0=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy0=iy*2*(par->n_grid/2+1);
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+iy0+iz0;
	  sigma2_thr+=par->grid_dens[index]*par->grid_dens[index];
	  mean_thr+=par->grid_dens[index];
	}
      }
    } //end omp for
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      mean_gauss+=mean_thr/ng_tot;
      sigma2_gauss+=sigma2_thr/ng_tot;
    } //end omp critical
  } //end omp parallel

#ifdef _HAVE_MPI
  double sigma2_gathered=0,mean_gathered=0;

  MPI_Allreduce(&sigma2_gauss,&sigma2_gathered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&mean_gauss,&mean_gathered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  mean_gauss=mean_gathered;
  sigma2_gauss=sigma2_gathered;
#endif //_HAVE_MPI

  sigma2_gauss-=mean_gauss*mean_gauss;
  print_info(1," <d>=%.3lE, <d^2>=%.3lE\n",mean_gauss,sqrt(sigma2_gauss));
  
  return sigma2_gauss;
}

void mk_HI_maps(ParamGetHI *par)
{
  double sigma2_gauss;
  int n_pix_ang=nside2npix(par->pmap->n_side);

  print_info(0,"*** Computing HI density and interpolating to pixelized maps\n");
  print_info(1," Computing Gaussian variance\n");
  sigma2_gauss=compute_sigma_dens(par);

  print_info(1," Computing HI mass and interpolating to skymaps\n");
  //Compute lognormal density, redshift distortion, HI mass and interpolate to pixels
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,sigma2_gauss,IThread0,n_pix_ang)
#endif //_HAVE_OMP
  {
    int iz;
    double *x_sub,*y_sub,*z_sub;
    double dx=par->l_box/par->n_grid;
#ifdef _HAVE_OMP
    unsigned long seed_thr=par->seed_rng+IThread0+omp_get_thread_num();
#else //_HAVE_OMP
    unsigned long seed_thr=par->seed_rng;
#endif //_HAVE_OMP
    Rng *rng_thr=init_rng(seed_thr);
    x_sub=my_malloc(N_SUBPART_HI*sizeof(double));
    y_sub=my_malloc(N_SUBPART_HI*sizeof(double));
    z_sub=my_malloc(N_SUBPART_HI*sizeof(double));
    for(iz=0;iz<N_SUBPART_HI;iz++) {
      x_sub[iz]=dx*(rand_real01(rng_thr)-0.5);
      y_sub[iz]=dx*(rand_real01(rng_thr)-0.5);
      z_sub[iz]=dx*(rand_real01(rng_thr)-0.5);
    }
    end_rng(rng_thr);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz0=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      double z=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy0=iy*2*(par->n_grid/2+1);
	double y=dx*(iy+0.5)-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int isub;
	  int inu=1;
	  lint index=ix+iy0+iz0;
	  double x=dx*(ix+0.5)-par->pos_obs[0];
	  double r=sqrt(x*x+y*y+z*z);
	  if((r>0.9*par->rmin) && (r<1.1*par->rmax)) {
	    double redshift=z_of_r(par,r);
	    double dgrowth=dgrowth_of_z(par,redshift)*bias_HI(redshift);
	    double vgrowth=vgrowth_of_z(par,redshift);
	    double hi_frac=fraction_HI(redshift);
	    double dz_rsd=par->grid_rvel[index]*vgrowth;
	    double dens_ln=exp(dgrowth*(par->grid_dens[index]-0.5*dgrowth*sigma2_gauss));
	    double mass_HI_sub=dx*dx*dx*hi_frac*dens_ln/N_SUBPART_HI;
	    for(isub=0;isub<N_SUBPART_HI;isub++) {
	      double xx=x+x_sub[isub];
	      double yy=y+y_sub[isub];
	      double zz=z+z_sub[isub];
	      double rr=sqrt(xx*xx+yy*yy+zz*zz);
	      double rred=z_of_r(par,rr)+dz_rsd;
	      double nu=NU_21/(1+rred);
	      inu=get_inu(par->pmap,nu,inu-1);
	      if((inu>0)&&(inu<par->pmap->n_nu)) {
		long ipix;
		double pos[3]={xx,yy,zz};
		
		vec2pix_ring(par->pmap->n_side,pos,&ipix);
#pragma omp atomic
		par->maps_HI[ipix+n_pix_ang*inu]+=mass_HI_sub;
	      }
	    }
	  }
	}
      }
    } //end omp for
    free(x_sub);
    free(y_sub);
    free(z_sub);
  }//end omp parallel

  print_info(1," Normalizing skymaps\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,n_pix_ang)
#endif //_HAVE_OMP
  {
    int inu;
    double m2t_factor=90.057156*par->OmegaB*par->hhub*n_pix_ang/(4*M_PI);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(inu=0;inu<par->pmap->n_nu;inu++) {
      long ipix;
      double dnu=par->pmap->nu_f[inu]-par->pmap->nu_0[inu];
      double nu=0.5*(par->pmap->nu_f[inu]+par->pmap->nu_0[inu]);
      double r=r_of_z(par,NU_21/nu-1);
      double prefac=m2t_factor/(dnu*r*r);
      for(ipix=0;ipix<n_pix_ang;ipix++)
	par->maps_HI[ipix+n_pix_ang*inu]*=prefac;
    } //end omp for
  } //end omp parallel

#ifdef _HAVE_MPI
  if(NodeThis==0) {
    MPI_Reduce(MPI_IN_PLACE,par->maps_HI,n_pix_ang*par->pmap->n_nu,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(par->maps_HI,NULL,n_pix_ang*par->pmap->n_nu,
	       FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
#endif //_HAVE_MPI
}
