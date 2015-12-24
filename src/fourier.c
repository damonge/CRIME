#include "common_gh.h"

static void fftw_wrap(int ng,dftw_complex *pin,flouble *pout)
{
#ifdef _SPREC
  fftwf_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftwf_mpi_plan_dft_c2r_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftwf_plan_dft_c2r_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftwf_execute(plan_ft);
  fftwf_destroy_plan(plan_ft);
#else //_SPREC
  fftw_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftw_mpi_plan_dft_c2r_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftw_plan_dft_c2r_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftw_execute(plan_ft);
  fftw_destroy_plan(plan_ft);
#endif //_SPREC
}

void init_fftw(ParamGetHI *par)
{
  ptrdiff_t dsize;
#ifdef _HAVE_MPI
  ptrdiff_t nz,iz0;

  //Initialize OpenMP fftw (if possible)
#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
    int stat;
#ifdef _SPREC
    stat=fftwf_init_threads();
#else //_SPREC
    stat=fftw_init_threads();
#endif //_SPREC
    if(!stat) {
      fprintf(stderr,"CRIME: Couldn't initialize FFTW threads \n");
      exit(1);
    }
  }
#endif //_HAVE_OMP

  //Initialize MPI fftw
#ifdef _SPREC
  fftwf_mpi_init();
#else //_SPREC
  fftw_mpi_init();
#endif //_SPREC

  //Plan OpenMP fftw (if possible)
#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
#ifdef _SPREC
    fftwf_plan_with_nthreads(omp_get_max_threads());
#else //_SPREC
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif //_SPREC
  }
#endif //_HAVE_OMP

  //Get MPI FFT bounds
#ifdef _SPREC
  dsize=fftwf_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#else //_SPREC
  dsize=fftw_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#endif //_SPREC
  par->nz_here=nz;
  par->iz0_here=iz0;

#else //_HAVE_MPI

#ifdef _HAVE_OMP
  int stat;
#ifdef _SPREC
  stat=fftwf_init_threads();
#else //_SPREC
  stat=fftw_init_threads();
#endif //_SPREC
  if(!stat) {
    fprintf(stderr,"CRIME: Couldn't initialize FFTW threads \n");
    exit(1);
  }
#ifdef _SPREC
  fftwf_plan_with_nthreads(omp_get_max_threads());
#else //_SPREC
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif //_SPREC
#endif //_HAVE_OMP

  dsize=(par->n_grid/2+1)*((lint)(par->n_grid*par->n_grid));
  par->nz_here=par->n_grid;
  par->iz0_here=0;
#endif //_HAVE_MPI

#ifdef _SPREC
  par->grid_dens_f=fftwf_alloc_complex(dsize);
#else //_SPREC
  par->grid_dens_f=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(par->grid_dens_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_dens=(flouble *)(par->grid_dens_f);

#ifdef _SPREC
  par->grid_vpot_f=fftwf_alloc_complex(dsize);
#else //_SPREC
  par->grid_vpot_f=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(par->grid_vpot_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_vpot=(flouble *)(par->grid_vpot_f);

  par->grid_rvel=my_malloc(2*dsize*sizeof(flouble));
#ifdef _HAVE_MPI
  par->slice_left=my_malloc(2*(par->n_grid/2+1)*par->n_grid*sizeof(flouble));
  par->slice_right=my_malloc(2*(par->n_grid/2+1)*par->n_grid*sizeof(flouble));
#endif //_HAVE_MPI
}

void end_fftw(void)
{
#ifdef _HAVE_MPI

#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
#ifdef _SPREC
    fftwf_cleanup_threads();
#else //_SPREC
    fftw_cleanup_threads();
#endif //_SPREC
  }
#endif //_HAVE_OMP

#ifdef _SPREC
  fftwf_mpi_cleanup();
#else //_SPREC
  fftw_mpi_cleanup();
#endif //_SPREC

#else //_HAVE_MPI

#ifdef _HAVE_OMP
#ifdef _SPREC
  fftwf_cleanup_threads();
#else //_SPREC
  fftw_cleanup_threads();
#endif //_SPREC
#endif //_HAVE_OMP

#endif //_HAVE_MPI
}

void create_density_and_radvel(ParamGetHI *par)
{
  print_info(0,"*** Generating Gaussian density and velocity fields\n");

  print_info(1," Generating Fourier-space quantities\n");

  //Create Fourier-space delta and velocity potential
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,IThread0)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double idk3=1./(dk*dk*dk);
#ifdef _HAVE_OMP
    unsigned long seed_thr=par->seed_rng+IThread0+omp_get_thread_num();
#else //_HAVE_OMP
    unsigned long seed_thr=par->seed_rng+IThread0;
#endif //_HAVE_OMP
    Rng *rng_thr=init_rng(seed_thr);
    double factor=par->fgrowth_0*par->hubble_0;
    double norm_factor=pow(sqrt(2*M_PI)/par->l_box,3);
    
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->nz_here;ii++) {
      int jj,ii_true;
      double kz;
      ii_true=par->iz0_here+ii;
      if(2*ii_true<=par->n_grid)
	kz=ii_true*dk;
      else
	kz=-(par->n_grid-ii_true)*dk;
      for(jj=0;jj<par->n_grid;jj++) {
	int kk;
	double ky;
	if(2*jj<=par->n_grid)
	  ky=jj*dk;
	else
	  ky=-(par->n_grid-jj)*dk;
	for(kk=0;kk<=par->n_grid/2;kk++) {
	  double kx;
	  double k_mod2;
	  lint index=kk+(par->n_grid/2+1)*((lint)(jj+par->n_grid*ii)); //Grid index for +k
	  if(2*kk<=par->n_grid)
	    kx=kk*dk;
	  else
	    kx=-(par->n_grid-kk)*dk; //This should never happen
	  
	  k_mod2=kx*kx+ky*ky+kz*kz;
	  
	  if(k_mod2<=0) {
	    par->grid_dens_f[index]=0;
	    par->grid_vpot_f[index]=0;
	  }
	  else {
	    double lgk=0.5*log10(k_mod2);
	    double sigma=sqrt(pk_linear0(par,lgk)*idk3*0.5);
	    double re=sigma*rand_gauss(rng_thr);
	    double im=sigma*rand_gauss(rng_thr);
	    par->grid_dens_f[index]=norm_factor*(re+I*im);
	    par->grid_vpot_f[index]=factor*par->grid_dens_f[index]/k_mod2;
	  }
	}
      }
    } //end omp for
    end_rng(rng_thr);
  } //end omp parallel

  //Transform to real space
  print_info(1," Transforming to real space\n");
  fftw_wrap(par->n_grid,par->grid_dens_f,par->grid_dens);
  fftw_wrap(par->n_grid,par->grid_vpot_f,par->grid_vpot);

  lint slice_size=2*(par->n_grid/2+1)*par->n_grid;
#ifdef _HAVE_MPI
  MPI_Status stat;

  //Pass rightmost slice to right node and receive left slice from left node
  MPI_Sendrecv(&(par->grid_vpot[(par->nz_here-1)*slice_size]),slice_size,FLOUBLE_MPI,NodeRight,1,
	       par->slice_left,slice_size,FLOUBLE_MPI,NodeLeft,1,MPI_COMM_WORLD,&stat);
  //Pass leftmost slice to left node and receive right slice from right node
  MPI_Sendrecv(par->grid_vpot,slice_size,FLOUBLE_MPI,NodeLeft,2,
	       par->slice_right,slice_size,FLOUBLE_MPI,NodeRight,2,MPI_COMM_WORLD,&stat);
#else //_HAVE_MPI
  par->slice_left=&(par->grid_vpot[(par->n_grid-1)*slice_size]);
  par->slice_right=par->grid_vpot;
#endif //_HAVE_MPI

  //Compute radial velocity as projected gradient of the velocity potential
  print_info(1," Differetiating velocity potential\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par)
#endif //_HAVE_OMP
  {
    double dx=par->l_box/par->n_grid;
    double idx=1./dx;
    int iz;
    int ngx=2*(par->n_grid/2+1);
    
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz_hi=iz+1;
      lint iz_lo=iz-1;
      lint iz_0=iz;
      double z=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      iz_hi*=ngx*par->n_grid;
      iz_lo*=ngx*par->n_grid;
      iz_0*=ngx*par->n_grid;
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy_hi=iy+1;
	lint iy_lo=iy-1;
	lint iy_0=iy;
	double y=dx*(iy+0.5)-par->pos_obs[1];
	if(iy==0) iy_lo=par->n_grid-1;
	if(iy==par->n_grid-1) iy_hi=0;
	iy_hi*=ngx;
	iy_lo*=ngx;
	iy_0*=ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  double vel[3],ur[3];
	  lint ix_hi=ix+1;
	  lint ix_lo=ix-1;
	  lint ix_0=ix;
	  double x=dx*(ix+0.5)-par->pos_obs[0];
	  double irr=1./sqrt(x*x+y*y+z*z);
	  if(ix==0) ix_lo=par->n_grid-1;
	  if(ix==par->n_grid-1) ix_hi=0;
	  
	  ur[0]=x*irr;
	  ur[1]=y*irr;
	  ur[2]=z*irr;
	  
	  vel[0]=0.5*idx*(par->grid_vpot[ix_hi+iy_0+iz_0]-par->grid_vpot[ix_lo+iy_0+iz_0]);
	  vel[1]=0.5*idx*(par->grid_vpot[ix_0+iy_hi+iz_0]-par->grid_vpot[ix_0+iy_lo+iz_0]);
	  if(iz==0)
	    vel[2]=0.5*idx*(par->grid_vpot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0]);
	  else if(iz==par->nz_here-1)
	    vel[2]=0.5*idx*(par->slice_right[ix_0+iy_0]-par->grid_vpot[ix_0+iy_0+iz_lo]);
	  else
	    vel[2]=0.5*idx*(par->grid_vpot[ix_0+iy_0+iz_hi]-par->grid_vpot[ix_0+iy_0+iz_lo]);

	  par->grid_rvel[ix_0+iy_0+iz_0]=vel[0]*ur[0]+vel[1]*ur[1]+vel[2]*ur[2];
	}
      }
    } // end omp for
  } // end omp parallel
}
