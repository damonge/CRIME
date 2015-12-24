#include "common_gh.h"

int main(int argc,char **argv)
{
  ParamGetHI *par;
  char fname_init[256];

  if(argc!=2) {
    fprintf(stderr,"Usage: ./GetHI param_file\n");
    exit(1);
  }
  sprintf(fname_init,"%s",argv[1]);

  mpi_init(&argc,&argv);

  par=read_run_params(fname_init);
  int ii;
  int nz=1000;
  FILE *fi=my_fopen("test_bg.txt","w");
  for(ii=0;ii<nz;ii++) {
    double z=5.*(ii+0.5)/nz;
    fprintf(fi,"%lE %lE \n",z,z-z_of_r(par,r_of_z(par,z)));
  }
  fclose(fi);
  create_density_and_radvel(par);
  mk_HI_maps(par);
  if(NodeThis==0)
    write_output(par);
  param_gethi_free(par);

#ifdef _HAVE_MPI
  MPI_Finalize();
#endif //_HAVE_MPI
  return 0;
}
