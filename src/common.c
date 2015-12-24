#include "common.h"

int NodeThis=0;
int NodeLeft=0;
int NodeRight=0;
int NNodes=1;
int IThread0=0;
int MPIThreadsOK=0;

void mpi_init(int* p_argc,char*** p_argv)
{
#ifdef _HAVE_MPI
  int ii,nthreads_this;
  int *nthreads_all;
#ifdef _HAVE_OMP
  int provided;
  MPI_Init_thread(p_argc,p_argv,MPI_THREAD_FUNNELED,&provided);
  MPIThreadsOK = provided >= MPI_THREAD_FUNNELED;
#else //_HAVE_OMP
  MPI_Init(p_argc,p_argv);
  MPIThreadsOK=0;
#endif //_HAVE_OMP

  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);
  if(NodeThis==0)
    NodeLeft=NNodes-1;
  else
    NodeLeft=NodeThis-1;
  if(NodeThis==NNodes-1)
    NodeRight=0;
  else
    NodeRight=NodeThis+1;
     

  nthreads_all=my_malloc(NNodes*sizeof(int));
#ifdef _HAVE_OMP
  nthreads_this=omp_get_max_threads();
#else //_HAVE_OMP
  nthreads_this=1;
#endif //_HAVE_OMP
  MPI_Allgather(&nthreads_this,1,MPI_INT,nthreads_all,1,MPI_INT,MPI_COMM_WORLD);
  if(NodeThis==0) {
    for(ii=0;ii<NNodes;ii++)
      printf("Node %d has %d threads\n",ii,nthreads_all[ii]);
  }
  IThread0=0;
  for(ii=0;ii<NodeThis;ii++)
    IThread0+=nthreads_all[ii];
  free(nthreads_all);

#else //_HAVE_MPI
  NodeThis=0;
  NNodes=1;
  IThread0=0;
#ifdef _HAVE_OMP
  MPIThreadsOK=1;
#else //_HAVE_OMP
  MPIThreadsOK=0;
#endif //_HAVE_OMP
#endif //_HAVE_MPI
#ifdef _DEBUG
  printf("Node %d, thread count starts at %d\n",NodeThis,IThread0);
#endif //_DEBUG
}

void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr,"Node %d, Fatal: %s",NodeThis,msg);
    exit(level);
  }
  else {
    fprintf(stderr,"Node %d, Warning: %s",NodeThis,msg);
  }
}

void print_info(int level,char *fmt,...)
{
  if(NodeThis==0) {
    if(level<=VERBOSITY) {
      va_list args;
      char msg[256];
      
      va_start(args,fmt);
      vsprintf(msg,fmt,args);
      va_end(args);
      
      printf("%s",msg);
    }
  }
}

void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

size_t my_fread(void *ptr,size_t size,size_t nmemb,FILE *stream)
{
  if(fread(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error freading\n");

  return nmemb;
}

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL)
    report_error(1,"Out of memory\n");

  return outptr;
}

FILE *my_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL)
    report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

int linecount(FILE *f)
{
  int i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

SplPar *spline_init(int n,double *x,double *y,double y0,double yf)
{
  SplPar *spl=(SplPar *)my_malloc(sizeof(SplPar));

  spl->x0=x[0];
  spl->xf=x[n-1];
  spl->y0=y0;
  spl->yf=yf;

  spl->intacc=gsl_interp_accel_alloc();
  spl->spline=gsl_spline_alloc(gsl_interp_cspline,n);
  //  spl->spline=gsl_spline_alloc(gsl_interp_linear,n);
  gsl_spline_init(spl->spline,x,y,n);

  return spl;
}

double spline_eval(double x,SplPar *spl)
{
  if(x<=spl->x0)
    return spl->y0;
  else if(x>=spl->xf) 
    return spl->yf;
  else
    return gsl_spline_eval(spl->spline,x,spl->intacc);
}

void spline_free(SplPar *spl)
{
  gsl_spline_free(spl->spline);
  gsl_interp_accel_free(spl->intacc);
  free(spl);
}

void read_nulist(char *fname,ParamMaps *pmap)
{
  int ii;
  FILE *fi=my_fopen(fname,"r");
  pmap->n_nu=linecount(fi);
  rewind(fi);
  pmap->nu_0=my_malloc(pmap->n_nu*sizeof(double));
  pmap->nu_f=my_malloc(pmap->n_nu*sizeof(double));
  for(ii=0;ii<pmap->n_nu;ii++) {
    int idum;
    int stat=fscanf(fi,"%d %lf %lf\n",&idum,&(pmap->nu_0[ii]),&(pmap->nu_f[ii]));
    if(stat!=3)
      report_error(1,"Error reading file %s, line %d\n",fname,ii+1);
  }
  fclose(fi);
}

int get_inu(ParamMaps *pmap,double nu,int inu_start)
{
  int gotit=0;
  int inu0;

  if(inu_start<0)
    inu0=0;
  else if(inu_start>=pmap->n_nu)
    inu0=pmap->n_nu-1;
  else
    inu0=inu_start;

  while(!gotit) {
    if((inu0==-1) || (inu0==pmap->n_nu))
      gotit=1;
    else {
      if(nu<pmap->nu_0[inu0])
	inu0--;
      else {
	if(nu>=pmap->nu_f[inu0])
	  inu0++;
	else
	  gotit=1;
      }
    }
  }

  return inu0;
}

void param_maps_free(ParamMaps *pmap)
{
  free(pmap->nu_0);
  free(pmap->nu_f);
  free(pmap);
}
