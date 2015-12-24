#include "common.h"
#include <fitsio.h>
#include <chealpix.h>
#include "sharp_almhelpers.h"
#include "sharp_geomhelpers.h"
#include "sharp.h"

#define MAX_SHT 32

long he_nalms(int lmax)
{
  return ((lmax+1)*(lmax+2))/2;
}

long he_indexlm(int l,int m,int lmax)
{
  if(m>0)
    return l+m*lmax-(m*(m-1))/2;
  else
    return l;
}

void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms)
{
  int nbatches,nodd,itrans;
  double time;
  sharp_alm_info *alm_info;
  sharp_geom_info *geom_info;

  sharp_make_triangular_alm_info(lmax,lmax,1,&alm_info);
  sharp_make_weighted_healpix_geom_info(nside,1,NULL,&geom_info);
  
  nbatches=ntrans/MAX_SHT;
  nodd=ntrans%MAX_SHT;

  for(itrans=0;itrans<nbatches;itrans++) {
    time=0;
#ifdef _SPREC
    sharp_execute(SHARP_ALM2MAP,0,&(alms[itrans*MAX_SHT]),
                  &(maps[itrans*MAX_SHT]),geom_info,
                  alm_info,MAX_SHT,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_ALM2MAP,0,&(alms[itrans*MAX_SHT]),
                  &(maps[itrans*MAX_SHT]),geom_info,
                  alm_info,MAX_SHT,SHARP_DP,&time,NULL);
#endif //_SPREC
    //  printf("  Took %lf s according to libsharp\n",time);
  }
  if(nodd>0) {
    time=0;
#ifdef _SPREC
    sharp_execute(SHARP_ALM2MAP,0,&(alms[nbatches*MAX_SHT]),
		  &(maps[nbatches*MAX_SHT]),geom_info,
		  alm_info,nodd,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_ALM2MAP,0,&(alms[nbatches*MAX_SHT]),
		  &(maps[nbatches*MAX_SHT]),geom_info,
		  alm_info,nodd,SHARP_DP,&time,NULL);
#endif //_SPREC
    //  printf("  Took %lf s according to libsharp\n",time);
  }

  sharp_destroy_geom_info(geom_info);
  sharp_destroy_alm_info(alm_info);
}

void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms)
{
  int nbatches,nodd,itrans;
  double time;
  sharp_alm_info *alm_info;
  sharp_geom_info *geom_info;

  sharp_make_triangular_alm_info(lmax,lmax,1,&alm_info);
  sharp_make_weighted_healpix_geom_info(nside,1,NULL,&geom_info);

  nbatches=ntrans/MAX_SHT;
  nodd=ntrans%MAX_SHT;

  for(itrans=0;itrans<nbatches;itrans++) {
    time=0;
#ifdef _SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms[itrans*MAX_SHT]),
                  &(maps[itrans*MAX_SHT]),geom_info,
                  alm_info,MAX_SHT,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms[itrans*MAX_SHT]),
                  &(maps[itrans*MAX_SHT]),geom_info,
                  alm_info,MAX_SHT,SHARP_DP,&time,NULL);
#endif //_SPREC
    //  printf("  Took %lf s according to libsharp\n",time);
  }
  if(nodd>0) {
    time=0;
#ifdef _SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms[nbatches*MAX_SHT]),
		  &(maps[nbatches*MAX_SHT]),geom_info,
		  alm_info,nodd,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms[nbatches*MAX_SHT]),
		  &(maps[nbatches*MAX_SHT]),geom_info,
		  alm_info,nodd,SHARP_DP,&time,NULL);
#endif //_SPREC
    //  printf("  Took %lf s according to libsharp\n",time);
  }

  sharp_destroy_geom_info(geom_info);
  sharp_destroy_alm_info(alm_info);
}

void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,
		int pol_1,int pol_2,
		flouble **cls,int nside,int lmax)
{
  double time;
  sharp_alm_info *alm_info;
  sharp_geom_info *geom_info;
  fcomplex **alms_1,**alms_2;
  int i1,index_cl;
  int lmax_here=3*nside-1;


  alms_1=my_malloc(nmaps_1*sizeof(fcomplex *));
  for(i1=0;i1<nmaps_1;i1++)
    alms_1[i1]=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));
  if(pol_1) {
    if(nmaps_1!=2)
      report_error(1,"Must provide 2 maps for polarization\n");

    sharp_make_triangular_alm_info(lmax,lmax,1,&alm_info);
    sharp_make_weighted_healpix_geom_info(nside,1,NULL,&geom_info);

    /*
    //Transform T
    time=0;
#ifdef _SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms_1[0]),&(maps_1[0]),geom_info,
		  alm_info,1,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_MAP2ALM,0,&(alms_1[0]),&(maps_1[0]),geom_info,
		  alm_info,1,SHARP_DP,&time,NULL);
#endif //_SPREC
    */
    //Transform Q,U
#ifdef _SPREC
    sharp_execute(SHARP_MAP2ALM,2,&(alms_1[0]),&(maps_1[0]),geom_info,
		  alm_info,1,0,&time,NULL);
#else //_SPREC
    sharp_execute(SHARP_MAP2ALM,2,&(alms_1[0]),&(maps_1[0]),geom_info,
		  alm_info,1,SHARP_DP,&time,NULL);
#endif //_SPREC

    sharp_destroy_geom_info(geom_info);
    sharp_destroy_alm_info(alm_info);
  }
  else {
    he_map2alm(nside,lmax,nmaps_1,maps_1,alms_1);
  }

  if(maps_1==maps_2)
    alms_2=alms_1;
  else {
    alms_2=my_malloc(nmaps_2*sizeof(fcomplex *));
    for(i1=0;i1<nmaps_2;i1++)
      alms_2[i1]=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));
    if(pol_2) {
      if(nmaps_2!=2)
	report_error(1,"Must provide 2 maps for polarization\n");
      
      sharp_make_triangular_alm_info(lmax,lmax,1,&alm_info);
      sharp_make_weighted_healpix_geom_info(nside,1,NULL,&geom_info);

      /*      
      //Transform T
      time=0;
#ifdef _SPREC
      sharp_execute(SHARP_MAP2ALM,0,&(alms_2[0]),&(maps_2[0]),geom_info,
		    alm_info,1,0,&time,NULL);
#else //_SPREC
      sharp_execute(SHARP_MAP2ALM,0,&(alms_2[0]),&(maps_2[0]),geom_info,
		    alm_info,1,SHARP_DP,&time,NULL);
#endif //_SPREC
      */
      //Transform Q,U
#ifdef _SPREC
      sharp_execute(SHARP_MAP2ALM,2,&(alms_2[0]),&(maps_2[0]),geom_info,
		    alm_info,1,0,&time,NULL);
#else //_SPREC
      sharp_execute(SHARP_MAP2ALM,2,&(alms_2[0]),&(maps_2[0]),geom_info,
		    alm_info,1,SHARP_DP,&time,NULL);
#endif //_SPREC

      sharp_destroy_geom_info(geom_info);
      sharp_destroy_alm_info(alm_info);
    }
    else {
      he_map2alm(nside,lmax,nmaps_2,maps_2,alms_2);
    }
  }

  index_cl=0;
  for(i1=0;i1<nmaps_1;i1++) {
    int i2;
    fcomplex *alm1=alms_1[i1];
    for(i2=0;i2<nmaps_2;i2++) {
      int l;
      fcomplex *alm2=alms_2[i2];
      for(l=0;l<=lmax;l++) {
	int m;
	cls[index_cl][l]=creal(alm1[he_indexlm(l,0,lmax)])*creal(alm2[he_indexlm(l,0,lmax)]);

	for(m=1;m<=l;m++) {
	  long index_lm=he_indexlm(l,m,lmax);
	  cls[index_cl][l]+=2*(creal(alm1[index_lm])*creal(alm2[index_lm])+
			       cimag(alm1[index_lm])*cimag(alm2[index_lm]));
	}
	cls[index_cl][l]/=(2*l+1.);
      }
      index_cl++;
    }
  }

  for(i1=0;i1<nmaps_1;i1++)
    free(alms_1[i1]);
  free(alms_1);
  if(alms_1!=alms_2) {
    for(i1=0;i1<nmaps_2;i1++)
      free(alms_2[i1]);
    free(alms_2);
  }
}

void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname)
{
  fitsfile *fptr;
  int ii,status=0;
  char *ttype[]={"T","Q","U"};
  char *tform[]={"1E","1E","1E"};
  char *tunit[]={"mK","mK","mK"};
  float *map_dum=my_malloc(nside2npix(nside)*sizeof(float));

  if((nfields!=1)&&(nfields!=3)) {
    fprintf(stderr,"CRIME: nfields must be 1 or 3\n");
    exit(1);
  }

  fits_create_file(&fptr,fname,&status);
  fits_create_tbl(fptr,BINARY_TBL,0,nfields,ttype,tform,
		  tunit,"BINTABLE",&status);
  fits_write_key(fptr,TSTRING,"PIXTYPE","HEALPIX","HEALPIX Pixelisation",
		 &status);

  fits_write_key(fptr,TSTRING,"ORDERING","RING",
		 "Pixel ordering scheme, either RING or NESTED",&status);
  fits_write_key(fptr,TLONG,"NSIDE",&nside,
		 "Resolution parameter for HEALPIX",&status);
  fits_write_key(fptr,TSTRING,"COORDSYS","G",
		 "Pixelisation coordinate system",&status);
  fits_write_comment(fptr,
		     "G = Galactic, E = ecliptic, C = celestial = equatorial",
		     &status);
  for(ii=0;ii<nfields;ii++) {
    lint ip;
    for(ip=0;ip<nside2npix(nside);ip++)
      map_dum[ip]=(float)(tmap[ii][ip]);
    fits_write_col(fptr,TFLOAT,ii+1,1,1,nside2npix(nside),map_dum,&status);
  }
  fits_close_file(fptr, &status);
}

flouble *he_read_healpix_map(char *fname,long *nside,int nfield)
{
  //////
  // Reads a healpix map from file fname. The map will be
  // read from column #nfield. It also returns the map's nside.
  int status=0,hdutype,nfound,anynul;
  long naxes,*naxis,npix;
  fitsfile *fptr;
  flouble *map,nulval;
  char order_in_file[32];
  int nested_in_file=0;

  fits_open_file(&fptr,fname,READONLY,&status);
  fits_movabs_hdu(fptr,2,&hdutype,&status);
  fits_read_key_lng(fptr,"NAXIS",&naxes,NULL,&status);
  naxis=my_malloc(naxes*sizeof(long));
  fits_read_keys_lng(fptr,"NAXIS",1,naxes,naxis,&nfound,&status);
  fits_read_key_lng(fptr,"NSIDE",nside,NULL,&status);
  npix=12*(*nside)*(*nside);
  if(npix%naxis[1]!=0) {
    fprintf(stderr,"CRIME: WTFFF\n");
    exit(1);
  }

  if (fits_read_key(fptr, TSTRING, "ORDERING", order_in_file, NULL, &status)) {
    fprintf(stderr, "WARNING: Could not find %s keyword in in file %s\n",
            "ORDERING",fname);
    exit(1);
  }
  if(!strncmp(order_in_file,"NEST",4))
    nested_in_file=1;

  map=my_malloc(npix*sizeof(flouble));
#ifdef _SPREC
  fits_read_col(fptr,TFLOAT,nfield+1,1,1,npix,&nulval,map,&anynul,&status);
#else //_SPREC
  fits_read_col(fptr,TDOUBLE,nfield+1,1,1,npix,&nulval,map,&anynul,&status);
#endif //_SPREC
  free(naxis);

  fits_close_file(fptr,&status);

  flouble *map_ring;
  if(nested_in_file) {
    long ipring,ipnest;

    printf("read_healpix_map: input is nested. Transforming to ring.\n");
    map_ring=my_malloc(npix*sizeof(flouble));
    for(ipnest=0;ipnest<npix;ipnest++) {
      nest2ring(*nside,ipnest,&ipring);
      map_ring[ipring]=map[ipnest];
    }
    free(map);
  }
  else
    map_ring=map;

  return map_ring;
}

int he_ring_num(long nside,double z)
{
  //Returns ring index for normalized height z
  int iring;

  iring=(int)(nside*(2-1.5*z)+0.5);
  if(z>0.66666666) {
    iring=(int)(nside*sqrt(3*(1-z))+0.5);
    if(iring==0) iring=1;
  }

  if(z<-0.66666666) {
    iring=(int)(nside*sqrt(3*(1+z))+0.5);
    if(iring==0) iring=1;
    iring=4*nside-iring;
  }

  return iring;
}

static void get_ring_limits(long nside,int iz,long *ip_lo,long *ip_hi)
{
  long ir;
  long ipix1,ipix2;
  long npix=12*nside*nside;
  long ncap=2*nside*(nside-1);

  if((iz>=nside)&&(iz<=3*nside)) { //eqt
    ir=iz-nside+1;
    ipix1=ncap+4*nside*(ir-1);
    ipix2=ipix1+4*nside-1;
  }
  else {
    if(iz<nside) { //north
      ir=iz;
      ipix1=2*ir*(ir-1);
      ipix2=ipix1+4*ir-1;
    }
    else { //south
      ir=4*nside-iz;
      ipix1=npix-2*ir*(ir+1);
      ipix2=ipix1+4*ir-1;
    }
  }

  *ip_lo=ipix1;
  *ip_hi=ipix2;
}

long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip)
{
  long *pixlist;
  double z_hi=cos(theta1);
  double z_lo=cos(theta2);
  int irmin,irmax;

  if((theta2<=theta1)||
     (theta1<0)||(theta1>M_PI)||
     (theta2<0)||(theta2>M_PI)) {
    fprintf(stderr,"CRIME: wrong strip boundaries\n");
    exit(1);
  }

  irmin=he_ring_num(nside,z_hi);
  irmax=he_ring_num(nside,z_lo);

  //Count number of pixels in strip
  int iz;
  long npix_in_strip=0;
  for(iz=irmin;iz<=irmax;iz++) {
    long ipix1,ipix2;
    get_ring_limits(nside,iz,&ipix1,&ipix2);
    npix_in_strip+=ipix2-ipix1+1;
  }
  *npix_strip=npix_in_strip;
  pixlist=my_malloc(npix_in_strip*sizeof(long));

  //Count number of pixels in strip
  long i_list=0;
  for(iz=irmin;iz<=irmax;iz++) {
    long ipix1,ipix2,ip;
    get_ring_limits(nside,iz,&ipix1,&ipix2);
    for(ip=ipix1;ip<=ipix2;ip++) {
      pixlist[i_list]=ip;
      i_list++;
    }    
  }

  return pixlist;
}

void he_ring2nest_inplace(flouble *map_in,long nside)
{
  long npix=12*nside*nside;
  flouble *map_out=my_malloc(npix*sizeof(flouble));

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(map_in,nside,npix,map_out)
#endif //_HAVE_OMP
  {
    long ip;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<npix;ip++) {
      long inest;
      ring2nest(nside,ip,&inest);
      
      map_out[inest]=map_in[ip];
    } //end omp for
  } //end omp parallel
  memcpy(map_in,map_out,npix*sizeof(flouble));
  
  free(map_out);
}

void he_nest2ring_inplace(flouble *map_in,long nside)
{
  long npix=12*nside*nside;
  flouble *map_out=my_malloc(npix*sizeof(flouble));

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(map_in,nside,npix,map_out)
#endif //_HAVE_OMP
  {
    long ip;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<npix;ip++) {
      long iring;
      nest2ring(nside,ip,&iring);

      map_out[iring]=map_in[ip];
    } //end omp for
  } //end omp parallel
  memcpy(map_in,map_out,npix*sizeof(flouble));

  free(map_out);
}

void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest)
{
  long npix_in=nside2npix(nside_in);
  long npix_out=nside2npix(nside_out);

  if(nside_in>nside_out) {
    long ii;
    long np_ratio=npix_in/npix_out;
    double i_np_ratio=1./((double)np_ratio);
    
    for(ii=0;ii<npix_out;ii++) {
      int jj;
      double tot=0;

      if(nest) {
	for(jj=0;jj<np_ratio;jj++)
	  tot+=map_in[jj+ii*np_ratio];
	map_out[ii]=tot*i_np_ratio;
      }
      else {
	long inest_out;

	ring2nest(nside_out,ii,&inest_out);
	for(jj=0;jj<np_ratio;jj++) {
	  long iring_in;
	  
	  nest2ring(nside_in,jj+np_ratio*inest_out,&iring_in);
	  tot+=map_in[iring_in];
	}
	map_out[ii]=tot*i_np_ratio;
      }
    }
  }
  else {
    long ii;
    long np_ratio=npix_out/npix_in;
    
    for(ii=0;ii<npix_in;ii++) {
      int jj;
      
      if(nest) {
	flouble value=map_in[ii];

	for(jj=0;jj<np_ratio;jj++)
	  map_out[jj+ii*np_ratio]=value;
      }
      else {
	long inest_in;
	flouble value=map_in[ii];
	ring2nest(nside_in,ii,&inest_in);
	
	for(jj=0;jj<np_ratio;jj++) {
	  long iring_out;
	  
	  nest2ring(nside_out,jj+inest_in*np_ratio,&iring_out);
	  map_out[iring_out]=value;
	}
      }
    }
  }
}

//Transforms FWHM in arcmin to sigma_G in rad:
//         pi/(60*180*sqrt(8*log(2))
#define FWHM2SIGMA 0.00012352884853326381 
double *he_generate_beam_window(int lmax,double fwhm_amin)
{
  long l;
  double sigma=FWHM2SIGMA*fwhm_amin;
  double *beam=my_malloc((lmax+1)*sizeof(double));

  for(l=0;l<=lmax;l++)
    beam[l]=exp(-0.5*l*(l+1)*sigma*sigma);
  
  return beam;
}

void he_alter_alm(int lmax,double fwhm_amin,fcomplex *alms,double *window)
{
  double *beam;
  int mm;

  if(window==NULL) beam=he_generate_beam_window(lmax,fwhm_amin);
  else beam=window;

  for(mm=0;mm<=lmax;mm++) {
    int ll;
    for(ll=mm;ll<=lmax;ll++) {
      long index=he_indexlm(ll,mm,lmax);

      alms[index]*=beam[ll];
    }
  }

  if(window==NULL)
    free(beam);
}

/*
flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed)
{
  fcomplex *alms;
  int lmax_here=lmax;
  long npix=12*((long)nside)*nside;
  flouble *map=my_malloc(npix*sizeof(flouble));
  dam_rng_state *rng=dam_init_rng(seed);

  if(lmax>3*nside-1)
    lmax_here=3*nside-1;
  alms=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));

  int ll;
  for(ll=0;ll<=lmax_here;ll++) {
    int mm;
    double sigma=sqrt(0.5*cl[ll]);
    double r1,r2;
    r1=dam_randgauss(rng);
    alms[he_indexlm(ll,0,lmax_here)]=(fcomplex)(M_SQRT2*sigma*r1);

    for(mm=1;mm<=ll;mm++) {
      r1=dam_randgauss(rng);
      r2=dam_randgauss(rng);
      alms[he_indexlm(ll,mm,lmax_here)]=(fcomplex)(sigma*(r1+I*r2));
    }
  }
  he_alm2map(nside,lmax_here,1,&map,&alms);
  free(alms);
  dam_end_rng(rng);

  return map;
}
*/
