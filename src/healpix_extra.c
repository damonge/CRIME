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

void he_write_healpix_map(float **tmap,int nfields,long nside,char *fname)
{
  fitsfile *fptr;
  int ii,status=0;
  char *ttype[]={"T","Q","U"};
  char *tform[]={"1E","1E","1E"};
  char *tunit[]={"mK","mK","mK"};

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
    fits_write_col(fptr,TFLOAT,ii+1,1,1,nside2npix(nside),tmap[ii],&status);
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
  naxis=(long *)malloc(naxes*sizeof(long));
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

  map=(flouble *)my_malloc(npix*sizeof(flouble));
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
    map_ring=(flouble *)my_malloc(npix*sizeof(flouble));
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
  pixlist=(long *)my_malloc(npix_in_strip*sizeof(long));

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

void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest)
{
  long npix_in=nside2npix(nside_in);
  long npix_out=nside2npix(nside_out);

  if(nside_in==nside_out) {
    long ii;
    for(ii=0;ii<npix_out;ii++)
      map_out[ii]=map_in[ii];
  }
  else if(nside_in>nside_out) {
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
  double *beam=(double *)my_malloc((lmax+1)*sizeof(double));

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
