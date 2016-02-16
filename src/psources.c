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

#define NL_PSOURCES 256
#define LLOGMIN -5. //Only luminosities between 10^17 and
#define LLOGMAX 6.  //10^28 W/Hz will be considered

//////
// User defined functions and macros

#define FLUXMIN 1E-2 //Minimum flux in uJy

#define LX_TOY 2.11
#define ALPHA_TOY -0.633
#define RHO_TOY 3.2E-4
#define NU0_PSOURCES 1420. //Frequency for the luminosity function
static double l_z_function(double L,double z)
{
  //////
  // Returns the luminosity function for point sources
  // at redshift z with L given in units of 10^22 W/Hz
  //   dn(z,L)/dlog10(L)

  double fz;
  //  double rho_l=M_LN10*RHO_TOY*pow(L/LX_TOY,ALPHA_TOY+1)*exp(-L/LX_TOY);
  double rho_l=2.5*M_LN10*RHO_TOY*pow(L/LX_TOY,ALPHA_TOY+1)*exp(-L/LX_TOY); //Use this if RHO_TOY is per mag.

  if(z<1.5) fz=pow((1+z),3.1);
  else fz=17.124;

  return fz*rho_l;
}

#define SED_NORM 0.114461
static double spec_ed(double nu)
{
  //////
  // Point sources' SED
  // It should be normalized to 1 for the frequency
  // corresponding to the provided luminosity function.
  double nu_GHz=nu*0.001;

  return SED_NORM*(pow(nu_GHz,-0.1)+10*pow(nu_GHz,-0.75));
}

double bias_psources(double z)
{
  return 1;
}
// End of user definable functions
//////


double n_of_z_psources(ParamGetHI *par,double z)
{
  //////
  // Returns the number density of point sources
  // as a function of redshift integrated over luminosities
  int iz=(int)(z*par->glob_inv_dz);
  
  if((iz>=NZ_PSOURCES)||(iz<0)) return -1;
  else if(iz==NZ_PSOURCES-1) return par->nz_psources_arr[NZ_PSOURCES-1];
  else {
    double zi=iz*par->glob_dz;
    double nz=par->nz_psources_arr[iz]+
      (par->nz_psources_arr[iz+1]-par->nz_psources_arr[iz])*
      (z-zi)*par->glob_inv_dz;
    return nz;
  }
}

static double l_distribution(ParamGetHI *par,double L,double z)
{
  //////
  // Returns the pdf P(log10(L)|z)
  return l_z_function(L,z)/n_of_z_psources(par,z);
}

void setup_psources(ParamGetHI *par)
{
  //////
  // Calculates arrays and parameters for the luminosity
  // function of point sources
  int ii;
  double dlogL=(LLOGMAX-LLOGMIN)/NL_PSOURCES;
  par->glob_dz=par->z_max/NZ_PSOURCES;
  par->glob_inv_dz=1./par->glob_dz;

  for(ii=0;ii<NZ_PSOURCES;ii++) {
    int jj;
    double max_dist;
    double z=ii*par->glob_dz;
    
    par->nz_psources_arr[ii]=0;
    max_dist=-1;
    for(jj=0;jj<NL_PSOURCES;jj++) {
      double logL=LLOGMIN+(jj+0.5)*dlogL;
      double lfunc=l_z_function(pow(10.,logL),z);
      
      if(lfunc>=max_dist) max_dist=lfunc;
      par->nz_psources_arr[ii]+=dlogL*lfunc;
    }
    if(par->nz_psources_arr[ii]>0)
      par->max_Lpdf_arr[ii]=1.1*max_dist/par->nz_psources_arr[ii];
    else
      par->max_Lpdf_arr[ii]=0;
  }
}

double draw_luminosity(ParamGetHI *par,double z,gsl_rng *rng)
{
  //////
  // Generates a random luminosity in units of 10^22 W/Hz
  // following the luminosity function of point sources
  if(z>=par->z_max) return 0;
  else {
    double u,logL,L;
    int accept=0;
    int iz=(int)(z*par->glob_inv_dz);
    double maxpdf=par->max_Lpdf_arr[iz];

    while(!accept) {
      u=rng_01(rng);
      logL=LLOGMIN+(LLOGMAX-LLOGMIN)*rng_01(rng);
      L=pow(10.,logL);
      if(u*maxpdf<l_distribution(par,L,z))
	accept=1;
    }

    return L;
  }
}

#define FLUX2TEMP 3.2548291E-2 //c^2*(1 microJy)/(2*kB*(1 MHz)^2)/[1 mK]
#define LUM2FLUX 8.35774E7 //(10^22 Joules)/(4*pi*(1 Mpc)^2)/[1 microJy]
#define PREFAC_L2T
double temp_of_l(ParamGetHI *par,double L0,double nu_obs,
		 double z,double r,double dOmega)
{
  //  double s_nu=LUM2FLUX*L0*spec_ed((1+z)*nu_obs)*hhub*hhub/(r*r*(1+z)); //Use this is Lum is given in W/Hz
  double s_nu=LUM2FLUX*4*M_PI*L0*spec_ed((1+z)*nu_obs)*par->hhub*par->hhub/(r*r*(1+z)); //Use this it Lum is given in W/Hz/sr

  //  if(s_nu<FLUXMIN)
  //    return 0;
  //  else
  return FLUX2TEMP*s_nu/(dOmega*nu_obs*nu_obs);
}
