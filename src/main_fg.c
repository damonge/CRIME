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
#include "common_fg.h"

static void run_foregrounds(ParamsForGet *pars)
{
  int ii;
  flouble **maps_q,**maps_u;
  flouble **maps_i=get_sck_maps(pars);

  if((pars->do_galaxy)&&(pars->do_polarization)) {
    do_polarization(pars,&maps_q,&maps_u);
  }
  else {
    maps_q=NULL;
    maps_u=NULL;
  }

  write_maps(pars->prefix_out,pars->nside,pars->n_nu,
	     maps_i,maps_q,maps_u);
  
  for(ii=0;ii<pars->n_nu;ii++)
    free(maps_i[ii]);
  free(maps_i);

  if((pars->do_galaxy)&&(pars->do_polarization)) {
    for(ii=0;ii<pars->n_nu;ii++) {
      free(maps_q[ii]);
      free(maps_u[ii]);
    }
    free(maps_q);
    free(maps_u);
  }
}

int main(int argc,char **argv)
{
  char fname_input[256];
  ParamsForGet *pars;

  if(argc!=2) {
    fprintf(stderr,"Usage: ./ForGet file_name\n");
    exit(1);
  }
  sprintf(fname_input,"%s",argv[1]);

  setbuf(stdout,NULL);
  printf("\n");
  printf("|-------------------------------------------------|\n");
  printf("|                      ForGet                     |\n");
  printf("|-------------------------------------------------|\n\n");
  
  timer(4);

  pars=read_input_params_ForGet(fname_input);

  run_foregrounds(pars);

  free_params_ForGet(pars);
  timer(5);

  printf("\n");
  printf("|-------------------------------------------------|\n\n");

  return 0;
}
