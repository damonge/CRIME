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
#include "common_jt.h"

ParamsJoinT *set_default_params_JoinT(void)
{
  ParamsJoinT *pars=(ParamsJoinT *)my_malloc(sizeof(ParamsJoinT));

  sprintf(pars->fname_input,"default");
  sprintf(pars->prefix_cosmo,"default");
  sprintf(pars->prefix_synchrotron,"default");
  sprintf(pars->prefix_galactic_freefree,"default");
  sprintf(pars->prefix_extragalactic_freefree,"default");
  sprintf(pars->prefix_point_sources,"default");
  sprintf(pars->prefix_custom,"custom");
  sprintf(pars->prefix_out,"default");
  sprintf(pars->fname_nutable,"default");
  sprintf(pars->fname_mask,"default");
  pars->n_nu=-1;
  pars->nside=-1;
  pars->do_polarization=0;
  pars->polarization_leakage=-1.0;
  pars->dish_diameter=-1;
  pars->n_dish=-1;
  pars->tsys=-1;
  pars->time_tot=-1;
  pars->mask=NULL;
  pars->fsky=-1;
  return pars;
}

void free_params_JoinT(ParamsJoinT *pars)
{
  int ii;
  for(ii=0;ii<3;ii++)
    free(pars->nutable[ii]);
  free(pars->nutable);
  free(pars->mask);
  free(pars);
}
