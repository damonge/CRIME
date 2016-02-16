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

ParamsForGet *set_default_params_ForGet(void)
{
  ParamsForGet *pars=(ParamsForGet *)my_malloc(sizeof(ParamsForGet));

  sprintf(pars->fname_input,"default");
  sprintf(pars->fname_nutable,"default");
  sprintf(pars->fname_haslam,"default");
  sprintf(pars->fname_specin,"default");
  sprintf(pars->prefix_out,"default");
  pars->lmin=10;
  pars->lmax=1000;
  pars->nside=512;
  pars->seed=1234;
  pars->do_galaxy=0;
  pars->amp=57.0;
  pars->beta=1.1;
  pars->alpha=2.07;
  pars->xi=1.0;
  pars->nu_ref=130.0;
  pars->lref=1000;
  pars->do_polarization=0;
  pars->xi_polarization=3.0;
  pars->beta_polarization=2.5;
  sprintf(pars->fname_faraday,"default");

  return pars;
}

void free_params_ForGet(ParamsForGet *pars)
{
  int ii;
  for(ii=0;ii<3;ii++)
    free(pars->nutable[ii]);
  free(pars->nutable);
  free(pars);
}
