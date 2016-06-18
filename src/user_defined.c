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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double fraction_HI(double z)
{
  return 0.008*pow(1+z,0.6);
}

double bias_HI(double z)
{
  return 0.904+0.135*pow(1+z,1.696);
}
