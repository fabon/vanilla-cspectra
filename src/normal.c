/**
 ** Copyright 2017 Fabon Dzogang
 **
 ** This file is part of vanilla-cspectra.
 ** vanilla-cspectra is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** any later version.
 ** vanilla-cspectra is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU General Public License for more details.
 ** You should have received a copy of the GNU General Public License
 ** along with vanilla-cspectra.  If not, see <http://www.gnu.org/licenses/>.
 **/
#include "normal.h"

double
random_uniform(unsigned int *seed)   /* uniform distribution, (0..1] */
{
  return (rand_r(seed)+1.0)/(RAND_MAX+1.0);
}

double
random_normal(unsigned int *seed)  /* normal distribution, centered on 0, std dev 1 */
{
  double t1=random_uniform(seed);
  double t2=random_uniform(seed);
  return sqrt(-2*log(t1)) * cos(2*M_PI*t2);
}



