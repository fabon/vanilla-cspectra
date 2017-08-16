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
#include <math.h>
#include <float.h>
#include "mem_utils.h"

void*
allocate_and_zero(int nb_bytes)
{
  void *p = malloc(nb_bytes);
  memset(p, 0, nb_bytes);
  return p;
}

void*
fft_allocate_and_zero(int nb_bytes)
{
  void *p = fftw_malloc(nb_bytes);
  memset(p, 0, nb_bytes);
  return p;
}

double*
fft_allocate_and_nan(int nb_bytes)
{
  double *p = (double *)fftw_malloc(nb_bytes);
  double f = NAN;
  memset(p, f, nb_bytes);
  return p;
}

char **
allocate_and_cpy_string_array(char **src, int start, int n)
{
  char **p = (char **)malloc(n * sizeof(char*));
  for (int i = start; i < (n+start); i++)
    p[i-start] = strdup(src[i]);
  return p;
}

void*
allocate_and_cpy(void *src, int nb_bytes)
{
  void *p = malloc(nb_bytes);
  memcpy(p, src, nb_bytes);
  return p;
}

void*
fft_allocate_and_cpy(void *src, int nb_bytes)
{
  void *p = fftw_malloc(nb_bytes);
  memcpy(p, src, nb_bytes);
  return p;
}



