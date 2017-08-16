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
#include <assert.h>
#include "config.h"
#include "gumbel.h"
#include "mem_utils.h"

double
average(double *array, int n)
{
  double avg = 0.0;
  for (int i = 0;
       i < n;
       i++)
    avg += array[i];
  return avg/n;
};

double
stdev(double *array, double avg, int n)
{
  double std = 0.0;
  for (int i = 0;
       i < n;
       i++)
    std += (avg - array[i]) * (avg - array[i]);
  return sqrt(std/n);
};

int
index_of_max_element(double *x, int n_samples)
{
  double max = -DBL_MAX;
  int index_max = 0;
  for (int i = 0; i < n_samples; i++)
    if (x[i] > max)
      {
	max = x[i];
	index_max = i;
      }
  return index_max;
}

int
index_of_second_max_element(double *x, int n_samples, int max_index)
{
  double max = -DBL_MAX;
  int index_max = 0;

  for (int i = 0; i < n_samples; i++)
    if (x[i] > max && i != max_index)
      {
	max = x[i];
	index_max = i;
      }
  return index_max;
}


int
block_maxima(double *x,
	     int n_samples,
	     int samples_per_block)
{
  assert(samples_per_block < n_samples &&
	 (n_samples % samples_per_block == 0));

  double block[samples_per_block];
  memset(block, 0, samples_per_block*sizeof(double));
  int curr_pos = 0;
  for(int i = 0;
      i < n_samples;
      i++)
    {
      if (i % samples_per_block == 0)
	{
	  memcpy(block, x + i, samples_per_block*sizeof(double));
	  int max_index = index_of_max_element(block, samples_per_block);
	  x[curr_pos] = block[max_index];
	  ++curr_pos;
	}
    }
  memset(x+curr_pos, 0, (n_samples-curr_pos)*sizeof(double));
  return curr_pos; //number of block maxima
}

int
block_maxima_spectra(double **x,
		     int lowest_periodic_component_index,
		     int highest_periodic_component_index,
		     int n_samples,
		     int samples_per_block)
{
  int curr_pos = 0;
  for (int omega = lowest_periodic_component_index;
       omega < highest_periodic_component_index;
       omega++)
    curr_pos =
      block_maxima(x[omega], n_samples, samples_per_block);
  return curr_pos; //number of block maxima
}

void
gumbel_estimate_parametric(double *mu,
			   double *beta,
			   double *mean,
			   double *std,
			   double *random_maxima,
			   int n_samples)
{
  *mean = average(random_maxima, n_samples);
  *std = stdev(random_maxima, *mean, n_samples);

  *beta=(*std)*sqrt(6)/PI;
  *mu = (*mean) - EULER*(*beta);
}

double
gumbel_cdf(double x, double mu, double beta)
{
  return exp(-exp(-(x-mu)/beta));
}



