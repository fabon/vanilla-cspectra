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
#ifndef GUMBEL_H
# define GUMBEL_H

double
average(double *array, int n);

double
stdev(double *array, double avg, int n);

int
index_of_max_element(double *x, int n_samples);

int
index_of_second_max_element(double *x, int n_samples, int max_index);

int
block_maxima(double *x,
	     int n_samples,
	     int samples_per_block);

int
block_maxima_spectra(double **x,
		     int lowest_periodic_component_index,
		     int highest_periodic_component_index,
		     int n_samples,
		     int samples_per_block);

void
gumbel_estimate_parametric(double *mu,
			   double *beta,
			   double *mean,
			   double *std,
			   double *random_maxima,
			   int n_samples);
double
gumbel_cdf(double x, double mu, double beta);

#endif /* ! GUMBEL_H */



