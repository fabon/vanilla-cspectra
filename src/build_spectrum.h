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
#ifndef BUILD_SPECTRUM_H
# define BUILD_SPECTRUM_H

#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <assert.h>

#include "complex.h"
#include "config.h"
#include "ioutils.h"

typedef struct t_spectral_coords
{
  int max_index;
  double energy;
  double period;
  int outreach;
  double pvalue;
  double bonferroni_pvalue;
  fftw_complex phase;
  double bandedEnergy;
  double spectrumEnergy;
} t_spectral_coords;

int
padd0s(double *x, int n_samples, int nfft);

int
save_series(double *x, char *signal_name, char *buff, char *out_path, int out_path_len, int n_samples);

int
save_spectrum(double *spectrum, char *signal_name, char *buff, char *out_path, int out_path_len, int real_size);

int
save_periodic_coord(FILE *f,
		    const char *signal_name,
		    t_spectral_coords *spectral_coord);
int
save_results(FILE *f,
	     char **signal_names,
	     int nb_signals_to_save,
	     t_spectral_coords *results);

int
save_vector(char *path,
	    double *data,
	    int n);

int
save_matrix(char *path,
	    double **data,
	    int n, int m);

/* int */
/* power_spectral_density(fftw_complex *X, double *psd, int real_size, int nfft); */

int
permute_data(double *x,
             int n_samples,
	     unsigned int *seed);

int
permute_data_with_vol(double *x,
		      double *vol,
		      int n_samples,
		      unsigned int *seed);

double
linear_detrend(double *x, int n_samples);

double
cste_atan_standardized_nonnan(double *x, int n_samples);

void
smooth_signal(double *x, int n_samples, int window_len);

void
interpolate_signal(double *y, int n_samples);

int
window_hann(double *x, int n_samples);

double
my_variance(double *x, int n_samples);

double
unit_energy(double *x, int n_samples);

void
high_pass_filter(double *x, int n_samples, double *x_buffer);

double
center_signal(double *x, int n_samples);

void
preprocess(double *x, int n_samples, double *x_buffer,
	   t_dump_info *dump_info);


double
extract_largest_frequency(double *spectrum,
			  fftw_complex *phases,
			  int lowest_periodic_component_index,
			  int highest_periodic_component_index,
			  t_spectral_coords *max_spectral_coords);

double
extract_nth_largest_frequency(int* indexes_largest,
			      double *spectrum,
			      fftw_complex *phases,
			      int lowest_periodic_component_index,
			      int highest_periodic_component_index,
			      t_spectral_coords *max_spectral_coords);


void
extract_time_process(double *x,
		     double *pattern,
		     double *median_buff,
		     int n_samples,
		     int pattern_period);

void
build_estimates(double *population_spectrum,
		double *population_variances,
		int lowest_periodic_component_index,
		int highest_periodic_component_index,
		int population_size);

int
build_fft(double *spectrum,
	  fftw_complex *phases,
	  double *x,
	  fftw_complex *X,
	  fftw_plan  *plan,
	  int lowest_periodic_component_index,
	  int highest_periodic_component_index,
	  double *bandedEnergy,
	  double *spectrumEnergy);

void
estimate_relative_frequency(double *x,
			    double *vol,
			    int n_samples);

double to2pi(double x);

double topi(double x);

double wrapto2pi(double x);


#endif /* !-- BUILD_SPECTRUM_H --! */



