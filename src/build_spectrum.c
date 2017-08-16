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
#include "build_spectrum.h"
#include "read_csv.h"

int
padd0s(double *x, int n_samples, int nfft)
{
  for (int i = n_samples; i < nfft; i++)
    x[i] = 0.0;
  return 1;
}

int
permute_data(double *x,
             int n_samples,
	     unsigned int *seed)
{
  for (int k = n_samples-1;
       k>0;
       k--)
    {
      int j = rand_r(seed) % (k+1);
      double temp = x[j];
      x[j] = x[k];
      x[k] = temp;
    }
  return 1;
};

int
permute_data_with_vol(double *x,
		      double *vol,
		      int n_samples,
		      unsigned int *seed)
{
  double order[n_samples];
  memset(order, 0, n_samples*sizeof(int));

  int n_non_zero_samples = 0;
  int min_vol = 0;
  for (int k = 0;
       k < n_samples;
       k++)
    if (vol[k] > min_vol)
      {
	order[n_non_zero_samples] = k;
	++n_non_zero_samples;
      }

  permute_data(order, n_non_zero_samples, seed);

  int pos = 0;
  for (int k = 0;
       k < n_samples;
       k++)
    {
      if (vol[k] > min_vol)
	{
	  double temp = x[k];
	  x[k] = x[(int)order[pos]];
	  x[(int)order[pos]]=temp;
	  ++pos;
	}
      else
      	x[k] = min_vol;
    }
  return 1;
};

static double
max(double *array, int n)
{
  double res = 0.0;
  for (int i = 0;
       i < n;
       i++)
    res = array[i] > res ? array[i] : res;
  return res;
};

double
create_permutation_threshold(double *max_random_components_energy)
{
  return max(max_random_components_energy, NPERMS);
}

int
save_series(double *x,
	    char *signal_name,
	    char *buff,
	    char *out_path,
	    int out_path_len,
	    int n_samples)
{
  sprintf(out_path+out_path_len, "%s_%i_series.csv",
  	  signal_name,
  	  NFFT);
  printf("saving [%s]\n", out_path);

  FILE *f = fopen(out_path, "w");
  if (!f)
    return 0;

  memset(buff, 0, NFFT*SUBBUFFSIZE*sizeof(char));
  int len = 0;
  for (int i = 0;
       i < n_samples;
       i++)
    {
      len += sprintf(buff + len, "%e\n", x[i]);
    }
  fwrite(buff, len * sizeof(char), 1, f);
  fclose(f);
  return 1;
}

int
save_spectrum(double *spectrum,
	      char *signal_name,
	      char *buff,
	      char *out_path,
	      int out_path_len,
	      int real_size)
{
  memset(out_path+out_path_len, 0, SUBBUFFSIZE-out_path_len);
  memcpy(out_path+out_path_len, signal_name, strlen(signal_name)*sizeof(char));

  strcat(out_path, "_");
  char b[SUBBUFFSIZE];
  sprintf(b, "%i", NFFT);
  strcat(out_path, b);
  strcat(out_path, "_spectrum.csv");
  FILE *f = fopen(out_path, "w");
  if (!f)
    return 0;

  memset(buff, 0, NFFT*SUBBUFFSIZE*sizeof(char));

  int len = 0;
  for (int i = 0;
       i < real_size;
       i++)
    len += sprintf(buff + len, "%i,%e\n",
		   i,
		   spectrum[i]);
  fwrite(buff, len * sizeof(char), 1, f);
  fclose(f);
  return 1;
}


int
save_vector(char *path,
	    double *data,
	    int n)
{
  FILE *f = fopen(path, "w");
  if (!f)
    return 0;

  for (int i = 0;
       i < n;
       i++)
    if (i < n - 1)
      fprintf(f, "%e,", data[i]);
    else
      fprintf(f, "%e\n", data[i]);
  fflush(f);
  fclose(f);
  return 1;
}

int
save_matrix(char *path,
	    double **data,
	    int n, int m)
{
  FILE *f = fopen(path, "w");
  if (!f)
    return 0;

  for (int i = 0;
       i < n;
       i++)
    {
      for (int j = 0;
	   j < m;
	   j++)
	{
	  if (j < m - 1)
	    fprintf(f, "%e,", data[i][j]);
	  else
	    fprintf(f, "%e\n", data[i][j]);
	}
    }
  fflush(f);
  fclose(f);
  return 1;
}

static void
convert_phase_to_time(char *date_str,
		      const char *init_date,
		      double phase,
		      double period)
{
  struct tm date;
  memset(&date, '\0', sizeof(date));
  setenv("TZ", "", 1);
  tzset();
  if (init_date != 0)
    strptime(init_date, DATE_FORMAT, &date);
  else
    strptime("01-01-2015 00:00:00", DATE_FORMAT, &date);

  double time = (1-wrapto2pi(phase)/(2*PI))*period;
  time = floor(time);
  add_hour(&date,time);
  strftime(date_str, SUBBUFFSIZE, DATE_FORMAT, &date);
}


int
save_results(FILE *f,
	     char **signal_names,
	     int nb_signals_to_save,
	     t_spectral_coords *results)
{
  char *buffer = (char*)malloc(_1MB*10 * sizeof(char));
  setbuf(f, buffer);

  char date_str[SUBBUFFSIZE];
  memset(date_str,0,SUBBUFFSIZE*sizeof(char));
  for (int j = 0; j < nb_signals_to_save; j++)
    {
      if (results[j].period == 365 || results[j].period == 182.5)
	{

	  convert_phase_to_time(date_str,
				"01-01-2010 00:00:00",
				results[j].phase[0],
				results[j].period);
	  fprintf(f, "\"%s\",%e,%e,%e,%i,%e,%e,%e,%e,%s\n",
		  signal_names[j],
		  results[j].energy,
		  results[j].period,
		  results[j].phase[0],
		  results[j].outreach,
		  results[j].pvalue,
		  results[j].bonferroni_pvalue,
		  results[j].bandedEnergy,
		  results[j].spectrumEnergy,
		  date_str);
	}
      else
	fprintf(f, "\"%s\",%e,%e,%e,%i,%e,%e,%e,%e\n",
		signal_names[j],
		results[j].energy,
		results[j].period,
		results[j].phase[0],
		results[j].outreach,
		results[j].pvalue,
		results[j].bonferroni_pvalue,
		results[j].bandedEnergy,
		results[j].spectrumEnergy);
    }
  fflush(f);
  free(buffer);
  return 0;
}

static double
linear_interpolate(double x_left, double y_left,
		   double x_right, double y_right,
		   double x)
{
  double a = (y_right - y_left)/(x_right - x_left);
  double b = y_right - a*x_right;
  return a*x+b;
}

static int
find_nonnan_left(double *y, int start)
{
  int i = start;
  while (isnan(y[i]) && i >= 0)
    --i;
  return i;
}

static int
find_nonnan_right(double *y, int start, int n_samples)
{
  int i = start;
  while (isnan(y[i]) && i <= n_samples - 1)
    ++i;
  if (i == n_samples)
    return -1;
  return i;
}

void
interpolate_signal(double *y, int n_samples)
{
  int x_left=0;
  int x_right=0;
  for (int x = 1; x < n_samples - 1; x++)
    if (isnan(y[x]))
      {
	x_left = find_nonnan_left(y, x);
	x_right = find_nonnan_right(y, x, n_samples);
	if (x_left == -1 || x_right == -1)
	  return;
	while (isnan(y[x]))
	  {
	    y[x] = linear_interpolate(x_left, y[x_left],
				      x_right, y[x_right],
				      x);
	    x++;
	  }
      }
}


double
linear_detrend(double *x, int n_samples)
{
  double sx = 0.0;
  double sy = 0.0;
  for (int i = 0; i < n_samples; i++)
    {
      sy += x[i];
      sx += i;
    }

  sy /= n_samples;
  sx /= n_samples;

  double syy = 0.0;
  double sxx = 0.0;
  double sxy = 0.0;
  for (int i = 0; i < n_samples; i++)
    {
      syy += (x[i]-sy)*(x[i]-sy);
      sxx += (i-sx)*(i-sx);
      sxy += (i-sx)*(x[i]-sy);
    }

  double slope = sxy/sxx;
  double intercept = sy - slope*sx;
  for (int i = 0; i < n_samples; i++)
    {
      x[i] -=  (i*slope + intercept);
    }

  return slope;
}

double
cste_atan_standardized_nonnan(double *x, int n_samples)
{
  double nonnan_mean = 0.0;
  double nonnan_stdev = EPS_NORMALISATION;
  int n_nonnan = 0;
  for (int i = 0;
       i < n_samples;
       i++)
    if (!isnan(x[i]))
      {
	nonnan_mean += x[i];
	++n_nonnan;
      }
  nonnan_mean /= n_nonnan;
  for (int i = 0;
       i < n_samples;
       i++)
    if (!isnan(x[i]))
      nonnan_stdev += (x[i]-nonnan_mean)*(x[i]-nonnan_mean);
  nonnan_stdev /= n_nonnan;
  nonnan_stdev = sqrt(nonnan_stdev);
  double zero_mean = 0.0;
  for (int i = 0;
       i < n_samples;
       i++)
    if (!isnan(x[i])) {
      x[i] = (x[i] - nonnan_mean)/(nonnan_stdev);
      zero_mean += x[i];
      x[i] = atan((1.0/SCALE_ATAN)*x[i]);
    }
  return zero_mean/n_nonnan;
}

int
window_hann(double *x, int n_samples)
{
  for (int i = 0; i < n_samples; i++)
    x[i] = x[i] * 0.5 * (1-cos(2 * M_PI * i * (1.0/(n_samples-1))));
  return 1;
}

double
my_variance(double *x, int n_samples)
{
  double mean = 0.0;
  for (int i = 0; i < n_samples; i++)
    mean += x[i];
  mean /= n_samples;

  double energy = 0.0;
  for (int i = 0; i < n_samples; i++)
    energy += (x[i]-mean)*(x[i]-mean);

  energy/=n_samples;
  energy = sqrt(energy);
  return energy;

}

double
unit_energy(double *x, int n_samples)
{
  double energy = 0.0;

  for (int i = 0; i < n_samples; i++)
    energy += x[i]*x[i];

  energy/=n_samples;
  energy = sqrt(energy);
  for (int i = 0; i < n_samples; i++)
    x[i] = x[i]/energy;
  return energy;
}

double
extract_largest_frequency(double *spectrum,
			  fftw_complex *phases,
			  int lowest_periodic_component_index,
			  int highest_periodic_component_index,
			  t_spectral_coords *max_spectral_coords)
{
  double E_max = -DBL_MAX;
  int E_max_index = 0;
  for (int omega = lowest_periodic_component_index;
       omega < highest_periodic_component_index;
       omega++)
    {
      if (spectrum[omega] > E_max)
	{
	  E_max = spectrum[omega];
	  E_max_index = omega;
	}
    }
  
  if (!max_spectral_coords)
    return E_max;
  max_spectral_coords->max_index = E_max_index;
  max_spectral_coords->energy = E_max;
  max_spectral_coords->period = ((1.0/(double)E_max_index)*NFFT);
  max_spectral_coords->phase[0]=phases != 0 ? phases[E_max_index][0] : NAN;
  max_spectral_coords->phase[1]=NAN;
  return E_max;
}

static int
find_in_array(int* indexes_largest,
	      int index)
{
  if (!indexes_largest)
    return 0;
  int* p = indexes_largest;
  while (p && *p != -1 &&
	 *p != index)
    ++p;
  return !(*p == -1);
}

double
extract_nth_largest_frequency(int* indexes_largest,
			      double *spectrum,
			      fftw_complex *phases,
			      int lowest_periodic_component_index,
			      int highest_periodic_component_index,
			      t_spectral_coords *max_spectral_coords)
{
  double E_max = -DBL_MAX;
  int E_max_index = 0;
  for (int omega = lowest_periodic_component_index;
       omega < highest_periodic_component_index;
       omega++)
    {
      if (!find_in_array(indexes_largest, omega) &&
	  spectrum[omega] > E_max)
	{
	  E_max = spectrum[omega];
	  E_max_index = omega;
	}
    }
  if (!max_spectral_coords)
    return E_max;
  max_spectral_coords->max_index = E_max_index;
  max_spectral_coords->energy = E_max;
  max_spectral_coords->period = ((1.0/(double)E_max_index)*NFFT);
  max_spectral_coords->phase[0]=phases != 0 ? phases[E_max_index][0] : NAN;
  max_spectral_coords->phase[1]=NAN;
  return E_max;
}

void
build_estimates(double *population_spectrum,
		double *population_variances,
		int lowest_periodic_component_index,
		int highest_periodic_component_index,
		int population_size)
{
    for (int omega = lowest_periodic_component_index;
    	 omega < highest_periodic_component_index;
    	 omega++)
      {
	population_spectrum[omega] /= population_size;
	population_variances[omega] /= population_size;
	population_variances[omega] -=
	  population_spectrum[omega]*population_spectrum[omega];

	population_variances[omega] = sqrt(population_variances[omega]);
      }
}

int
build_fft(double *spectrum,
	  fftw_complex *phases,
	  double *x,
	  fftw_complex *X,
	  fftw_plan  *plan,
	  int lowest_periodic_component_index,
	  int highest_periodic_component_index,
	  double *bandedEnergy,
	  double *spectrumEnergy)
{
  double periodogram[NFFT];
  memset(periodogram, 0, NFFT*sizeof(double));

  fftw_execute_dft_r2c(*plan, x, X);

  if (spectrumEnergy)
    {
      *spectrumEnergy = 0.0;
      for (int omega = 0;
	   omega < NFFT;
	   omega++)
	{
	  *spectrumEnergy +=
	    (2.0/((double)NFFT*NFFT))*(X[omega][0]*X[omega][0] + X[omega][1]*X[omega][1]);
	}
    }


  double Ebandlimited = 0.0;
  for (int omega = lowest_periodic_component_index;
       omega < highest_periodic_component_index;
       omega++)
    {
      periodogram[omega] =
	(2.0/((double)NFFT*NFFT))*(X[omega][0]*X[omega][0] + X[omega][1]*X[omega][1]);
      Ebandlimited += periodogram[omega];
    }
  for (int omega = lowest_periodic_component_index;
       omega < highest_periodic_component_index;
       omega++)
    {
      spectrum[omega] = periodogram[omega];
      if (phases)
	{
	  phases[omega][0] = atan2(X[omega][1], X[omega][0]);
	  phases[omega][1] = NAN;
	}
    }

  if (bandedEnergy)
    *bandedEnergy = Ebandlimited/(*spectrumEnergy);
  return 1;
}

void
estimate_relative_frequency(double *x,
			    double *vol,
			    int n_samples)
{
  if (PREPOC_LVL < RF)
    return;

  for (int j=0;
       j<n_samples;
       j++)
    {
      x[j] = vol[j] != 0 ? x[j]/vol[j] : NAN;
    }
}

void
smooth_signal(double *x, int n_samples, int window_len)
{
  assert(window_len % 2 == 1);
  double smoothed_x[n_samples];
  memset(smoothed_x, 0, n_samples*sizeof(double));
  double window[window_len];
  memset(window, 0, window_len*sizeof(double));

  int span = (window_len-1)/2;
  int left_span = 0;
  int right_span = 0;
  double sum_left = 0.0;
  double sum_right = 0.0;
  for (int i=1; i < span; i++)
    sum_right += x[i];
  for (int i = 0; i < n_samples; i++)
    {
      left_span = (i - span) > 0 ? span : i;
      right_span = (i + span) < n_samples ? span : (n_samples - i - 1);

      smoothed_x[i]=
	(sum_left + sum_right + x[i])/(left_span + right_span + 1);
      if (span > 0)
	sum_left += x[i];
      if (left_span == span)
	sum_left -= x[i - span];

      if (span > 0)
	sum_right -= x[i+1];
      if (i+span < n_samples)
	sum_right += x[i+span];
    }

  memcpy(x, smoothed_x, n_samples*sizeof(double));
}

void
high_pass_filter(double *x, int n_samples, double *x_buffer)
{
  memcpy(x_buffer, x, n_samples*sizeof(double));
  smooth_signal(x_buffer, n_samples, HIGH_PASS_WINDOW);
  for (int i = 0; i < n_samples; i++)
    x[i] -= x_buffer[i];
}

double
center_signal(double *x, int n_samples)
{
  double mean = 0.0;
  for (int i = 0; i < n_samples; i++)
    mean += x[i];
  mean /= n_samples;
  for (int i = 0; i < n_samples; i++)
    x[i] -= mean;
  return mean;
}

void
preprocess(double *x, int n_samples, double *x_buffer,
	   t_dump_info *dump_info)
{
  char name_buffer[SUBBUFFSIZE];
  if (dump_info)
    {
      sprintf(name_buffer, "%s_RF.csv",
  	      dump_info->dump_name);
      save_series(x,
  		  name_buffer,
  		  dump_info->output_buffer,
  		  dump_info->output_path,
  		  dump_info->export_dir_len,
  		  n_samples);
    }

  high_pass_filter(x, n_samples, x_buffer);


  center_signal(x, n_samples);

  unit_energy(x, n_samples);
  if (dump_info)
    {
      sprintf(name_buffer, "%s_INPUT.csv",
  	      dump_info->dump_name);
      save_series(x,
  		  name_buffer,
  		  dump_info->output_buffer,
  		  dump_info->output_path,
  		  dump_info->export_dir_len,
  		  n_samples);
    }
}

double
wrapto2pi(double x)
{
  return ((x > 0 && x < 2*PI) ? x : fmod(x,2*PI));
}

double
to2pi(double x)
{
  return (x > 0 ? x : (2*PI + x));
}

double
topi(double x)
{
  return (x <= PI ? x : (2*PI - x));
}



