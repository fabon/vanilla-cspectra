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
#include <ctype.h>
#include "run_thread.h"
#include "ioutils.h"
#include "read_csv.h"
#include "mem_utils.h"
#include "normal.h"

static char *
extract_data_from_split_files(char *data_buff_str,
                              int file_index,
                              char **paths)
{
  int content_size = 0;

  if (!(content_size = append_file_content_to_buff(data_buff_str,
						   content_size,
						   paths[file_index])))
      return 0;
  data_buff_str[content_size] = '\0';

  return data_buff_str;
}

static void
launch(t_thread_data *thread_data)
{
  char *output_buffer = (char*)malloc(NFFT*SUBBUFFSIZE*sizeof(char));

  double *time_series_data = (double*) fft_allocate_and_zero(sizeof(double) * MAX_LEN);
  for (int l = 0; l < MAX_LEN ; l++)
    time_series_data[l]=NAN;

  double *x = (double*) fft_allocate_and_zero(sizeof(double) * MAX_LEN);
  double *x_buffer = (double*) fft_allocate_and_zero(sizeof(double) * MAX_LEN);
  fftw_complex *X = (fftw_complex*) fft_allocate_and_zero(sizeof(fftw_complex) * NFFT);
  double *signal_spectrum = (double*) fft_allocate_and_zero(sizeof(double) * NFFT);
  fftw_complex *signal_phases = (fftw_complex*) fft_allocate_and_zero(sizeof(fftw_complex) * NFFT);
  double *noise_spectrum = (double*) fft_allocate_and_zero(sizeof(double) * NFFT);

  double *largest_population_samples = (double*)fft_allocate_and_zero(sizeof(double) * NPERMS);

  double *vol = thread_data->vol;
  char *data_buff_str = (char*)malloc(100*_1MB*sizeof(char));
  char buff[SUBBUFFSIZE];
  t_spectral_coords *results =
    (t_spectral_coords *)malloc(thread_data->nb_loaded_signals * sizeof (t_spectral_coords));

  t_spectral_coords **results_remaining =
    (t_spectral_coords **)malloc(10 * sizeof(t_spectral_coords*));
  for (int i = 0; i < 10;i++)
    results_remaining[i] =
      (t_spectral_coords *)malloc(thread_data->nb_loaded_signals * sizeof (t_spectral_coords));

  int file_index = 0;
  for (file_index = 0;
       file_index < thread_data->nb_loaded_signals;
       file_index++)
    {
      results[file_index].max_index = -1;
      results[file_index].energy = -1;
      results[file_index].period = -1;
      results[file_index].outreach = 0;
      results[file_index].pvalue = -1;
      results[file_index].phase[0] = -1;
      results[file_index].phase[1] = -1;
    }
  for (file_index = 0;
       file_index < thread_data->nb_loaded_signals;
       file_index++)
    {
      fprintf(stdout,
	      "thread-%i: %i/%i\n",
	      thread_data->thread_id,
	      file_index,
	      thread_data->nb_loaded_signals);

      extract_data_from_split_files(data_buff_str,
      				    file_index,
      				    thread_data->paths);
      memset(time_series_data, 0, MAX_LEN*sizeof(double));
      int n_samples =0;
      if (n_samples < MAX_LEN)
      	{
      	  n_samples = 0;
      	  memset(time_series_data, 0, MAX_LEN*sizeof(double));
      	  n_samples = read_data_from_buff(time_series_data,
      					  data_buff_str,
      					  &thread_data->from,
      					  &thread_data->to);
      	}

      printf("n samples %s [%i]\n", thread_data->signal_names[file_index], n_samples);
      if (n_samples == -1)
	continue;
      memcpy(x, time_series_data, MAX_LEN*sizeof(double));

      for (int i=0; i<n_samples;i++)
      	x[i] = (thread_data->missing[i]) ? NAN : x[i];
      interpolate_signal(x, n_samples);
      smooth_signal(x, n_samples, SMOOTHING_WINDOW);


      estimate_relative_frequency(x,
      				  vol,
      				  n_samples);



      memset(signal_spectrum, 0, NFFT*sizeof(double));
      memset(signal_phases, 0, NFFT * sizeof(fftw_complex));
      t_dump_info dump_info;
      dump_info.dump_name = thread_data->signal_names[file_index];
      dump_info.output_buffer = output_buffer;
      dump_info.output_path = thread_data->output_path;
      dump_info.export_dir_len = thread_data->export_dir_len;

      preprocess(x, n_samples, x_buffer, &dump_info);
      memcpy(time_series_data, x, MAX_LEN*sizeof(double));

      build_fft(signal_spectrum,
		signal_phases,
		x,
		X,
		thread_data->plan,
		thread_data->lowest_periodic_component_index,
		thread_data->highest_periodic_component_index,
		&(results[file_index].bandedEnergy),
		&(results[file_index].spectrumEnergy));

      extract_largest_frequency(signal_spectrum,
      				signal_phases,
      				thread_data->lowest_periodic_component_index,
      				thread_data->highest_periodic_component_index,
      				&(results[file_index]));
      for (int perm_index = 0;
           perm_index < NPERMS;
           perm_index++)
        {
      	  memcpy(x, time_series_data, MAX_LEN*sizeof(double));
      	  permute_data(x,
	  	       n_samples,
	  	       &(thread_data->seed));

	  if (perm_index%1000 ==0)
	    printf("%s %i permutations\n",
		   thread_data->signal_names[file_index],
		   perm_index);

	  build_fft(noise_spectrum,
		    0,
		    x,
		    X,
		    thread_data->plan,
		    thread_data->lowest_periodic_component_index,
		    thread_data->highest_periodic_component_index,
		    0,0);
	  /**
	   ** largest across frequencies
	   **/

	  t_spectral_coords dumny;
	  extract_largest_frequency(noise_spectrum,
				    0,
				    thread_data->lowest_periodic_component_index,
				    thread_data->highest_periodic_component_index,
				    &dumny);
	  double permutation_score = dumny.energy;

	  largest_population_samples[perm_index] =
	    permutation_score;

	  results[file_index].outreach +=
	    results[file_index].energy <= permutation_score;
	}
      results[file_index].pvalue =
	((double)results[file_index].outreach) / ((double)NPERMS);
      results[file_index].bonferroni_pvalue = results[file_index].pvalue*NB_CMPS;
    }

  sprintf(buff, "%s/%s_%i.csv",
  	  thread_data->export_dir,
  	  PERIODIC_COORDS_EXPORT_PATH,
  	  thread_data->thread_id);
  thread_data->out = fopen(buff, "w");

  save_results(thread_data->out,
  	       thread_data->signal_names,
	       file_index,
  	       results);

  fclose(thread_data->out);

  for (int i = 0 ;i < 10;i++)
    {
      sprintf(buff, "%s/%s_remain_%i_%i.csv",
	      thread_data->export_dir,
	      PERIODIC_COORDS_EXPORT_PATH,
	      i,
	      thread_data->thread_id);
      thread_data->out = fopen(buff, "w");

      save_results(thread_data->out,
		   thread_data->signal_names,
		   file_index,
		   results_remaining[i]);

      fclose(thread_data->out);
    }
  if (largest_population_samples)
    free(largest_population_samples);
  if (x)
    fftw_free(x);
  if (x_buffer)
    free(x_buffer);
  if (X)
    fftw_free(X);
  if (signal_spectrum)
    fftw_free(signal_spectrum);
  if (signal_phases)
    fftw_free(signal_phases);
  if (noise_spectrum)
    free(noise_spectrum);
  if (output_buffer)
    free(output_buffer);
}

void
clean_data(t_thread_data *data)
{
  if (data->vol)
    free(data->vol);
  if (data->missing)
    free(data->missing);
  if (data->export_dir)
    free(data->export_dir);
  if (data->signal_names)
    {
      for (int i = 0; i < data->nb_loaded_signals; i++)
	if (data->signal_names[i])
	  free(data->signal_names[i]);
      free(data->signal_names);
    }

  if (data->paths)
    {
      for (int i = 0; i < data->nb_loaded_signals; i++)
	if (data->paths[i])
	  free(data->paths[i]);
      free(data->paths);
    }
}

void* run_thread(void *arg_data)
{
  t_thread_data *thread_data = (t_thread_data *)arg_data;

  char buff[SUBBUFFSIZE];
  sprintf(buff, "%s/%s_%i.csv",
  	  thread_data->export_dir,
  	  PERIODIC_COORDS_EXPORT_PATH,
  	  thread_data->thread_id);
  launch(thread_data);
  clean_data(thread_data);
  return 0;
}



