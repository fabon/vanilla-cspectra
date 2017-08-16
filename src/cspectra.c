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
#include <stdlib.h>
#include <sys/time.h>
#include <fftw3.h>
#include <pthread.h>
#include <ctype.h>
#include <math.h>
#include "read_csv.h"
#include "build_spectrum.h"
#include "ioutils.h"
#include "mem_utils.h"
#include "config.h"
#include "run_thread.h"

static void
usage()
{
  fprintf(stderr,
          "./vanilla_cspectra_beta < NB_THREADS > < exportdir > < volpath > < inputdir >\n");
}

static int
load_volume(const char *volume_path,
	    char *output_path,
	    double *vol,
	    double *missing,
	    int export_dir_len)
{
  char *data_buff_str = (char*)malloc(100*_1MB*sizeof(char));
  char *head = data_buff_str;
  append_file_content_to_buff(data_buff_str,
			      0,
			      volume_path);
  char line[BUFFSIZE];
  int max_len = 0;

  while (max_len < MAX_LEN &&
	 1 < sgets(line, &data_buff_str, BUFFSIZE))
    {
      for (unsigned i = 0; i < strlen(line); ++i)
	if (!isdigit(line[i]) && line[i] != '+' && line[i] != 'e' && line[i] != '.')
	  return 0;
      double v = atof(line);
      vol[max_len] = (v <= MIN_SAMPLE_SIZE) ? NAN : v;
      missing[max_len] = (v <= MIN_SAMPLE_SIZE) ? 1 : 0;
      ++max_len;
    }

  interpolate_signal(vol, max_len);
  char *output_buffer = (char*)malloc(NFFT*SUBBUFFSIZE*sizeof(char));
  save_series(vol, "_loaded_volume_raw_", output_buffer, output_path, export_dir_len, max_len);
  smooth_signal(vol, max_len, SMOOTHING_WINDOW);

  memset(output_buffer, 0, NFFT*SUBBUFFSIZE*sizeof(char));
  save_series(vol, "_loaded_volume_", output_buffer, output_path, export_dir_len, max_len);
  free(head);
  free(output_buffer);
  return 1;
}

static int
assign_user_variables(char **argv,
                      int *nb_threads,
                      int *export_dir_len,
                      const char **export_dir,
                      const char **input_dir,
                      const char **volume_path,
                      struct tm *from,
                      struct tm *to)
{
  /*
  ** nb processing threads
  */
  *nb_threads = atoi(argv[1]);

  /*
  ** directory where results are stored
  */
  *export_dir = argv[2];
  if (!file_exists(*export_dir))
    {
      fprintf(stderr, "ERROR: export dir [%s] does not exist!\n", *export_dir);
      return 6;
    }

  *export_dir_len = strlen(*export_dir);
  if (*export_dir_len > SUBBUFFSIZE)
    {
      fprintf(stderr, "ERROR: export dir path [%s] is too long!\n", *export_dir);
      return 7;
    }

  /*
  ** path for the global volume timeline
  */
  *volume_path = argv[3];
  fprintf(stdout, "Volume timeline path: %s", *volume_path);

  /*
  ** paths of the input directory on disk
  */
  *input_dir = argv[4];
  fprintf(stdout, "\t input_dir: %s", *input_dir);

  /*
  ** begin and end date of the processed time window
  */
  set_time_interval(from, START_DATE,to, LAST_DATE);

  char tobuffer[SUBBUFFSIZE];
  char frombuffer[SUBBUFFSIZE];
  msg(stdout, "Time interval:", INFO);
  strftime(frombuffer, sizeof(frombuffer), DATE_FORMAT, from);
  strftime(tobuffer, sizeof(tobuffer), DATE_FORMAT, to);
  fprintf(stdout, "from %s to %s", frombuffer, tobuffer);
  return 0;
}

static int
init_environment(t_spectral_coords ***spectral_coords,
                 double **vol,
		 double **missing,
                 int *real_size,
                 int *lowest_periodic_component_index,
                 int *highest_periodic_component_index,
		 const char *export_dir)
{
  *vol = (double*)fft_allocate_and_zero(sizeof(double) * MAX_LEN);
  *missing = (double*)fft_allocate_and_zero(sizeof(double) * MAX_LEN);

  *spectral_coords = (t_spectral_coords**)allocate_and_zero(NB_FILES_LIMIT*sizeof(t_spectral_coords*));
  *real_size = NFFT/2+1;

  printf("EXPORT DIR [%s]\n", export_dir);
  printf("NFFT: %i\n", NFFT);
  printf("NPERMS: %i\n", NPERMS);
  printf("LONGEST_PERIOD: %i\n", LONGEST_PERIOD);
  printf("MINIMUM_PERIOD: %i\n", MINIMUM_PERIOD);
  printf("REAL SIZE: %i\n", *real_size);
  printf("SMOOTHING_WINDOW: %i\n", SMOOTHING_WINDOW);

  *lowest_periodic_component_index = ceil(NFFT/(float)(NFFT/2.0));
  *highest_periodic_component_index = ceil(NFFT/(float)MINIMUM_PERIOD);
  printf("lowest index: %i\n", *lowest_periodic_component_index);
  printf("highest index: %i\n", *highest_periodic_component_index);
  if (*highest_periodic_component_index > *real_size)
    {
      fprintf(stderr, "highest period %i is past the Nyquist frequency %i\n",
	      *highest_periodic_component_index,
	      *real_size);
      return 1;
    }
  return 0;
}

static void
clean_environment(t_spectral_coords **spectral_coords,
                  char **paths,
                  char **signal_names,
                  double *vol,
		  double *missing,
                  int nb_files)
{
  if (vol)
    fftw_free(vol);
  if (missing)
    fftw_free(missing);

  if (paths)
    {
      for (int i = 0;
           i < nb_files;
           i++)
        if (paths[i])
          free(paths[i]);
      free(paths);
    }
  if (signal_names)
    {
      for (int i = 0;
           i < nb_files;
           i++)
        if (signal_names[i])
          free(signal_names[i]);
      free(signal_names);
    }

  if (spectral_coords)
    {
      for (int i = 0;
           i < nb_files;
           i++)
        if (spectral_coords[i])
          free(spectral_coords[i]);
      free(spectral_coords);
    }
}

int main(int argc, char **argv)
{
  if (argc < 5)
    {
      usage();
      return 3;
    }

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  /**
   ** Read user/run settings
   **/

  printf("Reading user settings...\n");
  int nb_threads = 0;
  const char *export_dir = 0;
  int export_dir_len = 0;
  const char *input_dir = 0;
  const char *volume_path = 0;
  struct tm from, to;
  int code = 0;
  if (0 != (code = assign_user_variables(argv,
                                         &nb_threads,
                                         &export_dir_len,
                                         &export_dir,
                                         &input_dir,
                                         &volume_path,
                                         &from,
                                         &to)))
    return code;
  printf("Done reading user settings.\n");

  /**
   ** Traverse input dir and collect paths
   **/

  fprintf(stdout, "Traversing input directory %s ...\n", input_dir);
  char **paths = 0;
  char **signal_names = 0;
  int nb_files = 0;
  if (!traverse_dumps_hierarchy(input_dir,
                                &paths,
                                &signal_names,
                                &nb_files))
    return 4;
  printf("Done traversing.\n");

  /**
   ** Init computing environement
   **/

  printf("Initialising env...\n");
  double *vol = 0;
  double *missing = 0;
  t_spectral_coords **spectral_coords = 0;
  int real_size = 0;
  int lowest_periodic_component_index = 0;
  int highest_periodic_component_index = 0;
  init_environment(&spectral_coords,
                   &vol,
		   &missing,
                   &real_size,
                   &lowest_periodic_component_index,
		   &highest_periodic_component_index,
		   export_dir);
  printf("Done init.\n");

  /**
   ** Process Volume series for normalisation
   **/
  char output_path[SUBBUFFSIZE];
  memcpy(output_path, export_dir, export_dir_len*sizeof(char));
  output_path[export_dir_len] = '/';
  output_path[export_dir_len + 1] = '\0';
  fprintf(stdout, "Processing volume [%s]...\n", volume_path);

  if (!load_volume(volume_path,
		   output_path,
		   vol,
		   missing,
		   export_dir_len))
    {
      fprintf(stderr, "Cannot load vol[%s]\n", volume_path);
      return 4;
    }
  printf("Done volume.\n");

  struct timeval gstart;
  gettimeofday(&gstart, NULL);

  /**
   ** Prepare FFTW plan
   **/

  double *x = (double*) fft_allocate_and_zero(sizeof(double) * NFFT);
  fftw_complex *X = (fftw_complex*) fft_allocate_and_zero(sizeof(fftw_complex) * NFFT);

  fftw_plan  plan[nb_threads];
  struct timeval start,end;

  for (int i = 0; i < nb_threads; i++)
    {
      fprintf(stdout, "Preparing plan for thread-%i\n", i);
      gettimeofday(&start, NULL);
      plan[i] = fftw_plan_dft_r2c_1d(NFFT,
				     /* x[i], */
				     /* X[i], */
				     x,
				     X,
				     FFTW_MEASURE);

      gettimeofday(&end, NULL);
      fprintf(stdout, "Done. Elapsed time: %ldmsec(s)",
  	      1000*(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000);
    }

  /**
   ** Init and Run thread pool
   **/

  pthread_t thread[nb_threads];
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);


  int nb_signals_per_threads = ceil((double)nb_files/(double)nb_threads);
  t_thread_data *data = (t_thread_data *)malloc(nb_threads * sizeof(t_thread_data));
  for(int t=0; t<nb_threads; t++)
    {
      printf("Main: creating thread %i\n", t);

      data[t].plan = &(plan[t]);
      data[t].thread_id = t;
      data[t].nb_files = nb_files;
      memcpy(&(data[t]).from, &from, sizeof(struct tm));
      memcpy(&(data[t]).to, &to, sizeof(struct tm));
      data[t].seed = time(NULL);
      data[t].real_size = real_size;
      memcpy(data[t].output_path, output_path, export_dir_len*sizeof(char));
      data[t].export_dir = strdup(export_dir);
      data[t].lowest_periodic_component_index = lowest_periodic_component_index;
      data[t].highest_periodic_component_index = highest_periodic_component_index;
      data[t].export_dir_len = export_dir_len;
      data[t].vol = (double*)fft_allocate_and_cpy(vol, MAX_LEN*sizeof(double));
      data[t].missing = (double*)fft_allocate_and_cpy(missing, MAX_LEN*sizeof(double));

      int first_signal_index = t*nb_signals_per_threads;
      int nb_signals_to_load =
        (first_signal_index + nb_signals_per_threads) < nb_files ? nb_signals_per_threads : nb_files - first_signal_index;

      data[t].paths = allocate_and_cpy_string_array(paths,
						    first_signal_index,
						    nb_signals_to_load);
      fprintf(stdout, "first signal index %i\n", first_signal_index);
      fprintf(stdout, "nb_signals to load %i\n", nb_signals_to_load);
      data[t].signal_names = allocate_and_cpy_string_array(signal_names,
							   first_signal_index,
							   nb_signals_to_load);
      data[t].nb_loaded_signals = nb_signals_to_load;
      data[t].signal_names = allocate_and_cpy_string_array(signal_names,
							   first_signal_index,
							   nb_signals_to_load);

      data[t].first_signal_index = first_signal_index;
      fprintf(stdout, "Allocated %i signals to thread-%i\n", nb_signals_to_load, t);

      int rc = pthread_create(&thread[t], &attr, run_thread, (void *)&(data[t]));
      if (rc)
        {
          printf("ERROR. posix code: %d\n", rc);
          exit(-1);
        }
    }

  pthread_attr_destroy(&attr);
  for(int t=0;
      t<nb_threads;
      t++)
    {
      void *status;
      int rc = pthread_join(thread[t], &status);
      if (rc)
        {
          printf("ERROR. posix code: %d\n", rc);
          exit(-1);
        }
      printf("Main: joined thread-%i, status:%ldf\n",t,(long)status);
    }

  struct timeval gend;
  gettimeofday(&gend, NULL);

  char buff[SUBBUFFSIZE];
  sprintf(buff, "Elapsed time: %ldmsec(s)",
          1000*(gend.tv_sec - gstart.tv_sec) + (gend.tv_usec - gstart.tv_usec)/1000);
  msg(stdout, buff, INFO);

  /**
   ** Merge results
   **/

  sprintf(buff, "%s/%s_merged.csv",
  	  export_dir,
  	  PERIODIC_COORDS_EXPORT_PATH);
  FILE *f = fopen(buff, "w");
  fprintf(f, "%s\n", HEADERS);
  printf("Merging results...\n");
  for (int t = 0; t < nb_threads; t++)
    {
      sprintf(buff, "%s/%s_%i.csv",
	      export_dir,
	      PERIODIC_COORDS_EXPORT_PATH,
	      t);
      fprintf(stdout, "%s\n", buff);
      FILE *fin = fopen(buff, "rb");
      if (!fin)
	{
	  fprintf(stdout, "Could not open [%s]...\n", buff);
	  continue;
	}
      fseek(fin, 0, SEEK_END);
      long fsize = ftell(fin);
      fseek(fin, 0, SEEK_SET);
      char *content = (char*)malloc(fsize + 1);
      if (!fread(content, fsize, 1, fin))
	printf("The Merge failed results are accessible in the periodics_coords_* files\n");
      content[fsize] = '\0';
      fclose(fin);
      if (!(fprintf(f, "%s", content)))
	printf("The Merge failed results are accessible in the periodics_coords_* files\n");
      free(content);
    }
  printf("Done merging.\n");
  fflush(f);
  fclose(f);

  for (int i = 0;i <10;i++)
    {

      sprintf(buff, "%s/%s_merged_remain_%i.csv",
	      export_dir,
	      PERIODIC_COORDS_EXPORT_PATH,
	      i);
      FILE *f = fopen(buff, "w");
      fprintf(f, "%s\n", HEADERS);
      printf("Merging results...\n");
      for (int t = 0; t < nb_threads; t++)
	{
	  sprintf(buff, "%s/%s_remain_%i_%i.csv",
		  export_dir,
		  PERIODIC_COORDS_EXPORT_PATH,
		  i,
		  t);
	  fprintf(stdout, "%s\n", buff);
	  FILE *fin = fopen(buff, "rb");
	  if (!fin)
	    {
	      fprintf(stdout, "Could not open [%s]...\n", buff);
	      continue;
	    }
	  fseek(fin, 0, SEEK_END);
	  long fsize = ftell(fin);
	  fseek(fin, 0, SEEK_SET);
	  char *content = (char*)malloc(fsize + 1);
	  if (!fread(content, fsize, 1, fin))
	    printf("The Merge failed results are accessible in the periodics_coords_* files\n");
	  content[fsize] = '\0';
	  fclose(fin);
	  if (!(fprintf(f, "%s", content)))
	    printf("The Merge failed results are accessible in the periodics_coords_* files\n");
	  free(content);
	}
      printf("Done merging.\n");
      fflush(f);
      fclose(f);



    }


  /**
   ** Clean environment
   **/
  for (int i = 0; i < nb_threads; i++)
    {
      fftw_destroy_plan(plan[i]);
    }
  if (data)
    free(data);
  if (x)
    fftw_free(x);
  if (X)
    fftw_free(X);
  clean_environment(spectral_coords,
                    paths,
                    signal_names,
                    vol,
		    missing,
                    nb_files);

  pthread_exit(NULL);
  return 0;
}



