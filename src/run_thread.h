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
#ifndef RUN_THREAD_H
# define RUN_THREAD_H

#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <fftw3.h>
#include <string.h>
#include "config.h"
#include "build_spectrum.h"
#include "gumbel.h"

void* run_thread(void *data);

typedef struct t_thread_data
{
  int thread_id;
  int process_id;
  unsigned int seed;
  fftw_plan  *plan;
  char **paths;
  char **signal_names;
  char output_path[SUBBUFFSIZE];
  char *export_dir;
  double *vol;
  double *missing;
  int first_signal_index;
  int nb_files;
  int nb_loaded_signals;
  struct tm from;
  struct tm to;
  int real_size;
  int lowest_periodic_component_index;
  int highest_periodic_component_index;
  int export_dir_len;
  FILE *out;
} t_thread_data;
#endif /* ! RUN_THREAD_H */



