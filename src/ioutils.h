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
#ifndef IOUTILS_H
# define IOUTILS_H

#include <sys/types.h>
#include <dirent.h>
#include <time.h>
#include "config.h"

int
traverse_dumps_hierarchy(const char *input_dir,
			 char ***paths,
			 char ***signal_names,
			 int *nb_files);

int
append_file_content_to_buff(char *data_buff_str,
			    int buff_size,
			    const char *path);

void
msg(FILE *stream, const char *msg, t_verbose_lvl lvl);

void
msg_date(FILE *stream, const char *msg, struct tm *date, t_verbose_lvl lvl);

int file_exists(const char *path);

typedef struct t_dump_info
{
  char *dump_name;
  char *output_buffer;
  char *output_path;
  int export_dir_len;
} t_dump_info;


#endif /* !-- IOUTILS_H --! */



