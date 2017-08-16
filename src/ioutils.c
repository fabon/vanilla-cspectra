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
#include <stdio.h>
#include <string.h>
#include "ioutils.h"
#include <sys/stat.h>

static char*
extract_signal_name(const char *filename)
{
  int name_len = strlen(filename) - 4;
  char *name = (char*)malloc((name_len + 1)*sizeof(char));
  memcpy(name, filename, name_len*sizeof(char));
  name[name_len] = '\0';
  return name;
}

static char *
concat_path(const char *prefix, const char *filename)
{
  int len_prefix = strlen(prefix);
  int len_filename = strlen(filename);
  char *res = (char *)malloc((len_prefix + len_filename + 2) * sizeof(char));
  memcpy(res, prefix, len_prefix*sizeof(char));
  res[len_prefix] = '/';
  memcpy(res + len_prefix + 1, filename, len_filename*sizeof(char));
  res[len_prefix + len_filename + 1] = '\0';
  return res;
}

int
traverse_dumps_hierarchy(const char *input_dir,
			 char ***paths,
			 char ***signal_names,
			 int *nb_files)
{
  DIR *dp;
  struct dirent *ep;

  *nb_files = 0;
  int file_index = 0;
  dp = opendir(input_dir);
  if (dp == 0)
    {
      fprintf(stderr, "could not open the directory [%s]\n", input_dir);
      return 0;
    }
  file_index = 0;
  while ((ep = readdir(dp)))
    if (strstr(ep->d_name, "csv") && !strstr(ep->d_name, "~") && !strstr(ep->d_name, "#"))
      {
	printf("\t%i %s\n", file_index, ep->d_name);
	*signal_names = (char **)realloc(*signal_names, (file_index + 1) * sizeof(char*));
	(*signal_names)[file_index] = extract_signal_name(ep->d_name);
	*paths = (char**)realloc(*paths, (file_index + 1) * sizeof (char*));
	(*paths)[file_index] = concat_path(input_dir, ep->d_name);
	++(*nb_files);
	++file_index;
      }
  closedir(dp);
  printf("************************ %i signals\n", *nb_files);
  return 1;
}

int
append_file_content_to_buff(char *data_buff_str,
			    int buff_size,
			    const char *path)
{
  FILE *f = fopen(path, "r");
  fseek(f, 0, SEEK_END);
  long file_size = ftell(f);
  fseek(f, 0, SEEK_SET);

  size_t pos = fread(data_buff_str + buff_size, file_size, 1, f);
  fclose(f);
  if (pos == 0)
    return 0;
  return buff_size + file_size + 1;
}

void
msg(FILE *stream, const char *msg, t_verbose_lvl lvl)
{
  if (lvl < VERBOSE_LVL)
    return;
  time_t now;
  time(&now);
  struct tm *tm = localtime(&now);
  char buffer[BUFFSIZE];
  switch (lvl)
    {
    case VERBOSE:
      memcpy(buffer, "----  ", 6*sizeof(char));
      break;
    case INFO:
      memcpy(buffer, "info  ", 6*sizeof(char));
      break;
    case WARNING:
      memcpy(buffer, "warn  ", 6*sizeof(char));
      break;
    case ERROR:
      memcpy(buffer, "error ", 6*sizeof(char));
      break;
    default:
      break;
    };
  int pos = strftime(buffer + 6, SUBBUFFSIZE, DATE_FORMAT_LOG, tm);
  buffer[pos + 6] = ':';
  buffer[pos + 7] = ' ';
  buffer[pos + 8] = '\0';

  strncat(buffer, msg, (BUFFSIZE-1) - (pos+8));
  fprintf(stream, "%s\n", buffer);
}

void
msg_date(FILE *stream, const char *message, struct tm *date, t_verbose_lvl lvl)
{
  if (lvl < VERBOSE_LVL)
    return;

  char date_buffer[SUBBUFFSIZE];
  sprintf(date_buffer, "%s day[%i] month[%i] year[%i] hour[%i] min[%i] sec[%i]",
	  message,
	  date->tm_mday,
	  date->tm_mon,
	  date->tm_year,
	  date->tm_hour,
	  date->tm_min,
	  date->tm_sec);
  msg(stream, date_buffer, lvl);
}

int file_exists(const char *path)
{
  struct stat b;
  return stat(path, &b) == 0;
}



