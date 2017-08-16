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
#ifndef READ_CSV_H
# define READ_CSV_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "config.h"

typedef enum
  {
    JAN=0,
    FEB,
    MAR,
    APR,
    MAY,
    JUN,
    JUL,
    AUG,
    SEP,
    OCT,
    NOV,
    DEC
  } Month;
int
nextToken(const char*line,
          int *start,
          int *end);

int
parseLine(const char *line,
          char *date_time,
          double *score,
	  char *tok);

int
read_data_from_buff(double *all_data,
		    char*content,
		    struct tm *begin_date,
		    struct tm *end_date);

int
set_time_interval(struct tm *from,
		  const char *from_str,
		  struct tm *to,
		  const char *to_str);

int
add_hour(struct tm *date, int nb_hours);

int
sgets(char *line,
      char **content,
      int max_len);

#endif /*!-- READ_CSV_H --! */



