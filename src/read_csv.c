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
#include <assert.h>
#include "read_csv.h"
#include "ioutils.h"

int
nextToken(const char*line,
          int *start,
          int *end)
{
  if (!line || line[*end] == '\0')
    return 0;

  int i = *start;
  while (line[i] != '\0' && line[i] != ' ' &&
         line[i] != '\n')
    ++i;
  *end = i;
  return (*end - *start) > 0;
}

int
parseLine(const char *line,
          char *date_time,
          double *score,
          char *tok)
{
  int col = 0;
  int start = 0;
  int end = 0;

  char *time = 0;
  *score=0.0;
  int i = 0;

  char buff[SUBBUFFSIZE];

  while (nextToken(line, &start, &end))
    {
      memcpy(tok, &(line[start]), (end-start)*sizeof(char));
      tok[end - start] = 0;

      if (col > TIME_FIELD)
        {
          double val = atof(tok);

          sprintf(buff,
		  "just read [%s]=>[%e] (%ith token start:%i, end%i)",
		  tok, val, i++, start, end);
          msg(stdout, buff, VERBOSE);
          val = isnan(val) ? 0.0 : val;
          *score += atof(tok);
        }
      switch (col) {
      case DATE_FIELD:
        memcpy(date_time, tok, (end-start)*sizeof(char));
        date_time[end-start] = ' ';
        time = date_time + (end-start+1);
        break;
      case TIME_FIELD:
        if (time)
          {
            memcpy(time, tok, (end-start)*sizeof(char));
            time[end-start] = '\0';
          }
        else
          return 0;
        break;
      default:
        break;
      };
      start = end + 1;
      ++col;
    }
  return 1;
}

static int
date_after_or_equal_hour(struct tm *date1,
                         struct tm *date2)
{
  if (!date1 || !date2)
    return 0;
  if (date1->tm_year < date2->tm_year)
    return 0;
  if (date1->tm_year > date2->tm_year)
    return 1;
  if (date1->tm_mon < date2->tm_mon)
    return 0;
  if (date1->tm_mon > date2->tm_mon)
    return 1;
  if (DAY_STEPS)
    return date1->tm_mday >= date2->tm_mday;

  if (date1->tm_mday < date2->tm_mday)
    return 0;
  if (date1->tm_mday > date2->tm_mday)
    return 1;
  return date1->tm_hour >= date2->tm_hour;
}

static int
compare_dates_hour(struct tm *date1,
                   struct tm *date2)
{
  if (!date1 || !date2)
    return 0;

  if (DAY_STEPS)
    return (date1->tm_year == date2->tm_year &&
            date1->tm_mon == date2->tm_mon &&
            date1->tm_mday == date2->tm_mday);

  return (date1->tm_year == date2->tm_year &&
          date1->tm_mon == date2->tm_mon &&
          date1->tm_mday == date2->tm_mday &&
          date1->tm_hour == date2->tm_hour);
}

int
add_hour(struct tm *date, int nb_hours)
{
  int old = 0;
  int new = 0;
  if (DAY_STEPS)
    {
      old = date->tm_mday;
      date->tm_mday += nb_hours;
      new = date->tm_mday;
    }
  else
    {
      old = date->tm_hour;
      date->tm_hour += nb_hours;
      new = date->tm_hour;
    }

  int res=!!mktime(date);
  char buff[SUBBUFFSIZE];
  if (DAY_STEPS)
    sprintf(buff, "inc day[%i]=>[%i]|[%i] (status %i)", old, new, date->tm_mday, res);
  else
    sprintf(buff, "inc hour[%i]=>[%i]|[%i] (status %i)", old, new, date->tm_hour, res);
  msg(stdout, buff, VERBOSE);
  return res;
}

static int
remove_hour(struct tm *date, int nb_hours)
{
  int old = 0;
  int new = 0;
  if (DAY_STEPS)
    {
      old = date->tm_mday;
      date->tm_mday -= nb_hours;
      new = date->tm_mday;
    }
  else
    {
      old = date->tm_hour;
      date->tm_hour -= nb_hours;
      new = date->tm_hour;
    }
  int res=!!mktime(date);
  char buff[SUBBUFFSIZE];
  if (DAY_STEPS)
    sprintf(buff, "remove day[%i]=>[%i]|[%i] (status %i)", old, new, date->tm_mday, res);
  else
    sprintf(buff, "remove hour[%i]=>[%i]|[%i] (status %i)", old, new, date->tm_hour, res);
  msg(stdout, buff, VERBOSE);
  return res;
}

int
set_time_interval(struct tm *from,
                  const char *from_str,
                  struct tm *to,
                  const char *to_str)
{
  strptime(from_str, DATE_FORMAT, from);
  strptime(to_str, DATE_FORMAT, to);
  to->tm_isdst = 0;
  from->tm_isdst = 0;
  return date_after_or_equal_hour(to, from);
}


int
sgets(char *line,
      char **content,
      int max_len)
{
  if (!content || !*content)
    return 0;

  int i = 0;
  for (/*nothing*/;
                  i < max_len &&
                    **content != '\0' &&
                    **content != '\n';
                  i++)
    ++(*content);

  if (**content == '\0')
    return -1; 

  memcpy(line, *content - i, i*sizeof(char));

  line[i] = '\0';
  ++(*content);

  return i+1;
}

static int
fill_time_gaps(double *all_data,
               struct tm *date,
               struct tm *last_date,
               int *hour_index_buff)
{
  char buffer[SUBBUFFSIZE];
  char buff[SUBBUFFSIZE];
  int hour_index = *hour_index_buff;
  int last_hour_index = hour_index;


  while (!compare_dates_hour(date, last_date))
    {
      /**
       ** Skip February 29
       **/

      if (last_date->tm_mon == FEB && last_date->tm_mday == 29)
        {
          add_hour(last_date, 1);
          continue;
        }
      all_data[hour_index++] = 0.0;
      strftime(buffer, sizeof(buffer), DATE_FORMAT, last_date);
      sprintf(buff, "\tstored(nan): %s %e index:%i", buffer, NAN, hour_index - 1);
      msg(stdout, buff, VERBOSE);
      add_hour(last_date, 1);
    }
  *hour_index_buff = hour_index;
  return hour_index - last_hour_index;
}

int
read_data_from_buff(double *all_data,
                    char *content,
                    struct tm *begin_date,
                    struct tm *end_date)
{
  char line[BUFFSIZE];
  double score = 0.0;
  char date_time[SUBBUFFSIZE];
  char tok[SUBBUFFSIZE];
  char buffer[SUBBUFFSIZE];
  struct tm date;
  memset(&date, '\0', sizeof(date));
  struct tm last_date;

  setenv("TZ", "", 1);
  tzset();

  memcpy(&last_date, begin_date, sizeof(*begin_date));
  remove_hour(&last_date, 1);
  int hour_index = 0;
  int lineNumber = 0;
  int res_next = 0;
  msg(stdout, "reading data ...", VERBOSE);

  int currGap = 0;
  while (1 < (res_next = sgets(line, &content, BUFFSIZE)))
    {
      char t[SUBBUFFSIZE];
      msg(stdout, line, VERBOSE);
      sprintf(t, "res next(%i)", res_next);
      msg(stdout, t, VERBOSE);
      msg(stdout, "end", VERBOSE);

      parseLine(line, date_time, &score, tok);
      memset(&date, '\0', sizeof(date));
      date.tm_isdst=0;
      strptime(date_time, DATE_FORMAT, &date);

      /**
       ** Skip February 29
       **/

      if (date.tm_mon == FEB && date.tm_mday == 29)
        continue;

      /**
       ** Start date: skip data before that point
       **/

      if (!date_after_or_equal_hour(&date, begin_date))
        continue;

      /**
       ** End date: gapped already filled with NAN (e.g init_and_nan) exit at this point
       **/

      int store = 1;
      if (date_after_or_equal_hour(&date, end_date))
        {
          char buff[SUBBUFFSIZE];
          msg(stdout, "End date attained:", VERBOSE);
          strftime(buffer, sizeof(buffer), DATE_FORMAT, &date);
          sprintf(buff, "date %s", buffer);
          msg(stdout, buff, VERBOSE);
          strftime(buffer, sizeof(buffer), DATE_FORMAT, end_date);
          sprintf(buff, "end_date %s", buffer);
          msg(stdout, buff, VERBOSE);

	  store = !!compare_dates_hour(&date, end_date);
          memcpy(&date, end_date, sizeof(*end_date));
	  if (!store)
	    add_hour(&date,1);
        }


      char buff[SUBBUFFSIZE];
      sprintf(buff, "parsed: %s %e", date_time, score);
      msg(stdout, buff, VERBOSE);
      msg_date(stdout, "last date:", &last_date, VERBOSE);
      add_hour(&last_date, 1);
      msg_date(stdout, "last date:", &last_date, VERBOSE);
      struct tm save_date;
      memcpy(&save_date, &last_date, sizeof(last_date));
      if (LONGEST_PERIOD <= (currGap = fill_time_gaps(all_data,
						      &date,
						      &last_date,
						      &hour_index)))
	{
	  return -1;
	}
      if (!store)
	break;
      all_data[hour_index++] = isnan(score)? 0.0 : score;
      strftime(buffer, sizeof(buffer), DATE_FORMAT, &date);
      sprintf(buff, "\tstored: %s %e index:%i", buffer, score, hour_index - 1);
      msg(stdout, buff, VERBOSE);
      ++lineNumber;
    }

  /**
   ** End date not attained: fill gap until that point
   **/

  char buff[SUBBUFFSIZE];
  if (!date_after_or_equal_hour(&date, end_date))
    {
      msg(stdout, "Completing missing our until end date:", VERBOSE);
      strftime(buffer, sizeof(buffer), DATE_FORMAT, &date);
      sprintf(buff, "date %s", buffer);
      msg(stdout, buff, VERBOSE);

      strftime(buffer, sizeof(buffer), DATE_FORMAT, end_date);
      sprintf(buff, "end date %s", buffer);
      msg(stdout, buff, VERBOSE);

      add_hour(&date, 1);
      struct tm final_date;
      memcpy(&final_date, end_date, sizeof(*end_date));
      add_hour(&final_date, 1);

      if (LONGEST_PERIOD <= (currGap = fill_time_gaps(all_data,
						      &final_date,
						      &date,
						      &hour_index)))
	{
	  return -1;
	}
    }

  msg(stdout, "done.", VERBOSE);
  return hour_index;
}



