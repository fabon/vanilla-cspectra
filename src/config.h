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
#ifndef CONFIG_H
# define CONFIG_H

#define MAX_LEN 2556
#define MIN_SAMPLE_SIZE 0

#define START_DATE "01-01-2010 00:00:00"
#define LAST_DATE "31-12-2016 23:59:59"

#define DAY_STEPS 1
#define NFFT 365
#define LONGEST_PERIOD 365
#define MINIMUM_PERIOD 2

#define SMOOTHING_WINDOW 3
#define HIGH_PASS_WINDOW 731

#define _1MB 1000000//bytes
#define _1KB 1000//bytes

#define PI acos(-1.0)
#define EULER exp(1)
#define NB_THREADS 7
#define SCALE_ATAN 3.0
#define NB_FILES_LIMIT 200000
#define NPERMS 0
#define BUFFSIZE 2048
#define SUBBUFFSIZE 1024
#define SPECTRUM_MAX_SIZE (NFFT/2)*SUBBUFFSIZE
#define CSVDELIM " "
#define CSVHEADERS 0
#define HEADERS "NAME,E,T,PHASE,OUTREACH,PVAL,BONFERRONI_PVAL,bandedEnergy,spectrumEnergy"
#define DATE_FORMAT "%d-%m-%Y %H:%M:%S"
#define DATE_FORMAT_LOG "%D %T"
#define EPS_NORMALISATION 1e-6

#define NB_CMPS 1
#define DATE_FIELD 0
#define TIME_FIELD 1
#define COUNT_FIELD 2
#define VERBOSE_LVL WARNING
#define OUTPUT_LVL WARNING
#define PREPOC_LVL RF_ZSCORE
#define PERIODIC_COORDS_EXPORT_PATH "./periodics_coords_"
#define PERIODIC_NAMES_EXPORT_PATH "./periodics_names_"

typedef enum
  {
    VERBOSE=0,
    INFO=1,
    WARNING=2,
    ERROR=3
  } t_verbose_lvl;

typedef enum
  {
    RAW_COUNTS=0,
    RF=1,
    RF_PREZSCORE=2,
    RF_ATAN=3,
    RF_HIGHPASS=4,
    RF_CENTER=5,
    RF_ZSCORE=6,
  } t_prepoc_lvl;

#endif /* !-- CONFIG_H --! */



