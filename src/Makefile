##
## Copyright 2017 Fabon Dzogang
##
## This file is part of vanilla-cspectra.
## vanilla-cspectra is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
## vanilla-cspectra is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with vanilla-cspectra.  If not, see <http://www.gnu.org/licenses/>.
##

include ../Makefile.rules

SRC=    cspectra.c			\
	read_csv.c			\
	build_spectrum.c		\
	ioutils.c			\
	run_thread.c			\
	mem_utils.c			\
	gumbel.c			\
	normal.c

HEADER= build_spectrum.h		\
	read_csv.h			\
	ioutils.h			\
	run_thread.h			\
	mem_utils.h			\
	config.h			\
	gumbel.h			\
	normal.h

DEBUG_FLAG=""

OBJ=${SRC:.c=.o}

all: $(OUT)
	chmod 700 $(OUT) && cp $(OUT) ../

$(OBJ): $(HEADER)

.c.o : $< $(SRC) $(HEADER)
	$(CC) -I$(INCLUDES) $(CFLAGS) -std=c99 -Wall -Wextra -Wformat-security -Ofast -D_FORTIFY_SOURCE=2  -fstack-protector-all -D_GNU_SOURCE -c $< -o $@

$(OUT): $(OBJ)
	$(CC) -lpthread $(OBJ)  $(LIBS) -L$(LIBS_PATH) -o $(OUT) #-Wl,--verbose

clean:
	$(CLEAN)

distclean: clean
	rm -f ../$(OUT) $(OUT)



