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
include ./Makefile.rules

all:
	cd src/ && make

clean:
	$(CLEAN) && cd src/ && make clean

distclean: clean
	cd src/ && make distclean && cd ../ && rm -f Makefile.rules

dist: distclean
	cd ../ && tar czf cspectra.tar.gz --exclude=*.svn --exclude=*.git* cspectra/ &&	\
	chmod 705 cspectra.tar.gz

debug: clean
	cd src/ && make debug


