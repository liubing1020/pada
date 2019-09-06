######################################################################
### Stochastic Ranking Evolution Strategy                           ##
### Makefile for libsres.a                                          ##
###                                                                 ##
### For ACADEMIC RESEARCH, this is licensed with GPL license        ##
### For COMMERCIAL ACTIVITIES, please contact the authors           ##
###                                                                 ##
### Copyright (C) 2005 Xinglai Ji (jix1@ornl.gov)                   ##
###                                                                 ##
### This program is distributed in the hope that it will be useful, ##
### but WITHOUT ANY WARRANTY; without even the implied warranty of  ##
### MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the    ##
### GNU General Public License for more details.                    ##
###                                                                 ##
### You should have received a copy of the GNU General Public       ##
### License along with is program; if not, write to the Free        ##
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, ##
### MA 02111-1307, USA.                                             ##
###                                                                 ##
### Author: Xinglai Ji (jix1@ornl.gov)                              ##
### Date:   Mar 2, 2005; Mar 3, 2005                                ##
### Organization: Oak Ridge National Laboratory                     ##
######################################################################

######################################################################
### usage:                                                          ##
######################################################################

# general
CC = gcc
CFLAGS = -g -Wall
LD = gcc -w
AR = ar rv
INSTALL = install -c
RANLIB = echo no randlib

INCLUDE = /usr/local/include
LIBRARY = /usr/local/lib

MKDEPEND = makedepend
MKDIR = mkdir -p
CP = cp
RM = rm -rf
# end of general

# for SRES
SRES_INC = ../include
SRES_LIB = ../lib
# end of ESSR
                                                                                
OPTS = -I$(INCLUDE) -I$(SRES_INC)
DEPS = source/sharefunc.c include/sharefunc.h source/ESSRSort.c  \
       include/ESSRSort.h source/ESES.c include/ESES.h
OBJS = source/sharefunc.o source/ESSRSort.o source/ESES.o 
lib lib/libcvode.a: $(OBJS)
	(ar rcv lib/libsres.a $(OBJS); ranlib lib/libsres.a; rm source/*.o)
source/sharefunc.o: $(DEPS) 
	(cd source; $(CC) $(OPTS) -c sharefunc.c)
source/ESSRSort.o: $(DEPS)
	(cd source; $(CC) $(OPTS) -c ESSRSort.c)
source/ESES.o: $(DEPS)
	(cd source; $(CC) $(OPTS) -c ESES.c)

# end of file
