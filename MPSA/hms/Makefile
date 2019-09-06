######################################################################
### Stochastic Ranking Evolution Strategy                           ##
### Makefile for examples                                           ##
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
### Date:   Mar 3, 2005; Mar 11, 2005; Mar 14, 2005;                ##
### Organization: Oak Ridge National Laboratory                     ##
######################################################################

# general
CC = gcc
CFLAGS = -g -Wall -DDSFMT_MEXP=521
LD = g++ -w
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

# for CVODE
CVODE = $(HOME)/pe/CVODE
LIBODE = libcvode.a
ODE_INC = $(CVODE)/include
ODE_LIB = $(CVODE)/lib
ODE_HDRS = $(ODE_INC)llnltyps.h $(ODE_INC)/cvode.h $(ODE_INC)/cvdense.h $(ODE_INC)nvector.h $(ODE_INC)/dense.h
# end of CVODE


# combining options
INC_DIRS =  -I$(ODE_INC) -I$(INCLUDE)
LIB_DIRS =  -L$(ODE_LIB) -L$(LIBRARY)  -lcvode -lm -lstdc++ 
# end of options

# rules
EXS = hmsMPSA2

.SUFFIXES: .C
.C.o: $(ODE_HDRS)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $<

all: $(EXS)

$(EXS): %: %.o
	$(CC) $@.o -o $@ $(LIB_DIRS) 

clean:
	$(RM) $(EXS)
	$(RM) *.o *.bak core

# end of rules

# end of file
