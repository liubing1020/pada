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
CVODE = $../../../CVODE
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
EXS = toyMPSA2

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
