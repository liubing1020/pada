
# general
CC = gcc
CFLAGS = -g -fpic -mcmodel=medium
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

# for SRES
SRES = ..
LIBSRES = libsres.a
SRES_INC = $(SRES)/include
SRES_LIB = $(SRES)/lib
# end of SRES

# combining options
INC_DIRS = -I$(SRES_INC) -I$(INCLUDE)
LIB_DIRS = -L$(SRES_LIB) -L$(LIBRARY) -lsres -lm -lstdc++ -fpic -mcmodel=medium
# end of options

# rules
EXS = m88

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
