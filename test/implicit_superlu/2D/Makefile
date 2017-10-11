include ../../../superlu/make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../../../superlu/SRC
LIBS	+= -lm

EXE = SuperLU_2D

# double real
DF77EXM	= superlu_2d.o superlu_fortran_wrapper.o


all:	$(EXE) run

$(EXE): $(DF77EXM) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(DF77EXM) $(LIBS) -o $@ -g

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

%.o:%.F90
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE) -g

run: 	$(EXE)
	./$<

clean:	
	rm -f *.o $(EXE) fort.*

