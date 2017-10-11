include ../../../superlu/make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../../../superlu/SRC
LIBS	+= -lm

# double real
DF77EXM	= superlu_1d.o superlu_fortran_wrapper.o


all:	SuperLU_1D

SuperLU_1D: $(DF77EXM) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(DF77EXM) $(LIBS) -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

.F.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o SuperLU_1D

