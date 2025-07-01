include ./Makefile.inc

EXEDIR=../bin
	
.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) $(INC) -c $*.f90

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) -c $*.c

#.SUFFIXES: .o .f90
#.f90.o:
#       $(FC) -c $(FFLAGS) $(INC) $<

%.o: %.f90
	$(FC) $(OPTF) $(INC) -c -o $@ $<

all: compile

compile: TEAT3D

clean:
	rm *.o *.mod

SUB = TEAT3D.o gaussnoise.o subname.o modforward.o subresiduals.o subgradient.o submisfit.o modadjoint.o submodels.o scalegradient.o subsmooth.o
TEAT3D : $(SUB)
	$(F90) $(OPTF) $(SUB) $(LIB) -o $(EXEDIR)/TEAT3D
