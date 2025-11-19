SHELL = /bin/sh
Objects = constants.o potential.o force.o file_io.o dynamic.o main.o
Program = main
FCOMP = gfortran

DEBUG ?= 0

ifeq ($(DEBUG), 1)
FFLAGS = -fcheck=all -Wextra
else
FFLAGS = -O3 -ffast-math -fno-protect-parens -cpp
endif

all : $(Program)
.PHONY : all

$(Program): $(Objects) Makefile
	$(FCOMP) -o $@ $(FFLAGS) $(Objects)

constants.o:	constants.f90 Makefile
	$(FCOMP) $(FFLAGS) -c $<

file_io.o:	file_io.f90 constants.o Makefile
	$(FCOMP) $(FFLAGS) -c $<

potential.o:	potential.f90 constants.o Makefile
	$(FCOMP) $(FFLAGS) -c $<

force.o:		force.f90 constants.o potential.o Makefile
	$(FCOMP) $(FFLAGS) -c $<

dynamic.o:		dynamic.f90 constants.o force.o Makefile
	$(FCOMP) $(FFLAGS) -c $<

main.o:			main.f90 constants.o file_io.o potential.o force.o dynamic.o Makefile
	$(FCOMP) $(FFLAGS) -c $<


.PHONY: clean
clean:
	$(RM) $(Objects) $(Program) *.mod *__genmod.f90 *~

.PHONY: run
run: $(Program)
	./$(Program)