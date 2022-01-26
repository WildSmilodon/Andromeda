PROGRAM = Andromeda
VERSION = 1.5.1

MODOBJECTS = mod_common.o mod_mesh.o mod_trikotnik.o mod_gauss.o mod_par.o mod_eqns.o mod_profile.o mod_linFlowFields.o
OBJECTSFIX = lib_lsqr.o lib_solveFM.o 
OBJECTSF90 = Andromeda.o lib_ppLine.o lib_mesh.o lib_int.o lib_domain.o lib_sdBem.o lib_util.o  lib_CRS.o lib_Stokes.o lib_StokesFull.o lib_particle.o

INST = $(HOME)/bin
LINS = ../bin
HOST = $(shell hostname)

# compiler
FC = gfortran

# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o  

# options
FCOPTS = -ffixed-line-length-132  -fcray-pointer 

# Blas library

ifneq (,$(findstring hobbit,$(HOST)))
  LDOPTO = /usr/lib/libblas.so
else
  LDOPTO = /usr/lib/x86_64-linux-gnu/blas/libblas.so
endif


all: mod fixed f90
	$(FC) -o $(LINS)/$(PROGRAM) $(FCOPTS) $(MODOBJECTS) $(OBJECTSF90) $(OBJECTSFIX) $(LDOPTO)
	

mod: $(MODOBJECTS)
fixed: $(OBJECTSFIX)
f90: $(OBJECTSF90)

clean:
	 rm -f *.o *.mod .#*

$(OBJECTSFIX): %.o: %.f
	 $(FC) $(FCOPTS) $(OPT) -c $<

$(OBJECTSF90): %.o: %.f90
	 $(FC) $(FCOPTS) $(OPT) -c $<

$(MODOBJECTS): %.o: %.f90
	 $(FC) $(FCOPTS) $(OPT) -c $<


install:
	cp $(LINS)/$(PROGRAM) $(INST)/$(PROGRAM)-$(VERSION)-$(OPDE)
	 
#
# Instalacija BLAS in LAPAC knjiznic
# sudo apt-get update
# sudo apt-get install libblas-dev liblapack-dev
#	 
