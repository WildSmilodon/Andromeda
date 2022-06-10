PROGRAM = Andromeda
VERSION = 1.6

MODOBJECTS = mod_common.o mod_parallel.o mod_par.o mod_plis.o mod_mesh.o mod_trikotnik.o mod_gauss.o mod_eqns.o mod_profile.o mod_linFlowFields.o mod_singInt.o
OBJECTSFIX = lib_lsqr.o lib_solveFM.o 
OBJECTSF90 = Andromeda.o laplace.o stokes.o lib_ppLine.o lib_mesh.o lib_int.o lib_domain.o\
             lib_sdBem.o lib_util.o  lib_CRS.o lib_Stokes.o lib_StokesFull.o\
		 lib_particle.o lib_StokesLis.o lib_Pressure.o lib_intAna.o

INST = $(HOME)/bin
LINS = ../bin
SRCF = ../src
HOST = $(shell hostname)

# compiler
#FC = gfortran
FC = mpif90

# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o  

# options
FCOPTS = -ffixed-line-length-132  -fcray-pointer -cpp  -DUSE_MPI -I../inc

# libraries
ifneq (,$(findstring hobbit,$(HOST)))
LDOPTO = -lblas -llis -lmpi
else
#  LDOPTO = /usr/lib/x86_64-linux-gnu/blas/libblas.so
LDOPTO = -lblas -llis -lmpi
endif


all: mod fixed f90
	$(FC) -o $(LINS)/$(PROGRAM) $(FCOPTS) $(MODOBJECTS) $(OBJECTSF90) $(OBJECTSFIX) $(LDOPTO)
	cp *mod $(SRCF)
	

mod: $(MODOBJECTS)
fixed: $(OBJECTSFIX)
f90: $(OBJECTSF90)

clean:
	 rm -f *.o *.mod .#*

$(OBJECTSFIX): %.o: $(SRCF)/%.f
	 $(FC) $(FCOPTS) $(OPT) -c $<

$(OBJECTSF90): %.o: $(SRCF)/%.f90
	 $(FC) $(FCOPTS) $(OPT) -c $<

$(MODOBJECTS): %.o: $(SRCF)/%.f90
	 $(FC) $(FCOPTS) $(OPT) -c $<	 


install:
	cp $(LINS)/$(PROGRAM) $(INST)/$(PROGRAM)-$(VERSION)-$(OPDE)
	 
#
# Instalacija BLAS in LAPAC knjiznic
# sudo apt-get update
# sudo apt-get install libblas-dev liblapack-dev
#	 
