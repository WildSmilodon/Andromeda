# Program name and version
PROGRAM = Andromeda
VERSION = 1.7

# Directories
INST_DIR = $(HOME)/bin
SRC_DIR = .
BIN_DIR = ../bin
BUILD_DIR = ../build

# Compiler and flags
FC = mpif90


# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o 


FCOPTS = -ffixed-line-length-132 -fcray-pointer -cpp -DUSE_MPI -I../inc
FLAGS = -J$(BUILD_DIR) # Store .mod files in build directory

# Source files
SMODS = mod_common.f90 mod_parallel.f90 mod_par.f90 mod_plis.f90 mod_mesh.f90 \
        mod_trikotnik.f90 mod_gauss.f90 mod_eqns.f90 mod_profile.f90 \
        mod_linFlowFields.f90 mod_singInt.f90 mod_domainMesh.f90 mod_domainData.f90

SFIX = lib_lsqr.f lib_solveFM.f

SF90 = Andromeda.f90 laplace.f90 stokes.f90 lib_ppLine.f90 lib_mesh.f90 \
       lib_int.f90 lib_domain.f90 lib_sdBem.f90 lib_util.f90 lib_CRS.f90 \
       lib_Stokes.f90 lib_StokesFull.f90 lib_particle.f90 lib_StokesLis.f90 \
       lib_Pressure.f90 lib_intAna.f90 lib_meshRead.f90

# Object files
MODOBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SMODS))
FIXOBJ = $(patsubst %.f,$(BUILD_DIR)/%.o,$(SFIX))
F90OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SF90))

# Libraries
#LDOPTO = -lblas -llis
LDOPTO = -lblas  /Users/ravnik/Library/lis/lib/liblis.a


# Default target
all: $(BUILD_DIR) program

# Rule to create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Linking the final executable
program: $(MODOBJ) $(F90OBJ) $(FIXOBJ)
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -o $(BIN_DIR)/$(PROGRAM) $^ $(LDOPTO)

# Pattern rule to compile Fortran source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

# Clean up generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/$(PROGRAM)


install:
	cp $(BIN_DIR)/$(PROGRAM) $(INST_DIR)/$(PROGRAM)-$(VERSION)-$(OPDE)	