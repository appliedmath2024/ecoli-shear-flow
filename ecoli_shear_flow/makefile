# Source files
SRC     = initial_higdon.f90 calcForceMoment_mot2.f90 \
	sqrtm.f90 initial_body2.f90 frame_CSP.f90 \
	StSolve.f90 StSolve_wall.f90 StSolveMk.f90 AngSolve.f90 AngSolve_wall.f90 \
	updateTriad_mot.f90 rotation.f90 selfcontactForce.f90 \
	selfcontactForce_Fla_Body.f90

SRC_MAIN   = $(SRC) Ecoli_OMP.f90

# Module files
MOD_SHORTROD = nrtype.f90 commonval.f90 interface_rod.f90

# Compiler flags
FC      = ifort
DIR     = -L/usr/lib/
LIBS    = -llapack
FLAGS   = -O3 -xHost -qopenmp -fp-model strict -mcmodel=medium

all: shortrod clean

shortrod:
	$(FC) $(FLAGS) $(DIR) -o ecoli-shortrod.out  $(MOD_SHORTROD) $(SRC_MAIN) $(LIBS)

clean: 
	rm *.mod
