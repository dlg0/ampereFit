.SUFFIXES : .o .f90 .f .F90

EXEC := xAmpereFit
RUN := xRunAmpereFit

OBJ_DIR := obj
MOD_DIR := mod
SRC_DIR := src

# If you are not DLG, then do not change this ;)
PARALLEL := 0

INTEL:=1

ifeq (${INTEL},1)
	F90FLAGS:= -module ${MOD_DIR}
	F77FLAGS:= -module ${MOD_DIR}
	DEBUG_OPTS := -g -check bounds
else
	F90FLAGS:= -J${MOD_DIR} #-fno-underscoring
	F77FLAGS:= -J${MOD_DIR} #-fno-underscoring
	DEBUG_OPTS := -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid -fbounds-check
endif

INCLUDES := 
LINK :=
CPP_DIRECTIVES :=

ifeq (${PARALLEL},1)
	F90 := mpif90
	F77 := mpif90
	CPP_DIRECTIVES += -D_parallel_
else
ifeq (${INTEL},1)
	F90 := ifort
	F77 := ifort
	CPP_DIRECTIVES += -D_intel_
else
	F90 := gfortran
	F77 := gfortran
endif
endif


# Default library paths for ampre.newcastle.edu
	NETCDF_DIR := /opt/local
	NETCDF_MOD := /usr/lib64/gfortran/modules/
	F95_LAPACK_DIR := /home/davidg/lapack95/LAPACK95
	F95_LAPACK_LIB:=${F95_LAPACK_DIR}/lapack95.a
	LAPACK_DIR := /home/davidg/lapack/lapack-3.4.0/lib
	LAPACK_LIB:=${LAPACK_DIR}/liblapack.a
	#BLAS_DIR := /usr/lib64
	#BLAS_LIB := -L${BLAS_DIR} -l:libblas.so.3 -lpthread
	BLAS_DIR := /home/davidg/blas/GotoBLAS2
	BLAS_LIB := -L${BLAS_DIR} -l:libgoto2.so

	#BLAS_LIB:=
	#MKLROOT:=/home/davidg/intel/composer_xe_2011_sp1.7.256/mkl
	#F95ROOT:=${MKLROOT}
	#INCLUDES+=${F95ROOT}/include/intel64/ilp64
	#LAPACK_LIB:= $(F95ROOT)/lib/intel64/libmkl_lapack95_ilp64.a  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -fopenmp -lpthread -lm
	#CPP_DIRECTIVES += -I$(F95ROOT)/include/intel64/ilp64  -m64 -I$(MKLROOT)/include

	GSL_DIR := /home/davidg/gsl
	FGSL_DIR := /home/davidg/fgsl
	AACGM_DIR:=/home/davidg/rst

ifeq (${PARALLEL},1)
	SCALAPACK_DIR:=${HOME}/code/scalapack
	BLACS:= \
		${HOME}/code/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a \
		${HOME}/code/blacs/blacs_gnu64/LIB/blacsF77init_MPI-LINUX-0.a \
		${HOME}/code/blacs/blacs_gnu64/LIB/blacs_MPI-LINUX-0.a
endif


# Include other machines
include Makefile.dlg-hp
include Makefile.fitacf

# Netcdf
LINK += -L${NETCDF_DIR}/lib -lnetcdf -lnetcdf_c++ -lnetcdff 
INCLUDES += -I${NETCDF_MOD}

## FGSL (Fortran interface to the GSL)
#LINK += -L${FGSL_DIR}/lib -L${GSL_DIR}/lib -lfgsl_gfortran -lgsl -lgslcblas
#INCLUDES += -I${FGSL_DIR}/include/gfortran

# F95 interface to LAPACK

ifeq (${INTEL},1)
# Taken directly from the intel link line advisor page
LINK +=  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a  $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread -lm 
INCLUDES += -I${F95_LAPACK_DIR}/include/intel64/lp64
else
LINK += ${F95_LAPACK_LIB} ${LAPACK_LIB} ${BLAS_LIB}
INCLUDES += -I${F95_LAPACK_DIR}/lapack95_modules
endif

# AACGM 
ifeq (${INTEL},1)
LINK += -L ${AACGM_DIR}/lib -laacgm -lmlt -lastalg -lrtime
else
LINK += -L ${AACGM_DIR}/lib -l:libaacgm.1.a -l:libmlt.1.a -l:libastalg.1.a -l:librtime.1.a
endif

# ScaLAPACK
ifeq (${PARALLEL},1)
LINK := -L ${SCALAPACK_DIR} -l:libscalapack_gnu_4.3.2_64.a ${BLACS} ${LINK}
endif

OBJECTS := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*)))

${RUN}: ${EXEC} ${OBJ_DIR}/ampFit_run.o
	${F90} ${F90FLAGS} -o $@ ${OBJECTS} ${LINK} ${DEBUG_OPTS} ${CPP_DIRECTIVES}

${EXEC}: ${OBJECTS}
	${F90} ${F90FLAGS} ${OBJECTS} -o $@ ${LINK} ${DEBUG_OPTS} ${CPP_DIRECTIVES}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS} ${CPP_DIRECTIVES}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS} ${CPP_DIRECTIVES}

#${OBJ_DIR}/ampFit_solve.o: ${SRC_DIR}/ampFit_solve.F90
#	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS} ${CPP_DIRECTIVES} -fdefault-integer-8

${OBJ_DIR}/geopack.o: ${SRC_DIR}/geopack/geopack-2008_dp.f
	${F77} -c ${F70FLAGS} $< -o $@ ${DEBUG_OPTS} ${CPP_DIRECTIVES}

clean:
	@echo ${OBJECTS} 
	rm obj/*.o mod/*.mod ${EXEC} ${RUN}

include Makefile.deps

