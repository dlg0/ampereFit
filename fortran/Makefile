MACHINE = greendl
EXEC = xAmpereFit.${MACHINE}

OBJ_DIR = obj
MOD_DIR = mod
SRC_DIR = src

F90 = /opt/local/bin/gfortran-mp-4.4
F77 = /opt/local/bin/gfortran-mp-4.4

F90FLAGS	= -J${MOD_DIR} 
F77FLAGS	= -J${MOD_DIR} 

DEBUG_OPTS = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid -fbounds-check
INCLUDES = 
LINK =

NETCDF_DIR = /opt/local
LINK += -L${NETCDF_DIR}/lib -lnetcdf -lnetcdf_c++ -lnetcdff
INCLUDES += -I${NETCDF_DIR}/include

GSL_DIR = /Users/dg6/code/ampereFit/fortran/gsl
FGSL_DIR = /Users/dg6/code/ampereFit/fortran/fgsl
LINK += -L${FGSL_DIR}/lib -L${GSL_DIR}/lib -lfgsl_gfortran -lgsl -lgslcblas
INCLUDES += -I${FGSL_DIR}/include/gfortran

F95_LAPACK_DIR = /Users/dg6/code/lapack95/LAPACK95
LAPACK_DIR = /Users/dg6/code/lapack/lapack-3.3.0
BLAS_DIR = /Users/dg6/code/gotoBLAS/GotoBLAS2

INCLUDES += -I${F95_LAPACK_DIR}/lapack95_modules
LINK += ${F95_LAPACK_DIR}/lapack95.a ${LAPACK_DIR}/lapack_OSX.a ${BLAS_DIR}/libgoto2.a

OBJECTS := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*)))

${EXEC}: ${OBJECTS}
	${F90} ${F90FLAGS} ${OBJECTS} -o $@ ${LINK} ${DEBUG_OPTS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS}

${OBJ_DIR}/geopack.o: ${SRC_DIR}/geopack/geopack-2008_dp.f
	${F77} -c ${F70FLAGS} $< -o $@ ${DEBUG_OPTS}

clean:
	@echo ${OBJECTS} 
	rm obj/*.o mod/*.mod ${EXEC}


# Dependencies

${OBJ_DIR}/ampereFit.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/spherHarmFns.o \
	${OBJ_DIR}/basisFns.o \
	${OBJ_DIR}/ampFit_solve.o \
	${OBJ_DIR}/write_output.o \
	${OBJ_DIR}/ampFit_shift.o \
	${OBJ_DIR}/ampFit_rotate.o

${OBJ_DIR}/ampFit_data.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/geopack.o

${OBJ_DIR}/write_output.o: \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/constants.o

${OBJ_DIR}/basisFns.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/ampFit_data.o

${OBJ_DIR}/ampFit_solve.o: \
	${OBJ_DIR}/basisFns.o \
	${OBJ_DIR}/constants.o

${OBJ_DIR}/ampFit_shift.o: \
	${OBJ_DIR}/constants.o

${OBJ_DIR}/ampFit_rotate.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_data.o



