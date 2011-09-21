EXEC = xAmpereFit

OBJ_DIR = obj
MOD_DIR = mod
SRC_DIR = src

F90 = gfortran
F77 = gfortran

F90FLAGS= -J${MOD_DIR} #-fno-underscoring
F77FLAGS= -J${MOD_DIR} #-fno-underscoring

DEBUG_OPTS = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid -fbounds-check
INCLUDES = 
LINK =

NETCDF_DIR = /opt/local
NETCDF_MOD = /usr/lib64/gfortran/modules/
F95_LAPACK_DIR = /home/davidg/lapack95/LAPACK95
LAPACK_DIR = /home/davidg/lapack/lapack-3.3.1/lib
BLAS_DIR = /usr/lib64
GSL_DIR = /home/davidg/gsl
FGSL_DIR = /home/davidg/fgsl
AACGM_DIR:=/home/davidg/rst

ifeq ($(findstring dlg-hp,$(HOSTNAME)),dlg-hp)
F90 = /home/dg6/code/gcc/gcc-4.4.5/bin/gfortran
F77 = $(F90)
NETCDF_DIR = /home/dg6/code/netcdf/netcdf_4.1.3_gnu_64
NETCDF_MOD = /home/dg6/code/netcdf/netcdf_4.1.3_gnu_64/include
F95_LAPACK_DIR = /home/dg6/code/lapack95/LAPACK95
LAPACK_DIR = /home/dg6/code/lapack/lapack-3.1.1/lib
BLAS_DIR = /home/dg6/code/goto_blas/GotoBLAS2
GSL_DIR = /home/dg6/code/gsl
FGSL_DIR = /home/dg6/code/fgsl
AACGM_DIR:=/home/dg6/code/rst/rst
endif

LINK += -L${NETCDF_DIR}/lib -lnetcdf -lnetcdf_c++ -lnetcdff
INCLUDES += -I${NETCDF_MOD}

LINK += -L${FGSL_DIR}/lib -L${GSL_DIR}/lib -lfgsl_gfortran -lgsl -lgslcblas
INCLUDES += -I${FGSL_DIR}/include/gfortran

INCLUDES += -I${F95_LAPACK_DIR}/lapack95_modules
ifeq ($(findstring dlg-hp,$(HOSTNAME)),dlg-hp)
	LINK += ${F95_LAPACK_DIR}/lapack95.a ${LAPACK_DIR}/liblapack.a -L${BLAS_DIR} -lgoto2
else
	LINK += ${F95_LAPACK_DIR}/lapack95.a ${LAPACK_DIR}/liblapack.a -L${BLAS_DIR} -l:libblas.so.3
endif
LINK += -L ${AACGM_DIR}/lib -l:libaacgm.1.a -l:libmlt.1.a -l:libastalg.1.a -l:librtime.1.a

OBJECTS := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*)))

${EXEC}: ${OBJECTS}
	${F90} ${F90FLAGS} ${OBJECTS} -o $@ ${LINK} ${DEBUG_OPTS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS}

${OBJ_DIR}/geopack.o: ${SRC_DIR}/geopack/geopack-2008_dp.f
	${F77} -c ${F70FLAGS} $< -o $@ ${DEBUG_OPTS}

clean:
	@echo ${OBJECTS} 
	rm obj/*.o mod/*.mod ${EXEC}


# Dependencies

${OBJ_DIR}/ampFit.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/spherHarmFns.o \
	${OBJ_DIR}/ampFit_basisFns.o \
	${OBJ_DIR}/ampFit_solve.o \
	${OBJ_DIR}/ampFit_write_output.o \
	${OBJ_DIR}/ampFit_shift.o \
	${OBJ_DIR}/ampFit_rotate.o \
	${OBJ_DIR}/f_aacgm.o \
	$(OBJ_DIR)/ampFit_sort.o \
	${OBJ_DIR}/ampFit_geivec_aacgmvec.o

${OBJ_DIR}/ampFit_data.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/geopack.o

${OBJ_DIR}/ampFit_write_output.o: \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/constants.o

${OBJ_DIR}/ampFit_basisFns.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/ampFit_data.o

${OBJ_DIR}/ampFit_solve.o: \
	${OBJ_DIR}/ampFit_basisFns.o \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_rotate.o

${OBJ_DIR}/ampFit_shift.o: \
	${OBJ_DIR}/constants.o \
	$(OBJ_DIR)/ampFit_rotate.o

${OBJ_DIR}/ampFit_rotate.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o

${OBJ_DIR}/dlg.o: \
	${OBJ_DIR}/constants.o

${OBJ_DIR}/ampFit_geivec_aacgmvec.o: \
		${OBJ_DIR}/f_aacgm.o \
		${OBJ_DIR}/constants.o \
		${OBJ_DIR}/ampFit_sort.o



