MACHINE = greendl
EXEC = xAmpereFit.${MACHINE}

OBJ_DIR = obj
MOD_DIR = mod
SRC_DIR = src

F90 = /opt/local/bin/gfortran-mp-4.4
FFLAGS	= -J${MOD_DIR} 
DEBUG_OPTS = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid -fbounds-check
INCLUDES = 
LINK =

NETCDF_DIR = /opt/local
LINK += -L${NETCDF_DIR}/lib -lnetcdf -lnetcdf_c++ -lnetcdff
INCLUDES += -I${NETCDF_DIR}/include

OBJECTS := $(patsubst src/%,obj/%.o,$(basename $(wildcard src/*)))

${EXEC}: ${OBJECTS}
	${F90} ${FFLAGS} ${OBJECTS} -o $@ ${LINK} ${DEBUG_OPTS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${F90} -c ${FFLAGS} $< -o $@ ${INCLUDES} ${DEBUG_OPTS}

clean:
	@echo ${OBJECTS} 
	rm obj/*.o mod/*.mod ${EXEC}


# Dependencies

${OBJ_DIR}/ampereFit.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/ampFit_data.o \
	${OBJ_DIR}/dlg.o \
	${OBJ_DIR}/spherHarmFns.o

${OBJ_DIR}/ampFit_data.o: \
	${OBJ_DIR}/constants.o \
	${OBJ_DIR}/ampFit_nameList.o \
	${OBJ_DIR}/dlg.o
