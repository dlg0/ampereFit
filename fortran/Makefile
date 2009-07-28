MACHINE = dlghp
EXEC = xAmpereFit.${MACHINE}

OBJ = obj
MOD = mod
SRC = src

F90 = ifort
FFLAGS	= -openmp -p -g -warn all -module ${MOD} #-check bounds 
DISLIN = -I${HOME}/dislin/ifc -L${HOME}/dislin -ldislin
NETCDF = -I ${HOME}/netcdf/netcdf_intel/include -L ${HOME}/netcdf/netcdf_intel/lib -lnetcdf -lnetcdf_c++

IMSL_DIR = ${HOME}/imsl/imsl/fnl600/lnxin100e64
IMSL = -L${IMSL_DIR}/lib -Bdynamic -limsl -limslsuperlu -limslscalar -limslblas -limslmpistub -lm -Xlinker -rpath -Xlinker ${IMSL_DIR}/lib -I ${IMSL_DIR}/include

OBJECTS = ${OBJ}/dlg.o ${OBJ}/spherHarmFns.o ${OBJ}/read_nameList.o ${OBJ}/constants.o

ampereFit: ${OBJECTS}
	${F90} ${FFLAGS} ${SRC}/ampereFit.f90 -o ${EXEC} ${OBJECTS} ${IMSL} ${DISLIN} ${NETCDF}

${OBJ}/spherHarmFns.o: ${SRC}/spherHarmFns.f90 ${OBJ}/read_nameList.o
	${F90} ${FFLAGS} -c ${SRC}/spherHarmFns.f90 -o ${OBJ}/spherHarmFns.o ${IMSL} ${DISLIN}

${OBJ}/read_nameList.o: ${SRC}/read_nameList.f90 ${OBJ}/constants.o 
	${F90} ${FFLAGS} -c ${SRC}/read_nameList.f90 -o ${OBJ}/read_nameList.o 

${OBJ}/constants.o: ${SRC}/constants.f90
	${F90} ${FFLAGS} -c ${SRC}/constants.f90 -o ${OBJ}/constants.o 

${OBJ}/dlg.o: ${SRC}/dlg.f90
	${F90} ${FFLAGS} -c ${SRC}/dlg.f90 -o ${OBJ}/dlg.o ${NETCDF}

clean:
	rm obj/*.o mod/*.mod ${EXEC}
