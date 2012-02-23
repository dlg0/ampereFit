EXEC = xAmpereFit

OBJ_DIR = obj
MOD_DIR = mod
SRC_DIR = src
F90 = gfortran
# Create the mod files 1st
${MOD_DIR}/constants.mod: ${OBJ_DIR}/constants.o ${SRC_DIR}/constants.f90
	${F90} -c ${SRC_DIR}/constants.f90 -o ${MOD_DIR}/constants.mod
${OBJ_DIR}/constants.o: ${SRC_DIR}/constants.f90
	${F90} -c ${SRC_DIR}/constants.f90 -o ${OBJ_DIR}/constants.o
${MOD_DIR}/cnvtime.mod: ${OBJ_DIR}/cnvtime.o ${SRC_DIR}/cnvtime.f90
	${F90} -c ${SRC_DIR}/cnvtime.f90 -o ${MOD_DIR}/cnvtime.mod
${OBJ_DIR}/cnvtime.o: ${SRC_DIR}/cnvtime.f90
	${F90} -c ${SRC_DIR}/cnvtime.f90 -o ${OBJ_DIR}/cnvtime.o
# Remove *.obj, *.mod and the exe file
clean:
	@echo ${OBJECTS}
	rm obj/*.o mod/*.mod 
