
MOD_DIR = ../mod
F90FLAGS= -J${MOD_DIR} #-fno-underscoring
DEBUG_OPTS = -g -fbacktrace -fsignaling-nans -ffpe-trap=zero,invalid -fbounds-check
AACGM_DIR:=/home/davidg/rst
LINK += -L ${AACGM_DIR}/lib -l:libaacgm.1.a -l:libmlt.1.a -l:libastalg.1.a -l:librtime.1.a
F90 = gfortran

ifeq ($(findstring dlg-hp,$(HOSTNAME)),dlg-hp)
F90 = /home/dg6/code/gcc/gcc-4.3.2/bin/gfortran
F77 = $(F90)
AACGM_DIR:=/home/dg6/code/superdarn/rst
endif

all: test_aacgm test_time

test_aacgm: test_aacgm.f90
	${F90} ${F90FLAGS} test_aacgm.f90 -o $@ ${LINK} ${DEBUG_OPTS}

test_time: test_superDARN_time.f90
	${F90} ${F90FLAGS} test_superDARN_time.f90 -o $@ ${LINK} ${DEBUG_OPTS}

clean:
	rm test_aacgm test_time
