
CC=g++
FLAGS=-Wall -std=c++11 -shared -Werror -fPIC -pedantic-errors `python3 -m pybind11 --includes`
OBJS=matrix.hpp


UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
MKLROOT ?= /opt/intel/mkl
MKLEXT ?= a
CXXFLAGS :=
endif

ifeq ($(UNAME_S),Linux)
MKLROOT ?= ${HOME}/opt/conda
MKLEXT ?= so
CXXFLAGS := -Wl,--no-as-needed
endif


MKLLINKLINE := \
	${MKLROOT}/lib/libmkl_intel_lp64.${MKLEXT} \
	${MKLROOT}/lib/libmkl_sequential.${MKLEXT} \
	${MKLROOT}/lib/libmkl_core.${MKLEXT} \
	-lpthread -lm -ldl

CXXFLAGS := ${CXXFLAGS} \
	-std=c++17 -O3 -g -m64 \
	-I${MKLROOT}/include \
	${MKLLINKLINE} \
	${FLAGS}
	
all:
	${CC} ${CXXFLAGS} ${OBJS} -o matrix`python3-config --extension-suffix`


test:
	python3 q2_test.py


clean:
	echo Cleaning
	rm -rf *.o *.dSYM

