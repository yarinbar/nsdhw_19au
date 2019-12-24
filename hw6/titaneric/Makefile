FLAGS ?=

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

CXX = g++
CXXFLAGS := ${CXXFLAGS} \
	-std=c++17 -O3 -m64 \
	${FLAGS}

# For MKL

MKLLINKLINE := ${MKLROOT}/lib/libmkl_intel_lp64.${MKLEXT} \
	${MKLROOT}/lib/libmkl_sequential.${MKLEXT} \
	${MKLROOT}/lib/libmkl_core.${MKLEXT} \
	-lpthread -lm -ldl

INCLUDES = 	-I${MKLROOT}/include

# For pybind11

build: TARGET = _matrix`python3-config --extension-suffix`

build: INCLUDES += `python3 -m pybind11 --includes`

build: CXXFLAGS += -shared -Wall -fPIC

build: CXXFLAGS += -DPYTHON_LIB

build: matrix.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} -o ${TARGET} $< ${MKLLINKLINE}
	export PRELOAD_MKL=True

test: CXXFLAGS += -DPYTHON_TEST
test: build
	bash setup.sh

debug: CXXFLAGS += -g
debug: matrix.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} -o $@ $< ${MKLLINKLINE}

.PHONY: clean
clean:
	rm -rf *.so