CXX = g++

SRC = matrix.cc
SUFFIX = $(shell python3-config --extension-suffix)
TARGET = _matrix$(SUFFIX)
TEST = test_matrix.py
BENCH = bench.py

# pybind 11 configuration
PY_LIBS = $(shell python3 -m pybind11 --includes)
CPPFLAGS += $(PY_LIBS) -O3 -shared -Wall -std=c++17 -g -m64 -fPIC -fopenmp -mavx

# mkl configuration
MKLROOT ?= $(HOME)/opt/conda
MKLINC = -I$(MKLROOT)/include
LDFLAGS += -L$(MKLROOT)/lib -lmkl_rt -lpthread -lm -ldl

.PHONY: all clean bench test

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(MKLINC) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

test: $(TARGET) $(TEST)
	python3 -m pytest -q $(TEST)

bench: $(TARGET) $(BENCH)	
	python $(BENCH)

clean:
	$(RM) $(TARGET) performance.txt

