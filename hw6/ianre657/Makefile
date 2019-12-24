PYBIND = $(shell python3 -m pybind11 --includes)
CXX = g++
CXXFLAGS = --std=c++11 -Wall -O3 -fPIC


MKLROOT ?= ${HOME}/opt/conda
INCLUDE += -I$(MKLROOT)/include
LDFLAGS += -L$(MKLROOT)/lib -lmkl_rt -lpthread -lm -ldl

SUFFIX = $(shell python3-config --extension-suffix )
TARGET = _matrix$(SUFFIX)

all: $(TARGET)

.PHONY : clean run test

Matrix.o: Matrix.cpp Matrix.hpp
	$(CXX) $(CXXFLAGS) ${PYBIND} $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(TARGET): Matrix_pybind.cpp Matrix.o Matrix.hpp
	$(CXX) $(CXXFLAGS) ${PYBIND} $(INCLUDE) --shared `python3 -m pybind11 --includes` $^ -o $@ $(LDFLAGS)

test: $(TARGET)
	pytest -p no:cacheprovider test_Matrix.py

bench:
	./benchmark.py --times 3

clean:
	rm -f *.o $(TARGET)
	rm -rf __pycache__ .pytest_cache
