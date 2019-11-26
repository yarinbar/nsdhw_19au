
CXX := g++
FLAGS ?= -std=c++17 -O3 -g -m64 -Wall -shared -fPIC
FILE = matrix.cpp
SUFFIX = $(shell python3-config --extension-suffix)
TARGET = _matrix$(SUFFIX)


MKLPATH ?= ${HOME}/opt/conda
INCLUDES += -m64 -I${MKLPATH}/include \
	${MKLPATH}/lib/libmkl_intel_thread.so \
	${MKLPATH}/lib/libmkl_intel_lp64.so \
	${MKLPATH}/lib/libmkl_sequential.so \
	${MKLPATH}/lib/libmkl_avx2.so \
	${MKLPATH}/lib/libmkl_def.so \
	${MKLPATH}/lib/libmkl_core.so

LDFLAGS += -L${MKLPATH}/lib -Wl,--no-as-needed \
	-lmkl_rt -lpthread -lm -ldl

INCLUDES += `python -m pybind11 --includes`
CXXFLAGS += -shared -fPIC

.PHONY: all test bench clean
all: $(TARGET)

$(TARGET): $(FILE)
	$(CXX) $(FLAGS) $(INCLUDES) $^ -o $@ $(LDFLAGS)

test: $(TARGET) test.py
	pytest -sv test.py

bench: $(TARGET) bench.py
	python3 bench.py > performance.txt

clean:
	rm -f *.so
	rm -f performance.txt
