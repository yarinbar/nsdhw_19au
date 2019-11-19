
CXX := g++
FLAGS ?= -std=c++17 -O3 -Wall
FILE = matrix.cpp
SUFFIX = $(shell python3-config --extension-suffix)
TARGET = _matrix$(SUFFIX)


MKLROOT ?= ${HOME}/opt/conda
INCLUDES += -m64 -I${MKLROOT}/include
LDFLAGS += -L${MKLROOT}/lib -Wl,--no-as-needed \
	-lmkl_rt -lpthread -lm -ldl

INCLUDES += `python -m pybind11 --includes`
CXXFLAGS += -shared -fPIC

.PHONY: all 
all: $(TARGET)

$(TARGET): $(FILE)
	$(CXX) $(FLAGS) $(INCLUDES) $< -o $(TARGET) $(LDFLAGS)



.PHONY: clean
clean:
	rm -rf $(TARGET)