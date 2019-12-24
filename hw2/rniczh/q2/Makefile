CXX = g++

SRC_DIR = src
LIB_DIR = lib

SRC = $(wildcard $(SRC_DIR)/*.cc)
SUFFIX = $(shell python3-config --extension-suffix)
LIB = $(LIB_DIR)/pybind_vector2d$(SUFFIX)

PY_LIBS = $(shell python3 -m pybind11 --includes)
CPPFLAGS += -Iinclude $(PY_LIBS)
CFLAGS += -O3 -shared -Wall -std=c++14 -fPIC

.PHONY: all clean

all: $(LIB)

$(LIB_DIR)/pybind_vector2d$(SUFFIX): $(SRC)
	@mkdir -p lib
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ -o $@

test: all
	ln -f ./lib/*.so ./test/
	pushd ./test && python3 -m pytest -q test.py; popd

clean:
	$(RM) -r $(LIB_DIR)
