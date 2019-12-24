.PHONY: test

HEADER   = matrix.hpp

O_DIR    = obj

CXX      = g++
LDFLAGS  = -lm
CXXFLAGS = -O3 -std=c++11 -fPIC -g -Wno-deprecated -Wall -fopenmp
#CXXFLAGS = -O3 -std=c++11 -fPIC -g -Wno-deprecated -Wall -fopenmp -DNOMKL # on arch

MKLROOT ?= $(HOME)/opt/conda
INCLUDE += -I$(MKLROOT)/include
LDFLAGS += -L$(MKLROOT)/lib -lpthread -lm -ldl -lmkl_rt
# LDFLAGS += -L$(MKLROOT)/lib -lpthread -lm -ldl -lcblas # on arch

INCLUDE += $(shell python3-config --includes)
SUFFIX   = $(shell python3-config --extension-suffix)
LDFLAGS += $(shell python3-config --ldflags)
BINDTAR      = _matrix
BINDNAME     = binding
BINDSRCS     = ${BINDNAME:%=%.cc}
BINDLIB      = $(BINDTAR)$(SUFFIX)
BINDALIAS    = $(BINDTAR).so
BINDOBJ  = $(O_DIR)/$(BINDLIB)

$(BINDOBJ): $(O_DIR)
	$(CXX) $(INCLUDE) $(CXXFLAGS) -shared -o $(BINDOBJ) $(BINDSRCS) $(LDFLAGS)
	ln -fs $(BINDOBJ) $(BINDALIAS)

$(O_DIR) :
	mkdir -p $@

$(O_DIR)/%.o: %.cc
	$(CXX) $(INCLUDE) $(CXXFLAGS) -c -o $@ $<  $(LDFLAGS)

bench: $(BINDOBJ) $(BINDSRCS) $(HEADER) benchmark.py
	python benchmark.py

test: $(BINDOBJ) $(HEADER)
	python -m pytest testing.py -v -s

clean:
	rm -f $(TARGET) $(BINDALIAS)
	rm -rf  $(O_DIR) $(BINDOBJ) $(BINDALIAS)

