# from mkl makefile in nsd repo
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

MKLLINKLINE := \
	${MKLROOT}/lib/libmkl_intel_lp64.${MKLEXT} \
	${MKLROOT}/lib/libmkl_sequential.${MKLEXT} \
	${MKLROOT}/lib/libmkl_core.${MKLEXT} \
	-lpthread -lm -ldl

# for my wsl
# ifeq ($(UNAME_S),Linux)
# MKLROOT ?= /opt/intel/mkl
# MKLEXT ?= so
# CXXFLAGS := -Wl,--no-as-needed
# endif

# MKLLINKLINE := \
# 	${MKLROOT}/lib/intel64/libmkl_intel_lp64.${MKLEXT} \
# 	${MKLROOT}/lib/intel64/libmkl_sequential.${MKLEXT} \
# 	${MKLROOT}/lib/intel64/libmkl_core.${MKLEXT} \
# 	-lpthread -lm -ldl
#

# Compile Info
CC := g++
SRCEXT := cpp
CFLAGS := -O3 -std=c++17 -g -Wall -shared -fPIC -m64
INCDIR = -I include `python3 -m pybind11 --includes` -I${MKLROOT}/include
LIB := ${MKLLINKLINE}

# Source Info, target = cpplox, entry should be in cpplox.cpp
SRCDIR := src
BINDIR := bin
TARGET = _matrix`python3-config --extension-suffix`
SOURCES := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)" ! -iname "$(TARGET).$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BINDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Compile main
$(TARGET): $(OBJECTS)
	@echo "Linking..."
	@echo "$(CC) $^ $(CFLAGS) -o $(TARGET) $(LIB)"; $(CC) $^ $(CFLAGS) -o $(TARGET) $(LIB)
	export PRELOAD_MKL=True

$(BINDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BINDIR)
	@echo "$(CC) $(CFLAGS) $(INCDIR) -c -o $@ $<"; $(CC) $(CFLAGS) $(INCDIR) -c -o $@ $<

# Clean all binary files
clean:
	@echo " Cleaning..."; 
	@echo "$(RM) -r $(BINDIR) $(TARGET)"; $(RM) -r $(BINDIR) $(TARGET)
	@echo "$(RM) -r $(TESTBINDIR)"; $(RM) -r $(TESTBINDIR)

test: $(TARGET)
	@echo " Testing...";
	pytest matrix_unittest.py -v -s

bench: $(TARGET)
	@echo " Procuding bench mark...";
	chmod a+x ./matrix_bench.py
	./matrix_bench.py

# Declare clean as utility, not a file
.PHONY: clean test