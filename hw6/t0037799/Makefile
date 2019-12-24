 
SRC_DIR = src
CXXFILES = $(wildcard $(SRC_DIR)/*.cc)
CXXFLAGS = -O3 -Wall -Wextra -Werror 
ifeq ($(shell uname -r), 5.3.6-arch1-1-ARCH)
MKL_INC = /opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include
MKL_LIB_DIR = /opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64
else
MKL_INC = /home/ubuntu/opt/conda/include
MKL_LIB_DIR = /home/ubuntu/opt/conda/lib
endif

TARGET = _matrix.so

all:  $(TARGET)

$(TARGET) : $(CXXFILES)
	$(CXX) -o $@ -fPIC -shared -Iinclude -I$(MKL_INC) `python3 -m pybind11 --includes` $(CXXFLAGS)  $^  -L$(MKL_LIB_DIR) -lmkl_rt -lpthread -ldl

.PHONY: clean
clean:
	rm -rf  $(TARGET)
