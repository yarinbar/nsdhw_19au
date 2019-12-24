MKL_INC = /home/ubuntu/opt/conda/include
MKL_LIB_DIR = /home/ubuntu/opt/conda/lib

_matrix.*.so: matrix.cpp
	g++ -O3 -Wall -shared -std=c++11 \
		matrix.cpp \
		-fPIC `python3 -m pybind11 --includes` \
		-o _matrix`python3-config --extension-suffix` \
		-L$(MKL_LIB_DIR) -lmkl_rt -lpthread -ldl \
		-I$(MKL_INC)

test: _matrix.*.so
	@pytest -v test_matrix.py

bench: _matrix.*.so
	@python3 benchmark.py

clean:
	rm *.so __pycache__ .pytest_cache core -rf
