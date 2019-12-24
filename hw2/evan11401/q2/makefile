CXX = g++
FLAGS = --std=c++11 -O3 -Wall -shared -fPIC `python3 -m pybind11 --includes` 
PY = python


angle.so: angle.cpp
	$(CXX) $(FLAGS) $< -o $@

.PHONY: clean run
clean:
	rm -rf  *.so
run:
	$(PY) test.py
