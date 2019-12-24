
vec2angle.so: vec2angle.cpp
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` vec2angle.cpp -o vec2angle.so

test: vec2angle.so test_vec2angle.py
	python3 -m unittest test_vec2angle.py

clean: vec2angle.so
	rm -rf vec2angle.so __pycache__