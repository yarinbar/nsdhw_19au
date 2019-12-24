.PHONY: clean

build:
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` angle.cpp -o angle`python3-config --extension-suffix`

test: build
	python3 angle_test.py

clean:
	rm *.so