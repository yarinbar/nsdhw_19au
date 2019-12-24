CXX=g++

hw2q2.cpython-37m-x86_64-linux-gnu.so: hw2q2.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` hw2q2.cpp -o hw2q2`python3-config --extension-suffix`
 
.PHONY: test
test: test.py hw2q2.cpython-37m-x86_64-linux-gnu.so
	python3 $<

.PHONY: clean
clean:
	rm -rf hw2q2.cpython-37m-x86_64-linux-gnu.so