
main: vector.cpython-36m-x86_64-linux-gnu.so

vector.cpython-36m-x86_64-linux-gnu.so : vector.cpp
	g++ -O3 -Wall -shared -std=c++17 -fPIC `python3 -m pybind11 --includes` vector.cpp -o vector`python3-config --extension-suffix`

clean:
	@echo "Cleaning...";
	@echo "$(RM) *.so"; $(RM) *.so

test:
	nosetests -v -s pybind_test.py

.PHONY: clean run