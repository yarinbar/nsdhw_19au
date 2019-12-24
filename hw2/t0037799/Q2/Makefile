angle.so: angle.cc
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` $^ -o $@

.PHONY: clean
clean:
	rm -f angle.so
.PHONY: test
test: test_angle.py angle.so
	python3 $< -v

