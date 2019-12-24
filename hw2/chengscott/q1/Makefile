all: test

test: test.cpp line.hpp
	g++ -std=c++11 test.cpp -o test

.PHONY: run
run: test
	./test

.PHONY: clean
clean:
	rm -rf test
