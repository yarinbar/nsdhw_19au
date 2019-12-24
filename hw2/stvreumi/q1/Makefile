CXX=g++

hw2q1.o: hw2q1.cpp
	$(CXX) -Wall -std=c++11 $< -o $@
 
.PHONY: run
run: hw2q1.o
	./hw2q1.o

.PHONY: clean
clean:
	rm -rf *.o