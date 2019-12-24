CXX = g++
FLAGS = --std=c++11

main: main.o line.o
	$(CXX) $^ -o $@

line.o: line.cpp line.hpp
	$(CXX) -c $< -o $@ $(FLAGS)

main.o: main.cpp line.hpp
	$(CXX) -c $< -o $@ $(FLAGS)



.PHONY: clean run
clean:
	rm -rf main *.o
run:
	./main