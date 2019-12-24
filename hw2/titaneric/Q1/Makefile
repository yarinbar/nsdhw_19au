.PHONY: clean

CXX = g++

line: line.o main.o
	$(CXX) $^ -o $@

line.o: line.cpp
	$(CXX) -c $< -o $@

main.o:	Q1.cpp
	$(CXX) -c $< -o $@

run: line
	./line

clean:
	rm *.o ./line