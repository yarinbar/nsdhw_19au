CXX = g++ -Wall -O3

main: main.o Line.o
	$(CXX) $^ -o $@

main.o: main.cpp
	$(CXX) -c $^

Line.o: Line.h Line.cpp
	$(CXX) -c $^

run:	main
	@./main

clean:
	rm main *.o *.h.gch -f
