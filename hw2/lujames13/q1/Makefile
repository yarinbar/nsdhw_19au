
CC = g++

XXX: main.o Line.o
	$(CC) main.o Line.o -o XXX

main.o: main.cpp
	$(CC) -c main.cpp

Line.o: Line.cpp Line.h
	$(CC) -c Line.cpp

# do not check if output exist or not, and make output if it don't exist.
run: 
	./XXX

clean:
	rm -rf *.o XXX