CC = g++
CFLAGS = --std=c++11 -Wall -g
TARGET = main
.PHONY : clean run

$(TARGET): main.cpp Line.o Line.hpp
	$(CC) $(CFLAGS) main.cpp Line.o -o $@

Line.o : Line.cpp Line.hpp 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)

run: $(TARGET) 
	./$(TARGET)