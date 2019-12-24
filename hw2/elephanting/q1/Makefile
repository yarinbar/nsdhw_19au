CXX = g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
EXE = main
OUT = output

$(EXE): main.o Line.o
	@$(CXX) $(CXXFLAGS) $^ -o $@

main.o: main.cpp Line.h
	@$(CXX) -c $(CXXFLAGS) $< -o $@

Line.o: Line.cpp Line.h
	@$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:
	@rm -f *o $(EXE) $(OUT)

run: $(EXE)
	@./$(EXE)