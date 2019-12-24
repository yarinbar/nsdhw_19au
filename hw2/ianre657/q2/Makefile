CXX = g++
CXXFLAGS = --std=c++11 -Wall -O3

TARGET = PointVec.*.so
all: $(TARGET)
.PHONY : clean run test

PointVec.o: PointVec.cpp PointVec.hpp
	$(CXX) $(CXXFLAGS) -c  $< -o $@

$(TARGET): PointVec_pybind.cpp PointVec.o PointVec.hpp
	$(CXX) $(CXXFLAGS) --shared -fPIC `python3 -m pybind11 --includes` $^ -o PointVec`python3-config --extension-suffix`

test: $(TARGET)
	pytest PointVec_test.py

clean:
	rm -f *.o $(TARGET) 
	rm -rf __pycache__ .pytest_cache
