libcpp.cpython-36m-x86_64-linux-gnu.so: libcpp.cpp
	g++ -O3 -Wall -shared -std=c++11 \
		libcpp.cpp \
		-fPIC `python3 -m pybind11 --includes` \
		-o libcpp`python3-config --extension-suffix`

test: libcpp.cpython-36m-x86_64-linux-gnu.so
	@pytest test.py

clean:
	rm *.so 
	rm -rf __pycache__