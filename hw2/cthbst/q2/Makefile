mylib.cpython-36m-x86_64-linux-gnu.so: mylib.cpp
	g++ -O3 -Wall -shared -std=c++11 \
		mylib.cpp \
		-fPIC `python3 -m pybind11 --includes` \
		-o mylib`python3-config --extension-suffix`

test: mylib.cpython-36m-x86_64-linux-gnu.so
	@pytest main.py

clean:
	rm *.so __pycache__ .pytest_cache -rf 
