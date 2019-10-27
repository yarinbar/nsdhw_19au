
CC=g++
FLAGS=-Wall -std=c++11 -shared -Werror -fPIC -pedantic-errors `python3 -m pybind11 --includes`
OBJS=matrix.hpp

all:
	${CC} ${FLAGS} ${OBJS} -o _matrix`python3-config --extension-suffix`


test:
	python3 q2_test.py


clean:
	echo Cleaning
	rm -f multiply.o matrix.o
	rm -f multiply`python3-config --extension-suffix` __pycache__
	${CC} ${FLAGS} ${OBJS} -o _matrix`python3-config --extension-suffix`


