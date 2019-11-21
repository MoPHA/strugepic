
all: lib test

lib: 
	g++  -O3 -std=c++11 -fvisibility=hidden -fPIC -shared -Wall -Werror interpolation.cpp -o libinterpol.so

test:
	g++ -O3 -std=c++11 -Wall -Werror test_interpolation.cpp  -Wl,-rpath=$(PWD) -L/$(PWD) -linterpol -o test_main	 

run-test: test
	rm -f data*.out
	maxima -b symbolic.mc > /dev/null
	./test_main

.PHONY: clean

clean:
	rm -f test_main libinterpol.so
	rm -f data*.out
