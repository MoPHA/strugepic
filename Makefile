
all: lib test

lib: 
	g++  -O3 -std=c++11 -fvisibility=hidden -fPIC -shared -Wall -Werror interpolation.cpp -o libinterpol.so

test:
	g++ -O3 -std=c++11 -Wall -Werror test_interpolation.cpp  -Wl,-rpath=$(PWD) -L/$(PWD) -linterpol -o test_main	 

run-test: test_main
	@rm -f data*.out
	maxima  --very-quiet --batch-string "load(\"generate_testdata.mc\")$$"
	./test_main
	@rm -f data*.out

run-test-symbolic:
	maxima  --very-quiet --batch-string "load(\"verify_symbolic.mc\")$$"


.PHONY: clean

clean:
	rm -f test_main libinterpol.so
	rm -f data*.out
