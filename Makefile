
libinterpol.so: internal_interpolation.hpp interpolation.cpp interpolation.hpp 
	g++  -O3 -std=c++11 -fvisibility=hidden -fPIC -shared -Wall -Werror interpolation.cpp -o libinterpol.so



test_main: libinterpol.so test_interpolation.cpp test_interpolation.cpp
	g++ -O3 -std=c++11 -Wall -Werror test_interpolation.cpp  -Wl,-rpath=$(PWD) -L/$(PWD) -linterpol -o test_main	 

run-test: test_main
	@rm -f data*.out
	@maxima  --very-quiet --batch-string "load(\"generate_testdata.mc\")$$" >/dev/null
	@./test_main 2>err.out
	@rm -f data*.out

run-test-symbolic:
	@maxima  --very-quiet --batch-string "load(\"verify_symbolic.mc\")$$"


.PHONY: clean

clean:
	rm -f test_main libinterpol.so
	rm -f *.out
