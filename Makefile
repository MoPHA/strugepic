
all: lib test

lib: 
	g++ -std=c++11 -fPIC -shared -Wall -Werror interpolation.cpp -o libinterpol.so

test:
	
