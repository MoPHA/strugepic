export
INTERPOLATION_NAME=interpol
INTERPOLATION_PATH=interpolation/lib

INTERPOLATION_RANGE=2
INTERPOLATION_FUNC=poly8_range2

CXXFLAGS+= -std=c++14 -Ofast -mavx2  -march=native   
MPICXX=mpic++
CXX=g++
LNAME=strugepic



all:
	cd core && make
ifeq ($(INTERPOLATION_PATH),interpolation/lib)
	@echo "Compiling default interpolation library"
	cd interpolation && make
endif

tests:
	cd core && make tests
ifeq ($(INTERPOLATION_PATH),interpolation/lib)
	@echo "Compiling default interpolation tests (Maxima required to run)"
	cd interpolation && make tests
endif

install:
	cd core && make install
ifeq ($(INTERPOLATION_PATH),interpolation/lib)
	@echo "Installing default interpolation library"
	cp interpolation/lib/lib$(INTERPOLATION_NAME) $(PREFIX)/lib
endif

clean:
	cd core && make clean
	cd interpolation && make clean 

