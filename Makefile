export

INTERPOLATION_NAME=interpol
INTERPOLATION_PATH=interpolation/lib
INTERPOLATION_RANGE=2
INTERPOLATION_FUNC=poly8_range2


INTERPOLATION_FULL_PATH=$(shell readlink -f $(INTERPOLATION_PATH) )


CXXFLAGS+= -std=c++14 -Ofast -mavx2  -march=native   
MPICXX=mpic++
CXX=g++
LNAME=strugepic



all:
	cd core && $(MAKE)
ifeq ($(INTERPOLATION_PATH),interpolation/lib)
	@echo "Compiling default interpolation library"
	cd interpolation && make
endif

tests:
	cd core && $(MAKE) tests

install:
	cd core && $(MAKE) install
ifeq ($(INTERPOLATION_PATH),interpolation/lib)
	@echo "Installing default interpolation library"
	cp interpolation/lib/lib$(INTERPOLATION_NAME) $(PREFIX)/lib
endif

clean:
	cd core && $(MAKE) clean
	cd interpolation && $(MAKE) clean 

