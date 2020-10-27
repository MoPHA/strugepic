export

CXXFLAGS+= -std=c++14 -O3 -mavx2  -march=native   
MPICXX=mpic++




ifeq ($(BUILD),GPU)
  CXX=nvcc -x cu --expt-relaxed-constexpr --expt-extended-lambda --gpu-architecture=compute_70 --gpu-code=sm_70,compute_70
  NVE="
  NVB=--compiler-options="
else
  CXX=g++
endif


LNAME=strugepic

SRCS := $(wildcard src/*.cpp src/interpolation/*.cpp )
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INTERPOLATION_FUNC=P8R2

all:
	cd src && $(MAKE)
	cd src/interpolation && $(MAKE)
	mkdir -p lib
	$(MPICXX) -shared  $(OBJS) -o lib/lib$(LNAME).so

install:
	mkdir $(PREFIX)/lib
	mkdir $(PREFIX)/include
	cp include/*.hpp $(PREFIX)/include
	cp lib/lib$(LNAME).so $(PREFIX)/lib


.PHONY: test
test:
	cd test && $(MAKE)


.PHONY: examples
examples:
	cd examples && $(MAKE)
.PHONY: examples-clean
examples-clean:
	cd examples && $(MAKE) clean

.PHONY: test-clean
test-clean:
	cd test && $(MAKE) clean

.PHONY: clean
clean:
	cd src && $(MAKE) clean
	cd src/interpolation && $(MAKE) clean 
.PHONY: clean-all
clean-all: clean test-clean examples-clean
