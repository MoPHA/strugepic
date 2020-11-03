export


CPU_CXX= -std=c++14 -O3 -mavx2  -march=native

ifeq ($(BUILD),GPU)
  CXX=nvcc
  CXXFLAGS+= --Werror ext-lambda-captures-this -dc -x cu -std=c++14 -maxrregcount=255 --ptxas-options=-O3 --use_fast_math --expt-relaxed-constexpr --expt-extended-lambda --gpu-architecture=compute_70 --gpu-code=sm_70,compute_70
  CXXFLAGS+= --compiler-options="$(CPU_CXX)"
  LFLAGS+= --gpu-architecture=compute_70 --gpu-code=sm_70,compute_70 -lib
  LNAME=libstrugepic.a
else
  CXX=g++
  CXXFLAGS+= $(CPU_CXX) -fPIC -c 
  LNAME=libstrugepic.so
  LFLAGS= -shared
endif


STRUGEPIC=$(shell readlink -f .)
AMREX_LIBRARY_HOME ?= ${AMREX_INSTALL_ROOT}
LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)

AMREX_CFLAGS := -I$(INCDIR)  -I $(STRUGEPIC)/include
AMREX_LFLAGS := -L$(LIBDIR) $(COMPILE_LIB_FLAGS) -L $(STRUGEPIC)/lib -lstrugepic  --linker-options=-rpath,$(STRUGEPIC)/lib





SRCS := $(wildcard src/*.cpp src/interpolation/*.cpp )
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INTERPOLATION_FUNC=P8R2
#INTERPOLATION_FUNC=PWL
all:
	cd src && $(MAKE)
	cd src/interpolation && $(MAKE)
	mkdir -p lib
	$(CXX) $(LFLAGS) $(OBJS) -o lib/$(LNAME)

install:
	mkdir $(PREFIX)/lib
	mkdir $(PREFIX)/include
	cp include/*.hpp $(PREFIX)/include
	cp lib/lib$(LNAME)* $(PREFIX)/lib


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
	rm -rf lib
.PHONY: clean-all
clean-all: clean test-clean examples-clean
