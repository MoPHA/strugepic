export


CPU_CXX= -std=c++14 -O3 -mavx2  -march=native



STRUGEPIC=$(shell readlink -f .)
AMREX_LIBRARY_HOME ?= ${AMREX_INSTALL_ROOT}
LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc )
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc | sed 's/-pthread//g' | sed 's/-Wl,-rpath,/--linker-options=-rpath,/g' | sed 's/-Wl,-rpath -Wl,/--linker-options=-rpath,/g')

AMREX_CFLAGS := -I$(INCDIR)  -I $(STRUGEPIC)/include $(COMPILE_CPP_FLAGS)
AMREX_LFLAGS := -L $(STRUGEPIC)/lib -lstrugepic -L$(LIBDIR) $(COMPILE_LIB_FLAGS)   

ifeq ($(BUILD),GPU)
  CXX=nvcc
  LFLAGS+= -lib
  CXXFLAGS=$(AMREX_CFLAGS) #Mahti workaround -ccbin=/appl/spack/v014/install-tree/gcc-4.8.5/gcc-9.3.0-3cdxud/bin/g++
  LNAME=libstrugepic.a
  AMREX_LFLAGS+=--linker-options=-rpath,$(STRUGEPIC)/lib
else
  CXX=g++
  CXXFLAGS+= $(CPU_CXX) -fPIC -c 
  LNAME=libstrugepic.so
  LFLAGS= -shared
  AMREX_LFLAGS+= -Wl,-rpath=$(STRUGEPIC)/lib
endif






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
	cp lib/$(LNAME)* $(PREFIX)/lib


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
