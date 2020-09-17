export

CXXFLAGS+= -std=c++14 -O3 -mavx2  -march=native   
MPICXX=mpic++
LNAME=strugepic

SRCS := $(wildcard src/*.cpp src/interpolation/*.cpp )
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INTERPOLATION_FUNC=P8R2

all:
	cd src && $(MAKE)
	cd src/interpolation && $(MAKE)
	mkdir -p lib
	$(MPICXX) -shared  $(OBJS) $(CXXFLAGS) -o lib/lib$(LNAME).so

install:
	mkdir $(PREFIX)/lib
	mkdir $(PREFIX)/include
	cp include/*.hpp $(PREFIX)/include
	cp lib/lib$(LNAME).so $(PREFIX)/lib


.PHONY: test
test:
	cd test && $(MAKE)

.PHONY: test-clean
test-clean:
	cd test && $(MAKE) clean

.PHONY: clean
clean:
	cd src && $(MAKE) clean
	cd src/interpolation && $(MAKE) clean 
.PHONY: clean-all
clean-all: clean test-clean
