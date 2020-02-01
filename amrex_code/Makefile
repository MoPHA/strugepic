# Assuming that amrex is properly in path
# -fno-omit-frame-pointer

LIB_NAME=interpol
LIB_DIR=../lib
DEPS=propagators.o amrex_util.o main_amrex.o
HDEPS=amrex_util.hpp
FLAGS= -std=c++14 -O3 -march=native# -g -fno-omit-frame-pointer
MPICXX=mpic++
LINK_FLAGS=-lamrex -lgfortran -Wl,-rpath=$(LIB_DIR) -L$(LIB_DIR) -l$(LIB_NAME)

all: amrex_main
	
amrex_main: $(DEPS) $(HDEPS)
	$(MPICXX)  $(DEPS) $(FLAGS) $(LINK_FLAGS)  -o amrex_main


%.o : %.cpp
	$(MPICXX) $(FLAGS) $<  -o $@ -c 

.PHONY: clean
clean: co
	rm -f amrex_main
	rm -f *.o
.PHONY: co
co:	
	rm -f plt* -r
	rm -f Backtrace*
