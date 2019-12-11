# Assuming that amrex is properly in path


LIB_NAME=interpol
LIB_DIR=../lib

all: amrex_main
	
amrex_main: hello_amrex.cpp propagators.o propagators.hpp amrex_util.o
	mpic++ -std=c++14 hello_amrex.cpp propagators.o amrex_util.o  -lamrex -Wl,-rpath=$(LIB_DIR) -L$(LIB_DIR) -l$(LIB_NAME) -o amrex_main

propagators.o: propagators.cpp amrex_util.hpp
	mpic++ -std=c++14 -o propagators.o -c propagators.cpp  -lamrex

amrex_util.o: amrex_util.cpp
	mpic++ -std=c++14 -o amrex_util.o -c amrex_util.cpp -lamrex

.PHONY: run
run: amrex_main
	./amrex_main


.PHONY: clean
clean: co
	rm -f amrex_main
	rm -f *.o
.PHONY: co
co:	
	rm -f plt_* -r
	rm -f Backtrace*
