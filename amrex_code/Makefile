# Assuming that amrex is properly in path


all: amrex_main
	
amrex_main: hello_amrex.cpp
	mpic++ -std=c++14 hello_amrex.cpp -lamrex -o amrex_main

.PHONY: run
run: amrex_main
	./amrex_main


.PHONY: clean
clean: co
	rm -f amrex_main
.PHONY: co
co:	
	rm -f plt_* -r
	rm -f Backtrace*
