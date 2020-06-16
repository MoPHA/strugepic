LFLAGS=-lamrex -lstrugepic -L ../../lib   -Wl,-rpath=$(INTERPOLATION_FULL_PATH) -L$(INTERPOLATION_FULL_PATH) -l$(INTERPOLATION_NAME) -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib 

all:
	$(MPICXX) $(LFLAGS) $(CXXFLAGS) *.cpp -o single_particle_main
clean: cexe cdata

cexe:
	rm -f single_particle_main	
cdata:
	rm -f Backtrace* *.out *.data
