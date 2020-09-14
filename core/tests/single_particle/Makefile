LFLAGS= -L ../../lib   -Wl,-rpath=$(INTERPOLATION_FULL_PATH) -L$(INTERPOLATION_FULL_PATH) -l$(INTERPOLATION_NAME) -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib -lgfortran  -lstrugepic -lamrex

all:
	$(MPICXX)  *.cpp -o single_particle_main $(CXXFLAGS) $(LFLAGS) 
clean: cexe cdata

cexe:
	rm -f single_particle_main	
cdata:
	rm -f Backtrace* *.out *.data
