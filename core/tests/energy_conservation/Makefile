LFLAGS= -L ../../lib   -Wl,-rpath=$(INTERPOLATION_FULL_PATH) -L$(INTERPOLATION_FULL_PATH) -l$(INTERPOLATION_NAME) -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib -lgfortran -lamrex -lstrugepic 
all:
	$(MPICXX) *.cpp -o energy_main $(CXXFLAGS) $(LFLAGS)
clean: cexe cdata

cexe:
	rm -f energy_main
cdata:
	rm -f Backtrace* *.out *.data
