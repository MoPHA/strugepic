LIB_NAME=interpol
LIB_DIR=../../../interpolation/lib
LFLAGS=-lamrex -lstrugepic -L ../../lib   -Wl,-rpath=\$$ORIGIN/$(LIB_DIR) -L$(LIB_DIR) -l$(LIB_NAME) -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib 

all:
	$(MPICXX) $(LFLAGS) $(CXXFLAGS) *.cpp -o single_particle_main
clean:
	rm -f single_particle_main
	rm -f Backtrace*
