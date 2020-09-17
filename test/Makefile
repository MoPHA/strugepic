

LFLAGS= -L ../../lib -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib -lgfortran -lstrugepic -lamrex  

ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: energy single

energy:
	cd energy_conservation && $(MPICXX) *.cpp -o main $(CXXFLAGS) $(LFLAGS)

single: 
	cd single_particle && $(MPICXX) *.cpp -o main $(CXXFLAGS) $(LFLAGS)

clean:
	cd  energy_conservation && rm -f main Backtrace* *.out *.data
	cd single_particle && rm -f main Backtrace* *.out *.data
