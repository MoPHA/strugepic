ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: energy single

energy:
	cd energy_conservation && $(CXX) $(CXXFLAGS) $(AMREX_CFLAGS) *.cpp -o main.o
	cd energy_conservation && $(CXX) $(AMREX_LFLAGS) main.o -o main

single: 
	cd single_particle && $(CXX) $(CXXFLAGS) $(AMREX_CFLAGS) *.cpp -o main.o
	cd single_particle && $(CXX) $(AMREX_LFLAGS) main.o -o main

clean:
	cd  energy_conservation && rm -f main Backtrace* *.out *.data
	cd single_particle && rm -f main Backtrace* *.out *.data
