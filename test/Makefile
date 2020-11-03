ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: energy single

energy:
	cd energy_conservation && $(CXX) *.cpp -o main $(CXXFLAGS) $(AMREX_CFLAGS) $(AMREX_LFLAGS)

single: 
	cd single_particle && $(CXX) *.cpp -o main $(CXXFLAGS) $(AMREX_CFLAGS) $(AMREX_LFLAGS)

clean:
	cd  energy_conservation && rm -f main Backtrace* *.out *.data
	cd single_particle && rm -f main Backtrace* *.out *.data
