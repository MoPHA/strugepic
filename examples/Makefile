
ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: field bernstein 

field:
	cd field_only && $(CXX) $(CXXFLAGS) $(AMREX_CFLAGS) *.cpp -o main.o
	cd field_only && $(CXX) $(AMREX_LFLAGS) main.o -o main
bernstein:
	cd full && $(CXX) $(CXXFLAGS) $(AMREX_CFLAGS) *.cpp -o main.o
	cd full && $(CXX) $(AMREX_LFLAGS) main.o -o main

clean:
	cd field_only && rm -rf main Backtrace* *.out *.data Field_Only_Data_out
	cd full && rm -rf main Backtrace* *.out *.data Bernstein_Data_out
