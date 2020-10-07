

LFLAGS= -L ../../lib -I ../../include -Wl,-rpath=\$$ORIGIN/../../lib -lgfortran -lstrugepic -lamrex  

ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: field bernstein 

field:
	cd field_only && $(MPICXX) *.cpp -o main $(CXXFLAGS) $(LFLAGS)
bernstein:
	cd full && $(MPICXX) *.cpp -o main $(CXXFLAGS) $(LFLAGS)

clean:
	cd field_only && rm -rf main Backtrace* *.out *.data Field_Only_Data_out
	cd full && rm -rf main Backtrace* *.out *.data Bernstein_Data_out
