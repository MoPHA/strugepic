
ifeq ($(INTERPOLATION_FUNC),P8R2)
TEST_WRANGE=2
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
TEST_WRANGE=1
endif

CXXFLAGS+= -DWRANGE=$(TEST_WRANGE)

all: field bernstein langmuir 

field:
	cd field_only && $(CXX) *.cpp -o main.o  $(CXXFLAGS) $(AMREX_CFLAGS) 
	cd field_only && $(CXX) main.o -o main $(AMREX_LFLAGS) 
bernstein:
	cd full && $(CXX) bernstein_main.cpp -o bernstein_main.o  $(CXXFLAGS) $(AMREX_CFLAGS) 
	cd full && $(CXX) bernstein_main.o -o bernstein_main $(AMREX_LFLAGS) 
langmuir:
	cd full && $(CXX) langmuir_main.cpp -o langmuir_main.o  $(CXXFLAGS) $(AMREX_CFLAGS) 
	cd full && $(CXX) langmuir_main.o -o langmuir_main $(AMREX_LFLAGS) 

clean:
	cd field_only && rm -rf main Backtrace* *.out *.data Field_Only_Data_out
	cd full && rm -rf main Backtrace* *.out *.data *Data_out
