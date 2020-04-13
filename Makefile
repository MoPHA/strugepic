
all: core

interpolation/lib/libinterpol.so:
	cd interpolation && make 

core: interpolation/lib/libinterpol.so 
	cd amrex_code && make

clean:
	cd amrex_code && make clean
	cd interpolation && make clean 

