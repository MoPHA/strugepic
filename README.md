# Structure preserving PIC algorithm implemented using AMReX

The code is still a bit messy and in development.


## Algorithm

The implemented core algorithm is from Phys. Plasmas 22, 112504 (2015) Explicit high-order non-canonical symplecticparticle-in-cell algorithms for Vlasov-Maxwell systems_ <https://doi.org/10.1063/1.4935904>
by Jianyuan Xiao, Hong Qin, Jian Liu, Yang He, Ruili Zhang, and Yajuan Sun.

## Installation

Install AMReX with `--dim=3`, MPI and particles enabled (Should be on by default).
A reasonably new compiler should work (at least `-std=c++11` support).
Currently, the main program is in `amrex_code/main_amrex.cpp` , and the idea is just to modify that to suit your needs (A better interface / usage is planned). Just run the makefile
`amrex_code/Makefile` to compile.


## Interpolation
By default The interpolation functions are linked in from `interpolation/lib/libinterpol.so`,
which is compiled by running `interpolation/Makefile` 


