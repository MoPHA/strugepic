# Structure preserving PIC algorithm implemented using AMReX

The code is still a bit messy and in development.


## Algorithm

The implemented core algorithm is from Phys. Plasmas 22, 112504 (2015) Explicit high-order non-canonical symplecticparticle-in-cell algorithms for Vlasov-Maxwell systems <https://doi.org/10.1063/1.4935904> by Jianyuan Xiao, Hong Qin, Jian Liu, Yang He, Ruili Zhang, and Yajuan Sun.

## Dependencies 

### AMReX

Install AMReX with `--dim=3`, MPI and particles enabled (Should be on by default).
A reasonably new compiler should work (at least `-std=c++11` support).


## Interpolation function

By default, the program uses a piecewise 
