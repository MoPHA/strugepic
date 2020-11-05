# StrugePIC


A Structure preserving PIC algorithm implemented using AMReX.


## Algorithm

The implemented core algorithm is from Phys. Plasmas 22, 112504 (2015) Explicit high-order non-canonical symplecticparticle-in-cell algorithms for Vlasov-Maxwell systems <https://doi.org/10.1063/1.4935904> by Jianyuan Xiao, Hong Qin, Jian Liu, Yang He, Ruili Zhang, and Yajuan Sun.

## Dependencies 

### AMReX

Install AMReX with `--dim=3`, MPI and particles enabled (Should be on by default).
A reasonably new compiler should work (at least `-std=c++11` support).
```
./configure --prefix=AMREX_INSTALL_DIR
make -j 4
make install
```

AMReX also has a working spack package.

### MPI

No direct mpi calls, so any mpi implementations which works with AMReX should work.

## Installation

```
git clone https://github.com/MoPHA/strugepic
make -j 4
make install PREFIX=/path/to/install
```
### GPU version
Build with `make BUILD=GPU`, AMReX also has to be built with GPU support (CUDA)
Currently make file is hardcoded to use `nvcc` for compiling gpu. All GPU interface is done through
AMReX, so in theory, HIP could work if the Makefile is edited and AMReX was built with HIP support.


**Note:** Makefile assumes that AMReX is properly in path, i.e. LIBRARY_PATH and CPATH.
Building the library also does not link against AMReX, this only occurs when building the final program.

By default the library uses a piece-wise 8th order polynominal defined on the range [-2,2] as the interpolation function.
`make -j 4 INTERPOLATION_FUNC=PWL` will build the library using a piece-wise linear function defined on [-1,1].
The interpolation function can be overridden during compiling the final program, make sure the range stays the same. 
Interpolation function interface is defined in `include/strugepic_w.hpp`. The default implementations are weak
symbols in the compiled dynamic library.

## Usage

1. Build and install the library 
2. Write your application, see example in `examples/full` 
3. Link and compile against AMReX and library

## Features

- Propagate system with 1,2 or 4:th order sympletic algorithm (skeleton in code for 2n:th order)
- Checkpoint System state at set intervals 
- Write fields and particle density to disk
- Soft plane-wave source.
- MABC boundary condition (currently only in X direction, but easy to modify) 

## Other

- Particles currently carry charge and mass information.
- AMReX default IO in use, every process writes to their own file 
- Open issue on AMReX about 2D domains.

