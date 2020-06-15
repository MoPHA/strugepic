#ifndef particledefs
#define particledefs
#include "AMReX_Particles.H"
#include<AMReX.H>


// We are not going to support multiple particle structure 
// in this code (I.e all particles have the same data fields)
// So here are the particle structure definitions

// CParticle is a general charged particle with the following info + position.

// mass 
// charge
// vx
// vy
// vz


// Se amrex documentation for more information
// https://amrex-codes.github.io/amrex/docs_html/Particle.html
#define C_NUM_REALS 5 
#define C_NUM_INTS 0 
#define C_NUM_SOA_REALS 0 
#define C_NUM_SOA_INTS 0 

// Pos variable indices
#define X 0
#define Y 1
#define Z 2

// 
#define M 0
#define Q 1
#define VX 2
#define VY 3
#define VZ 4

// Constant for non-periodic particle boundary
// Particles are reflected before the boundary as
// not to affect / be affected by the boundary.
// =2 would mean that particles are reflected 2 cells before the
// actual boundary
// This would affect the segment calculation and the position calculations
#define PART_BOUND 5 

// 



typedef amrex::ParIter< C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParIter;
typedef amrex::ParticleContainer<C_NUM_REALS,C_NUM_INTS,C_NUM_SOA_INTS,C_NUM_SOA_REALS> CParticleContainer;
typedef amrex::Particle<C_NUM_REALS,C_NUM_INTS> CParticle;
typedef amrex::ArrayOfStructs<C_NUM_REALS,C_NUM_INTS> CParticles;
typedef amrex::ParticleTile<C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParticleTile;

#endif
