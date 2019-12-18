#ifndef particledefs
#define particledefs
#include "AMReX_Particles.H"
#include<AMReX.H>
#include <AMReX_NeighborParticles.H>


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


typedef amrex::ParIter< C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParIter;
typedef amrex::NeighborParticleContainer<C_NUM_REALS,C_NUM_INTS> CParticleContainer;
typedef amrex::ParticleContainer<C_NUM_REALS, C_NUM_INTS, 0, 0>::ParticleVector CNParticles;
typedef amrex::Particle<C_NUM_REALS,C_NUM_INTS> CParticle;
typedef amrex::ArrayOfStructs<C_NUM_REALS,C_NUM_INTS> CParticles;
typedef amrex::ParticleTile<C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParticleTile;

#endif
