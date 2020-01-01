





## Polynominal evaluation

The methods for evaluating the polynominals are not
optimized! Just using straightforwards Horner's rule

We leave the normalization of the arguments to the simulation side,

## Code

Top-level `include`, `src` and `lib` are for the interpolating function

## Amrex example

The code using Amrex assumes that include paths and library paths 
point to the install location for Amrex. It is also assumed that there are some
mpi-wrappers available. The code is tested using open-mpi 


# Current status

Code is not working, as the particle velocities grown exponentially when running a complete simulation. 
I think this has something to do with the particle self-interaction.

## Following parts of the code has been verified.

### Particle velocity and position pusher in $\Theta_x, \Theta_y$ and $\Theta_z$
- Tested by disabling all the other update functions and placing one charged particle in a constant magnetic field
    - Circular motion is observed in the $xy, yz$ and $zx$ plane, as well as randomly selected combinations of orthogonal $B$ and $v$ vectors.
    - Adding a scaling factor was needed, I think. 
    - Motion with this is invariant to grid selection.

### Pure field propagation in $\Theta_E$ and $\Theta_B$
- Tested by having no particles in the simulation and just running the curl propagators.
    - Initial $B$ field $\sin$ stays as $\sin$ (stading wave) and system energy is conserved, $E$ also has correct form $\cos$
    - `B(i,j,k,1) = sin( 2*(geom.ProbLo(X) + i*geom.CellSize(X) +1)*3.1415962);`
    - Single point perturbation is also stable
    - Forward difference is unstable and diverges

### $E$ field change caused by particles and $E$ field velocity update 
- If there is no $B$ field, harmonic motion should occur for a single particle moving in a straight light (Why though?)
$$ 
    \frac{\partial \mathbf{E}}{\partial t}= -\mathbf{J} = -q\mathbf{v},\quad q\mathbf{E} = m\mathbf{\dot{v}}  \Rightarrow
$$

$$
m\frac{\partial^2\mathbf{v}}{\partial t^2}=\frac{\partial \mathbf{E}}{\partial t}q=-q^2\mathbf{v} 
$$


### Misc
`sgn(p.rdata(2+comp))*coef*res` in e_pusher for sine motion of particle

