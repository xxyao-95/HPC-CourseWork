# HPC-CourseWork
HPC Course Work

The objective of this coursework is to write a parallel numerical code for solving a two-dimensional smoothed
particle hydrodynamic (SPH) formulation of the Navier-Stokes equations. SPH models the behaviour of fluids by
approximating the fluid properties using smooth kernel density functions. Each particle therefore only influences
neighbouring particles within a given radius and properties of the fluid are smoothed between particles.
The position and velocity of each particle are updated at each time step using an explicit time-integration scheme.
We calculate these state variables by updating our approximation to the fluid density, pressure and viscous forces
on each particle and use these to calculate the acceleration of the particles.
