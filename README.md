# Inviscid Vortex System Analysis
Tracking trajectories of the vortices in arbitrary 2D vortex systems. The classical RK4 method is deployed for predicting the spatial advance of the vortices with a fixed time step.

# Language and Environment
Written in Fortran90 compatible with the Linux GNU gfortran compiler.

# Input (input.in)
  - NOV: # of vortices in the system
  - X0, Y0: Initial x and y coordinates of each vortex
  - GAMMA: Vortex strength (or circulation) of each vortex
  - NTST: # of steps for a single simulation
  - DT: Fixed time step
  
# Output (history.dat)
Every 50 steps, corresponding time and (x, y) positions of all vortices will be recorded.
