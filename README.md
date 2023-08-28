# LAMMPS-LBM-PSM

## About


## Adding this feature to LAMMPS
The current version has been tested and used with "lammps-stable_29Sep2021_update3".
Simply copy all files from the src folder into your LAMMPS src folder and compile your LAMMPS executable.

## How to use the coupled LBM-DEM feature
The LBM-DEM coupling can be used by employing some simple LAMMPS commands (fix and pair_style) in your LAMMPS input script.

### fix lbm-psm
This fix employs the LBM method which updates the LBM distribution functions on each lattice node, computes the velocity field, and the hydrodynamic forces and torques acting on the particles as well as the stresslet for each particle.

fix ID group-ID lbm-psm every value Nlc value lc value rho value nu value Re value tau value Fext value_x value_y value_z

+ ID, group-ID are documented in fix command
+ lb-psm: Style name of this fix command
+ nevery value: Call the LBM methods to update the fluid fields and hydrodynamic interaction forces every this many timesteps
+ Nlc value: Number of grid nodes over the characteristic length
+ lc value: The characteristic lenght in the system (units of [length]). (E.g. the particle diameter)
+ rho value: The fluid density (units of [mass/length^3])
+ nu value: The kinematic viscosity of the fluid (units of [length^2/time]).
+ Re value: The Reynolds number defined as Uc*lc/nu. In the current implementation, the charactristic velocity Uc is computed and used to set the velocity on some boundaries (depending on the boundary conditions)
+ tau value: The LBM BGK (single) relaxation time parameter (optional, the default is tau = 0.65. Don't change it unless you know what you do)
+ Fext value_x value_y value_z: Sets a volume force, based on Guo et al. (2002), in the LBM method where each value_x/y/z defines the volume force component in x/y/z-direction, respectively. (Optional input. Units of mass/length^3*length/time^2, i.e. density*acceleration)

Example: fix lbDYN all lbm-psm every 1 Nlc 16 lc 0.0005 rho 1000.0 nu 0.0001 Re 0.1 tau 0.65 Fext -9810.0 0.0 0.0

### fix lbm-psm-vtk
Writes a vtk output file containing information on the fluid velocity field, fluid pressure field, and solid fractions of the lattice nodes.

fix ID group-ID lbm-psm-vtk every value

+ ID, group-ID are documented in fix command
+ lb-psm-vtk: Style name of this fix command
+ every value: Write the vtk output file every this many timesteps

Example: fix lbVTK all lbm-psm-vtk every 20000

### fix lbm-psm-bc
Sets a combination of boundary condition. If this fix is not used, all domain boundary conditions are periodic.\
If on one or more boundaries a velocity is set, then the magnitude is the characteristic velocity computed by the Reynolds number defined in the fix lbm-psm command (described above).\
At some point in the future, this command will be refactored so that a general user input is possible (i.e. defined pressure and/or velocity boundary conditions specifically for a domain boundary). 

fix ID group-ID lbm-psm-bc keyword

+ ID, group-ID are documented in fix command
+ lb-psm-bc: Style name of this fix command
+ keyword: Defines a combination of boundary conditions depending on the keyword which can be:
  + shear: Shear with shear gradient in y-direction. Boundaries in x-direction are periodic.
  + xFlow:  Fluid flow in positive x-direction by setting velocity on inlet boundary (boundary at smaller domain coordinate in x-direction)
      and constant density on outlet boundary (boundary at larger domain coordinate in x-direction).
      Lateral boundaries (in y-direction) have the imposed inlet velocity in x-direction.
  + channel: Two parallel plates in y-direction which are not moving. Can be used for channel flow driven by external force/pressure gradient.
  + channel-velIn-pressOut: Enclosed channel with velocity inlet, pressure outlet, and surrounding no-slip walls (can be used for fluidised bed simulation)
  + closedBox: Closed box with no-slip boundary conditions all around

Example: fix lbBC all lbm-psm-bc channel-velIn-pressOut

### fix lbRestart
Sets a combination of boundary condition. If this fix is not used, all domain boundary conditions are periodic.\
If on one or more boundaries a velocity is set, then the magnitude is the characteristic velocity computed by the Reynolds number defined in the fix lbm-psm command (described above).\
At some point in the future, this command will be refactored so that a general user input is possible (i.e. defined pressure and/or velocity boundary conditions specifically for a domain boundary). 

fix ID group-ID lbm-psm-restart every value write value read value

+ ID, group-ID are documented in fix command
+ lbm-psm-restart: Style name of this fix command
+ every value: Write the binary restart output file for the LBM distribution functions every this many timesteps
+ write value: Write the binart restart file when the value is set to 1 (does not write if set to 0)
+ read value: Read the written bindary restart file when the value is set to 1 (does not write if set to 0)

Example: fix lbRestart all lbm-psm-restart every 10000 write 1 read 0

### pair_style lubricate/GRM/LBDEM
Imposes a lubrication force and torque correction for particles in close proximity where the gap distance is below a defined threshold.

pair_style lubricate/GRM/LBDEM mu ci co componet calibrationFactor

+ style = lubricate/GRM/LBDEM
+ mu = The dynamic viscosity of the fluid (unit of [mass/(length*time)])
+ ci = Inner cut-off below which the correction force/torque correction is computed based on the inner cut-off value as gap distance
+ co = Outer cut-off below which the lubrication force/torque correction is computed (value should be set to the lattice cell width which is the theoretical resolution limit)
+ component: Integer value which defines which lubrication correction component is computed
  + 1: Squeezing terms X_A (leads to a normal force for a relative particle-particle motion along the center-to-center line)
  + 2: Shearing terms Y_A (leads to a tangential force due to a relative tangential particle-particle motion)
  + 3: Shearing terms Y_B (leads to a torque due to a relative tangential particle-particle motion)
  + 4: Terms Y_B for forces from rotations
  + 5: Rotational Y_C terms (leads to a torque due to particle rotations)
+ calibrationFactor: Modifies the outer cut-off by multiplicaiton with this factor for calibration purposes (more details listed in one of the below publications [Najuch and Sun (2023)])

Example: pair_style lubricate/GRM/LBDEM 0.1 5.0e-11 3.3333333333333335e-05 1 0.75


## References to cite if using this package
If you use this LAMMPS extension, then please cite one of the following papers and possibly this repository (which has a .cff file):

Paper analysing the underlying methodology. Could be cited in the methodology section of your paper:

+ Tim Najuch and Jin Sun, "Analysis of two partially-saturated-cell methods for lattice Boltzmann simulation of granular suspension rheology", Computers & Fluids, Volume 189, 15 July 2019, Pages 1-12, 
Url: https://www.sciencedirect.com/science/article/abs/pii/S0045793019301458

Paper describing the lubrication force and torque correction. Could be cited if you use the lubrication corrections:

+ Tim Najuch and Jin Sun, "Lubrication force correction and calibration for a partially-saturated-cell lattice Boltzmann method", Computers & Fluids (Submitted 2023)
