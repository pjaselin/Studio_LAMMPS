# Studio LAMMPS

Studio LAMMPS is primarily an educational tool for introducing molecular dynamics (MD) simulations into the classroom without students requiring a knowledge of the simulation code. It is built as an R Shiny interactive dashboard that uses the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS: http://lammps.sandia.gov/) as a backend. In such a way, MD simulations running on LAMMPS can be introduced to the classroom without students actually needing to understand the simulation script.

There are several key technologies that are taken advantage of: LAMMPS as a highly
LAMMPS <-> Python <-> R <-> Shiny UI

## Objectives
There are several key objectives fulfilled by Studio LAMMPS:
  - Replaces Java applets that demonstrate rudimentary MD simulations (future-proofing and security advantages)
  - Total customization of LAMMPS simulations towards learning outcomes
  - Interactive inputs to LAMMPS simulations with measured results visualized

## Operation/Functionality

## Installation
### Requirements
### Shiny Server Setup


## Current Status

## Future Work 
  - Change Shiny app layout into Shinydashboard
  - Instead of enumerating options, create multiple types of simulations and provide unique interactivity upon these:
    * Crystal structure explorer
    * Tensile test
    * Diffusion
    * Phase change and spinodal decomposition

## References

Note this code was initially developed as a solo course project for MTLE-4500: Computational Materials Design at Rensselaer Polytechnic Insitute under Dr. Liping Huang.
