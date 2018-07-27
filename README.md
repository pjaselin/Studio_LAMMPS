# Studio LAMMPS

Studio LAMMPS is an R Shiny application that enables interactive control of the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS: http://lammps.sandia.gov/) for running molecular dynamics (MD) simulations. The primary interest is to leverage LAMMPS as a powerful research tool in order to develop an interactive experience for students to explore concepts through MD simulation without actually requiring an understanding of the underlying simulation script. This allows educators to develop LAMMPS simulations that are representative of key concepts in chemistry, physics, and materials science and present these in a way that allows students to engage and explore these topics.

What may be unusual is that because this is an R Shiny app, it does not run on a desktop, but instead runs on a server. For those familiar with R and RStudio, it can be run out of RStudio locally for a single user and only temporarily. However, when the application is hosted on a Shiny Server with the necessary packages installed and shared libraries correctly placed, it can be accessed at a web address by anyone on the same network as the hosting server.

In terms of hosting the application on a Shiny Server, there are several considerations. When someone accesses the Shiny app, a separate R thread is created for each user. This may cause a number of issues depending upon the scale of the simulations possible, the available resources on the server computer, and the number of users that may connect.


What is unique to this kind of application is that it can be placed on a Shiny Server and with the correct 



There are several key technologies that are taken advantage of: LAMMPS as a highly
LAMMPS <-> Python <-> R <-> Shiny UI

## Objectives
There are several key objectives fulfilled by Studio LAMMPS:
  - Replaces Java applets that demonstrate rudimentary MD simulations (future-proofing and security advantages)
  - Total customization of LAMMPS simulations towards learning outcomes
  - Application is completely standalone, stores nothing locally, all calculations and data are stored in memory
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
  - Parallel processing

## References

Note this code was initially developed as a solo course project for MTLE-4500: Computational Materials Design at Rensselaer Polytechnic Insitute under Dr. Liping Huang.
