# Studio LAMMPS

Studio LAMMPS is an educational tool, built as an R Shiny application, that enables interactive control of LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator: http://lammps.sandia.gov/). By leveraging the power of LAMMPS as a proven research tool, a powerful interactive experience has been developed for students to explore key concepts in chemistry, physics, and materials science through molecular dynamics (MD) simulations without requiring an understanding of the underlying simulation script. With this framework, educators can develop LAMMPS simulations that are representative of important topics, convert them into the format used in this app, and add Shiny interactive inputs in a way that allows students to engage and explore the content while ensuring that the simulations remain feasible.

A secondary innovation, produced in the course of this project, is that LAMMPS has now been exposed to the R universe, introducing an entirely new computational environment for the MD package.


![screenshot of start page](/markdown/FullScreen.png)


### Table of Contents
**[Motivation, Objectives](#motivation)**<br>
**[Technology Stack](#technology-stack)**<br>
**[Software Design](#software-design)**<br>
**[Installation, Requirements, Shiny Server Setup](#installation)**<br>
**[Current Status](#current-status)**<br>
**[Future Work](#future-work)**<br>

## Motivation
This code was initially developed as a course project for MTLE-4500: Computational Materials Design, taught by Dr. Liping Huang at Rensselaer Polytechnic Insitute in Troy, NY. The purpose at the onset was to develop an educational software to introduce MD simulations and that could replace Java applets currently available online. The primary concern with present solutions is that few modern browsers allow Java (i.e. Chrome, Safari, MS Edge) and it can be challenging to configure Internet Explorer to run them. Additionally, since these applets were developed years ago and haven't been regularly maintained, there are questions as to whether these will be available in the future let alone whether students' computers will be compatible. 

With the knowledge of the Python wrapper to LAMMPS, it seemed possible to develop a means of interacting with and reading data out of the simulation package, which could be controlled with some graphical interface. While a number of applications exist to control LAMMPS locally, the specific intent was to develop a web-hosted application that students could access at a server address on the school's network. An important caveat with this difference is that while the traditional use of LAMMPS (sending a text input file to the executable) requires that data be exported to local files, all data read out of the simulation via the Python wrapper is stored in memory. The advantage of this design is that it avoids unwanted behaviors arising from multiple user sessions interacting with files in a shared location.  The following primary and additional (personally desired) objectives were defined:

### Primary Objectives
-	Dashboard style interface for MD simulations.
- Ability to change at least several input parameters.
-	Visualizations of the coordinate system, RDF, MSD, and thermodynamic properties.
-	Multiple simulation configurations available.
-	The whole stack should process simulations and return data in a reasonably short time.
### Additional Objectives
-	All desired data should be stored in memory and no data should be stored in local files.
-	The architecture should be designed in such a way that allows future developers to produce custom experiences.
-	Plots should be interactive and provide time series animations.

## Technology Stack
  - LAMMPS: http://lammps.sandia.gov/
  - Python and C-level wrappers for LAMMPS: https://lammps.sandia.gov/doc/Section_python.html
  - Python NumPy library: http://www.numpy.org/
  - R Reticulate package to run Python in R: https://rstudio.github.io/reticulate/
  - R Shiny package for UI and server-behaviors: https://shiny.rstudio.com/
  - R Plotly package for interactive plotting in the Shiny app: https://plot.ly/r/
  - R Shiny Server for hosting the app: https://www.rstudio.com/products/shiny/shiny-server/ 

## Software Design
At the onset of this project there were certain objectives
There are a number of peculiarities with this application

What may be unusual is that because this is an R Shiny app, it does not run on a desktop, but instead runs on a server. For those familiar with R and RStudio, it can be run out of RStudio locally for a single user and only temporarily. However, when the application is hosted on a Shiny Server with the necessary packages installed and shared libraries correctly placed, it can be accessed at a web address by anyone on the same network as the hosting server.

In terms of hosting the application on a Shiny Server, there are several considerations. When someone accesses the Shiny app, a separate R thread is created for each user. This may cause a number of issues depending upon the scale of the simulations possible, the available resources on the server computer, and the number of users that may connect.


What is unique to this kind of application is that it can be placed on a Shiny Server and with the correct 



There are several key technologies that are taken advantage of: 
  - LAMMPS as a highly optimized and validated MD engine
  - Python and C-level wrappers for LAMMPS
  - Reticulate R package to run Python in R
  - R Shiny UI framework
LAMMPS <-> Python <-> R <-> Shiny UI


## Installation
### Requirements
### Shiny Server Setup

## Customization
### Learning Hurdles
developing LAMMPS scripts
accessing LAMMPS data in Python depending upon type
reading this data into R

## Current Status

## Future Work 
  - Change Shiny app layout into Shinydashboard
  - Improve plot formatting, especially velocity histogram
  - Instead of enumerating options, create multiple types of simulations and provide unique interactivity upon these:
    * Crystal structure explorer
    * Tensile test
    * Diffusion
    * Phase change and spinodal decomposition
  - Parallel processing
  - Look into using pymatgen for new modules
