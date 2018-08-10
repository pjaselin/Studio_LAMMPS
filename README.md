# Studio LAMMPS

Studio LAMMPS is an educational tool, built as an R Shiny application, that enables interactive control of LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator: http://lammps.sandia.gov/). By leveraging LAMMPS as a proven research tool, a powerful interactive experience has been developed for students to explore key concepts in chemistry, physics, and materials science through molecular dynamics (MD) simulations without requiring an understanding of the underlying simulation script. With this framework, educators can develop LAMMPS simulations corresponding with learning modules, convert them into the format used in this app, and add Shiny interactive inputs in a way that allows students to engage and explore the content.

![screenshot of start page](/markdown/FullScreen.png)

### Table of Contents
**[Motivation, Objectives](#motivation)**<br>
**[Technologies Involved](#technologies-involved)**<br>
**[Software Design](#software-design)**<br>
**[How to Customize](#how-to-customize)**<br>
**[Installation, Requirements, Shiny Server Setup](#installation)**<br>
**[Current Status](#current-status)**<br>
**[Future Work](#future-work)**<br>

## Motivation
This code was first developed as a course project for MTLE-4500: Computational Materials Design, taught by Dr. Liping Huang at Rensselaer Polytechnic Insitute in Troy, NY. The purpose at the onset was to develop an educational software to introduce MD simulations and that could replace Java applets currently available online. The primary concern with present solutions is that few modern browsers allow Java (i.e. Chrome, Safari, MS Edge) and it can be challenging to configure Internet Explorer to run them. Additionally, since these applets were developed years ago and haven't been regularly maintained, there are questions as to whether these will be available in the future let alone whether students' computers will be compatible. 

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

## Technologies Involved
  - LAMMPS: http://lammps.sandia.gov/
  - Python wrapper class for Python and C-level control of LAMMPS: https://lammps.sandia.gov/doc/Section_python.html
  - Python NumPy library: http://www.numpy.org/
  - R Reticulate package to run Python in R: https://rstudio.github.io/reticulate/
  - R Shiny package for UI and server-behaviors: https://shiny.rstudio.com/
  - R Plotly package for interactive plotting: https://plot.ly/r/
  - R Shiny Server for app hosting: https://www.rstudio.com/products/shiny/shiny-server/ 

## Software Design
The application presently consists of three files:
- app.R: Main application file
- lammps_helper.py: Python function for reading RDF data with C-level LAMMPS-Python wrapper
- /www/bgtexturered.png: Background image for header

There are four layers in this software: the LAMMPS engine at the backend, the Python wrapper interfacing with LAMMPS, R running the Python wrapper via the Reticulate package, and R Shiny providing the graphical interface. This ordering is the direction in which the software was developed. Beginning with a functional LAMMPS script, the Python wrapper was fully explored with the ability to access data out of the LAMMPS instance. When it was found that it might be difficult to develop a Python web app, the R Reticulate package was used to convert the Python scripts to R. Finally the web app with its interactive control of the LAMMPS commands and visuals were created using R Shiny.

N.B.: It is assumed that developers have a working knowledge of LAMMPS, have studied the Python wrapper to LAMMPS, and are willing to learn the R Shiny framework. For the latter two, reading the available documentation and playing around with the code should be quite enlightening.

### Application Flow
1. User makes desired simulation settings.
2. User presses the "Submit to LAMMPS" button.
3. Desired settings are passed to placeholder variables in the LAMMPS script.
4. LAMMPS is launched and runs in steps of 1/20 of the total simulation time.
5. At every step, the helper functions are used to read desired data out of the simulation and store this in data frames in memory.
6. At the end of the simulation, LAMMPS is closed.
7. The data is passed back to the user in the form of visualization and is made available for download as a CSV.

### LAMMPS Scripts
Comparing how these scripts look:
![script comparison](/markdown/ScriptComparison.png)

The native LAMMPS script on the left is the simple 3D Lennard-Jones melting simulation as provided in the LAMMPS example files. 

In the Python script version, the first peculiarity to note is how a lammps() class (C-level LAMMPS-Python wrapper) is defined and then wrapped with a PyLammps() class (Python-level LAMMPS-Python wrapper). Since Python classes work with LAMMPS in different ways, they each store their own sets of data. By using them together, users can access nearly nearly all data produced during the simulation, while using the fastest means possible to read this information (its better to read coordinate, RDF, MSD data from C-level and simulation data at the Python-level). Also, note that in Python, the LAMMPS commands become Python class methods.

In the R script version and because of the Reticulate package, the new code looks much like the Python script, replacing the periods with the extract symbol ($). Unlike the Python script, LAMMPS commands cannot take a comma-separated values and instead require these lists to be contained within a paste function, passing a string to the command. Also note that R does not interpret integer values in the same way as Python. Thus, R interprets the value 1 as 1.0 and not as an integer. To ensure that a value is passed as an integer, the command as.integer(1) should be used or by appending the letter L (1L) after the number (as.integer() is preferred when using placeholder variables for interactivity and if integers are passed in the LAMMPS commands, they can written as strings (i.e."1")).

### Reading Data from LAMMPS
At the beginning of the app.R file as well as the lammps_helper.py file, functions are created to read data from LAMMPS using the Python wrappers. The methods implemented here are largely covered in the LAMMPS documentation. Note that the C-level Python wrapper (lmp) is used to quickly read atom coordinates, velocities, types, ids, MSD, and RDF, while the Python-level Python wrapper (L) is used to access time step # and all data stored by the LAMMPS thermo_style command. Also note that when reading data from the C-level as presently done, the data read are the current values whereas in using the Python-level wrapper as implemented, the last element is the new piece of data.

### Enabling Interactivity to LAMMPS
Interactivity is introduced when placeholder variables are inserted into the LAMMPS commands, thus creating a parameter that can be changed based upon a Shiny input. See the previous section about implementing the LAMMPS script in R for more notes on this and look at the implemented LAMMPS code where the server side of the app stores all simulation data in the "results" variable. From the app.R code, it can be seen that a high degree of interactivity can be introduced by adding if statements, choosing certain commands based upon inputs and thus changing the type of simulation.

### Unique Features
There are several unique features this application provides:
- In Plotly visualizations of 2D simulation boxes and at boundaries, atoms do not disappear and reappear at the matching boundary at the opposite side of the box. Instead, they slide across the screen, wrapping around to the corresponding boundary. This makes it easy to understand and track atoms moving across periodic boundaries.
- Since all data in the simulation is stored in memory, it can be downloaded to a local computer as a CSV for further analysis.
- Plotly allows for PNG's of graphs to be downloaded, making it easy for students to complete homework assignments and lab reports (less memory involved with PNG's of graphs rather than Excel graphs).

## Installation
### Requirements
- Install LAMMPS (GitHub method recommended): https://lammps.sandia.gov/download.html
- Enable optional packages as desired: https://lammps.sandia.gov/doc/Section_start.html
- Build LAMMPS as a shared library: https://lammps.sandia.gov/doc/Section_python.html#py-3
- Install Python wrapper as a system and personal library: https://lammps.sandia.gov/doc/Section_python.html#py-4
- Copy LAMMPS shared library (liblammps.so) from the lammps/src/ directory to the system /lib/ directory.

### Shiny Server Setup
- Install Shiny Server: https://www.rstudio.com/products/shiny/download-server/
- Place the files that comprised Studio LAMMPS in a folder using the same structure as here
- Place this repository under /srv/shiny-server/ (i.e. /srv/shiny-server/StudioLAMMPS/: app.r, lammps_helper.py, www/bgtexturered.png)

## How to Customize
### Learning Hurdles
Before embarking on developing new teaching modules, it is necessary to understand how native LAMMPS scripts translate into R Reticulate scripts along with a working knowledge of Shiny. Between these two sides of the program, understanding how to access desired data from LAMMPS to store in R data frames can be quite tricky and may require a great deal of sandboxing to nail down.

### Procedure for Creating New Modules
- Develop the native LAMMPS script and identify which values are to be manipulated
- Convert the LAMMPS script to Python and ensure that the script still works (some additions may be necessary i.e. the command "atom_modify map array" is necessary in order to read atom coordinates)
- Determine the appropriate Python wrapper means of accessing data from LAMMPS
- Convert the Python-LAMMPS script to R Reticulate and convert the Python wrapper communications to R (it may be necessary to write some of these in native Python and access those functions in R)
- Implement the R Reticulate LAMMPS scripts in the application along with the functions to read data from LAMMPS every few iterations.
- Add plots and other visualizations as necessary

## Current Status
Unfortunately there are several bugs with the ranges of parameters currently available, causing LAMMPS to crash. This is really the result of overextending the code and doing too much (opening a Pandora's box). The immediate next step will be changing the format into a Shinydashboard or similar layout to provide a nicer design and cleaner way to add more features/configurations.

## Future Work 
  - Change Shiny app layout into Shinydashboard (v2)
  - Improve plot formatting, especially velocity histogram
  - Instead of enumerating options, create multiple types of simulations and provide unique features for these:
    * Crystal structure explorer
    * Tensile test
    * Diffusion
    * Phase change
    * Spinodal decomposition
  - Parallel processing
  - Look into using pymatgen for new modules

## Support
If anyone is looking to implement this code on a server for education, I'm sure you'll run into issues. Hopefully these explanations should get you pretty far, but some troubleshooting is inevitable. It is my hope that this can develop into a mostly contained (yet still customizable) educational package to replace and improve upon the Java applets out there. Suggestions are also welcome!
