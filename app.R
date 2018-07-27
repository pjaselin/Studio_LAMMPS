## Required R library imports
library(shiny) # web app framework
library(plotly) # interactive plotting
library(reticulate) # Python -> R wrapper
library(gsubfn) # utilities for function arguments

# Load python helper functions
source_python("lammps_helper.py")

# Python library loads
np <- import("numpy")
lammps <- import("lammps")


accumulate_by <- function(dat, var) {
  # function to collect data frames by time frame in order
  # to provide cumulative plots wrt time
  # source: https://plot.ly/r/cumulative-animations/
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(0, x)],], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

#' @name  per_atom_prop
#' @description Function to read per atom data (x y z coordinates, velocity magnitude, atom IDs, and atom types) using the C-level LAMMPS-Python interface.
#' @param lmp C-type Python LAMMPS class
#' @param time_step Current simulation iteration
#' @param natoms Number of atoms
#' @param check_nstep Number of steps as interval for recording data
#'
#' @return Data frame with (time step, atom id, atom type, atom coordinates, velocities)
per_atom_prop <- function(lmp, time_step, natoms, check_nstep){
  # collect atom coordinates
  c_coords <- lmp$gather_atoms("x", 1L, 3L)
  coords <- np$reshape(np$array(c_coords),list(as.integer(natoms),3L))
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x","y","z")
  
  #collect velocties and find the magnitudes
  velocities <- lmp$gather_atoms("v",  1L, 3L)
  velocities <- np$reshape(np$array(velocities), list(as.integer(natoms), 3L))
  vel_mag <- np$sqrt(np$sum(velocities^2, axis = 1L))
  
  # collect atom IDs and types
  types <- np$array(lmp$gather_atoms("type",0L,1L))
  atom_id <- np$array(lmp$gather_atoms("id",0L,1L))
  
  # combine in one data frame
  coords_w_ids <- cbind.data.frame(time_step =rep(time_step*check_nstep, natoms), 
                                   atom_id = atom_id, types = types, 
                                   coords, velocities = vel_mag)
  # return single data frame
  return(coords_w_ids)
}

#' @name current_RDF
#' @description Function to extract the current RDF using the C-level LAMMPS-Python interface
#' @param lmp C-type LAMMPS class
#' @param time_step Current simulation iteration
#' @param check_nstep Number of steps as interval for recording data
#'
#' @return current RDF 
current_RDF <- function(lmp, time_step, check_nstep){
  rdfone <- lmp$extract_compute("myRDF",0L,2L)
  RDF <- as.data.frame(c_to_np(rdfone))
  colnames(RDF) <- c("Distances", "RDF")
  RDF <- cbind.data.frame(time_step =rep(time_step*check_nstep, dim(RDF)[1]),RDF)
  return(RDF)
}

#' @name current_thermo
#' @description Function to collect all thermodynamic data from thermo_style command
#' @param L PyLammps instance wrapper
#' @param time_step Current simulation iteration
#' @param check_nstep Number of steps as interval for recording data
#'
#' @return current thermodynamic data
current_thermo <- function(L,time_step,check_nstep){
  thermo_results <- data.frame(
    time_step = time_step*check_nstep,
    ke = tail(L$runs[[time_step+1]]$thermo$KinEng, n=1),
    pe = tail(L$runs[[time_step+1]]$thermo$PotEng, n=1),
    tote = tail(L$runs[[time_step+1]]$thermo$TotEng, n=1),
    temperature = tail(L$runs[[time_step+1]]$thermo$Temp, n=1),
    volume = tail(L$runs[[time_step+1]]$thermo$Volume, n=1),
    pressure = tail(L$runs[[time_step+1]]$thermo$Press, n=1),
    enthalpy = tail(L$runs[[time_step+1]]$thermo$Enthalpy, n=1),
    density = tail(L$runs[[time_step+1]]$thermo$Density, n=1)
  )
  return(thermo_results)
}

#' @name current_MSD
#' @description Function to return the current scalar MSD value from compute
#' @param lmp lammps wrapper class
#' @param time_step Current simulation iteration
#' @param check_nstep Number of steps as interval for recording data
#'
#' @return current MSD of system
current_MSD <- function(lmp, time_step, check_nstep){
  MSD_result <- data.frame(time_step = time_step*check_nstep, 
                           MSD=lmp$extract_compute("msdcal", 0L, 1L)$contents$value
                           )
  return(MSD_result)
}


#####################
# Studio LAMMPS App

# Define UI
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .shiny-text-output {
      background-color:#fff;
    }
  "))),
  
  # Application title
  #titlePanel("Studio LAMMPS"),
  h1(span("Studio", style = "font-weight: 300"), "LAMMPS", 
     style = "font-family: 'Source Sans Pro';
        color: #fff; text-align: center;
        background-image: url('bgtexturered.jpg');
        padding: 20px"),
  br(),
   
   
   # Sidebar
   sidebarLayout(
     position = "left",
     sidebarPanel(
       div(style = "display: inline-block;vertical-align:top; width: 150px;",
           radioButtons("dimensions", label = "Dimensions",
                    choices = list("2D" = 2, "3D" = 3), 
                    selected = 3, inline =TRUE)
           ),
       # dimensions options and lattice style
       div(style ="display: inline-block;vertical-align:top; width: 175px;",
           conditionalPanel(
             condition = "input.dimensions == 3",
             radioButtons("lattice3d", "Lattice Style",
                     choices = list("FCC"="fcc", "BCC"="bcc", "SC"="sc"), inline = TRUE, selected="fcc")),
           conditionalPanel(
             condition = "input.dimensions == 2",
             radioButtons("lattice2d", "Lattice Style",
                      choices = list("SQ"="sq", "SQ2"="sq2"), inline = TRUE, selected = 'sq')
       )),
       
       # lattice constance
       
       div(style = "display: inline-block;vertical-align:top; width: 150px;",
           sliderInput("lattice_constant", "Lattice Constant", 
                       min = 0.6, max = 3.2, value = 0.84, step = 0.01)
       ),
       # change number of atoms
       div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
       div(style ="display: inline-block;vertical-align:top; width: 125px;",
           selectInput('natoms', 'Number of Atoms', 
                         choices = list("4" = 1, "14" = 1.5, "32" = 2, "63" = 2.5, "108" = 3, 
                                        "172" = 3.5, "256" = 4, "365" = 4.5, "500" = 5), 
                       selected = 4, selectize = FALSE)),
       
       #boundary conditions
      tags$h5(tags$b("Boundary Conditions")),
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           selectInput("xbc", "X", choices=list("Periodic"="p", "Shrink-wrapped"="s"), selected = "p", selectize = FALSE)
       ),
       div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           selectInput("ybc", "Y",choices=list("Periodic"="p", "Shrink-wrapped"="s"), selected = "p", selectize = FALSE)
       ),
       div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           conditionalPanel(
             condition = "input.dimensions == 3",
             selectInput("zbc", "Z", choices=list("Periodic"="p", "Shrink-wrapped"="s"), selected = "p", selectize = FALSE)
       )),
       
       #mass
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           sliderInput("mass", "Mass", 
                       min = 0.5, max =10, value = 1.0, step = 0.1)
       ),
      #sigma
       div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           sliderInput("sigma", "σ", 
                       min = 0.5, max =10, value = 1.0, step = 0.1)
       ),
      #epsilon
       div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
       div(style = "display: inline-block;vertical-align:top; width: 100px;",
           sliderInput("epsilon", "ε", 
                       min = 0.5, max =10, value = 1.0, step = 0.1)
       ),
       # velocity
       sliderInput("vel", "Velocity",
                    min = 0, max = 25,
                    value = 1.44, step = 0.01),
      # ensemble choice
      selectInput("ensemble", "Ensemble",choices=list("NVE"="nve", "NVT"="nvt"), selected = "nve", selectize = FALSE),
      conditionalPanel(
        condition = "input.ensemble == 'nvt'",
        div(style = "display: inline-block;vertical-align:top; width: 100px;",
            sliderInput("nvtstarttemp", "Start Temp", 
                        min = 0.01, max =5, value = 0.1, step = 0.01)
        ),
        div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
        div(style = "display: inline-block;vertical-align:top; width: 100px;",
            sliderInput("nvtendtemp", "End Temp", 
                        min = 0.01, max =5, value = 0.1, step = 0.01)
        ),
        div(style="display: inline-block;vertical-align:top; width: 25px;", HTML("<br>")),
        div(style = "display: inline-block;vertical-align:top; width: 100px;",
            sliderInput("nvtq", "Q", 
                        min = 0, max =3, value = 0, step = 0.1)
        )
      ),
      
      #time step size
      sliderInput("timestep", "Time Step",
                  min = 0.00001, max = 0.01,
                  value = 0.005, step = 0.00001),
       # number of total simulation time
       sliderInput("total_steps", "Simulation Time",
                   min = 0, max = 100000,
                   value = 1000, step = 100),
      # run simulation button
       actionButton("runsim", "Submit to LAMMPS", 
                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      ),
      
      
      # main panel on the right layout
      mainPanel(
        # Output: Tabset w/ plots
        tabsetPanel(type = "tabs",
                    tabPanel("System", plotlyOutput("plot")),
                    tabPanel("RDF", plotlyOutput("RDFplot")),
                    tabPanel("MSD", plotlyOutput("MSDplot")),
                    tabPanel("Thermo", plotlyOutput("thermoplot")),
                    tabPanel("Velocity", plotlyOutput("velocityplot")),
                    tabPanel("V/P/H/ρ",
                             fluidRow(
                               column(6, plotlyOutput("Volume")),
                               column(6, plotlyOutput("Pressure"))
                             ),
                             fluidRow(
                               column(6, plotlyOutput("Enthalpy")),
                               column(6, plotlyOutput("Density"))
                             )
                             ),
                    tabPanel("Downloads",
                             br(),
                             fluidRow(downloadButton("downloadcoords", "Atom Property History")),
                             br(),
                             fluidRow(downloadButton("downloadRDF", "RDF History")),
                             br(),
                             fluidRow(downloadButton("downloadMSD", "MSD History")),
                             br(),
                             fluidRow(downloadButton("downloadthermo", "Thermo History"))
                             )
        )
      )
      
   )
)


# Define server logic
server <- function(input, output, session) {
  #update number of atom choices by lattice structure
  observe({
    if (input$dimensions == 3) {
      if (input$lattice3d == 'fcc') {
        updateSelectInput(session,'natoms', 'Number of Atoms', 
                           choices = list("4" = 1, "14" = 1.5, "32" = 2, "63" = 2.5, "108" = 3, 
                                          "172" = 3.5, "256" = 4, "365" = 4.5, "500" = 5),
                          selected=4)
      }
      if (input$lattice3d == 'bcc') {
        updateSelectInput(session,'natoms', 'Number of Atoms', 
                          choices = list("9" = 1.5, "16" = 2, "35" = 2.5, "54" = 3, 
                                         "91" = 3.5, "128" = 4, "189" = 4.5, 
                                         "250" = 5, "341"=5.5,"432" = 6, "559"=6.5, "686"=7), 
                          selected = 5)
      }
      if (input$lattice3d == 'sc') {
        updateSelectInput(session,'natoms', 'Number of Atoms', 
                          choices = list("8"=2, "27" = 3, "64"=4, "125"=5, 
                                         "216" = 6, "343"=7, "512"=8, "729"=9), 
                          selected = 8)
      }
    }
    
    if (input$dimensions == 2) {
      if (input$lattice2d == 'sq') {
        updateSelectInput(session,'natoms', 'Number of Atoms', 
                          choices = list("9" = 3, "16"=4, "36"=6, 
                                         "64" = 8, "100"=10, "144"=12, 
                                         "196"=14,"256"=16, "324"=18,"400"=20,"484"=22, "576"=24), 
                          selected = 20)
      }
      
      if (input$lattice2d == 'sq2') {
        updateSelectInput(session,'natoms', 'Number of Atoms', 
                          choices = list("8"=2, "32" = 4, "72"=6, "128"=8, 
                                         "200" = 10, "288"=12, "392"=14, "512"=16),selected = 14)
      }
    }
  })
  
  results <- eventReactive(input$runsim, {
    # Create a Progress object
    progress <- shiny::Progress$new(style = "notification")
    progress$set(message = "Running LAMMPS", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    updateProgress <- function(value = NULL, detail = NULL) {
      #if (is.null(value)) {
      #  value <- progress$getValue()
      #  value <- value + (progress$getMax() - value) / 5
      #}
      progress$set(value = value, detail = detail)
    }
    
    # LAMMPS input file
    
    lmp <- lammps$lammps()
    L <- lammps$PyLammps(ptr = lmp)
    
    if (input$dimensions == 2){
      L$dimension("2")
      L$boundary(paste(input$xbc, input$ybc, "p"))
    } else {
      L$boundary(paste(input$xbc, input$ybc, input$zbc))
    }
    
    L$units("lj")
    L$atom_style("atomic")
    L$atom_modify("map array") # required to read coordinates
    if (input$dimensions == 2){
      L$lattice(paste(input$lattice2d, input$lattice_constant))
    } else {
      L$lattice(paste(input$lattice3d, input$lattice_constant))
    }
    
    if (input$dimensions == 2){
      L$region(paste("box block 0", input$natoms, "0", input$natoms, "-0.25 0.25"))
    } else {
      L$region(paste("box block 0", input$natoms, "0", input$natoms, "0", input$natoms))
    }
    
    L$create_box("1 box")
    L$create_atoms("1 box")
    L$mass(paste("1", input$mass))
    
    L$velocity(paste("all create", input$vel, "87287 loop geom"))
    
    L$pair_style("lj/cut 2.5")
    L$pair_coeff(paste("1 1", input$sigma, input$epsilon,"2.5"))
    L$neighbor("0.3 bin")
    L$neigh_modify("delay 0 every 20 check no")
    
    if (input$ensemble =='nve'){
      L$fix("1 all nve")
    }
    if (input$ensemble == 'nvt'){
      L$fix(paste("1 all nvt temp", input$nvtstarttemp, input$nvtendtemp, input$nvtq))
    }
    
    if (input$dimensions == 2){
      L$fix("2 all enforce2d")
    }
    L$compute("msdcal all msd") #Compute msd
    L$compute("myRDF all rdf 50") # Compute rdf
    L$thermo_style("custom step temp ke pe etotal vol press enthalpy density")
    L$timestep(paste(input$timestep))
    L$run("0") # run for 0 steps to initialize compute results
    
    # Simulation
    number_atoms <- L$atoms$natoms
    check_nstep <- input$total_steps*0.05
    
    # Initial values lmp, time_step, natoms
    atom_history <- per_atom_prop(lmp, 0, number_atoms, check_nstep)
    RDF_history <- current_RDF(lmp, 0, check_nstep)
    MSD_history <- current_MSD(lmp, 0, check_nstep)
    thermo_history <- current_thermo(L,0, check_nstep)
    
    # Perform simulation
    for (step in 1:(input$total_steps/check_nstep)){
      # Run next n iterations
      L$run(paste(check_nstep))
      # current coordinates
      atom_history <- rbind.data.frame(
        atom_history,
        per_atom_prop(lmp, step, number_atoms, check_nstep))
      # current RDF
      RDF_history <- rbind.data.frame(
        RDF_history, 
        current_RDF(lmp, step, check_nstep))
      # current MSD
      MSD_history <- rbind.data.frame(
        MSD_history, 
        current_MSD(lmp, step, check_nstep))
      # current thermo
      thermo_history <- rbind.data.frame(
        thermo_history, 
        current_thermo(L,step, check_nstep))
      # update progress bar
      updateProgress(value =(step*check_nstep)/input$total_steps, 
                     detail=paste0("Time Step: ",step*check_nstep,"    Progress to Completion: ",round((step*check_nstep)/input$total_steps*100.0,digits=2),"%"))
    }
    
    # set up data frames for cumulative plots over time
    MSD_history <- MSD_history %>% accumulate_by(~time_step)
    thermo_history <- thermo_history %>% accumulate_by(~time_step)
    
    # shutdown LAMMPS
    L$close()
    lmp$close()
    
    # return computations for further use
    results <- list(atom_history, RDF_history, MSD_history, thermo_history)
  })
  
  # all plots
  
  output$plot <- renderPlotly({
    input$runsim
    dims <- isolate(input$dimensions)
    if (dims == 3){
      plot_ly(data = results()[[1]], x=~x, y=~y,z=~z,
              frame=~time_step, ids=~atom_id,type="scatter3d",
              mode="markers", 
              marker=list(size=16,line = list(color = 'rgba(152, 0, 0, .8)', width = 5)),
              showlegend=F)%>% 
        animation_opts(frame = 1000, transition=500, redraw=F, easing="linear")%>%
       layout(autosize = F, height = 450) %>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
    } else {
      plot_ly(data = results()[[1]], x=~x, y=~y,
              frame=~time_step, ids=~atom_id,type="scatter",
              mode="markers", 
              marker=list(size=16,line = list(color = 'rgba(152, 0, 0, .8)', width = 2)),
              showlegend=F)%>% animation_opts(frame = 1000, transition=500,redraw=F,easing="linear")%>%
        layout(autosize = F, height = 450) %>% 
        animation_slider(
          currentvalue = list(
            prefix = "Time Step "
          ))
    }
   })
   
   output$RDFplot <- renderPlotly({
     plot_ly(data = results()[[2]], x=~Distances, y=~RDF,
             frame=~time_step, ids=~Distances,type="scatter",
             mode="markers+lines",showlegend=F)%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "
         ))
   })
   
   output$MSDplot <- renderPlotly({
     plot_ly(data = results()[[3]], x=~time_step, y=~MSD,
             frame=~frame, type="scatter",
             mode="markers+lines",showlegend=F)%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "
         )) %>%
       layout(yaxis = list(type = "log"))
   })
   
   output$thermoplot <- renderPlotly({
     plot_ly(data = results()[[4]], x=~time_step,
             frame=~frame)%>% 
       add_trace(y = ~pe, name = 'PotEng',mode = 'markers+lines')%>%
       add_trace(y = ~ke, name = 'KinEng',mode = 'markers+lines')%>%
       add_trace(y = ~tote, name = 'TotEng',mode = 'markers+lines')%>%
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   })
   
   output$velocityplot <- renderPlotly({
     atomdata <- results()[[1]]
     velocity_data <- as.numeric(atomdata$velocities)
     toprange <- max(velocity_data)
     lowrange <- min(velocity_data)
     nbins <- 75
     
     all_timesteps <- atomdata$time_step
     unique_timesteps <- unique(all_timesteps)
     bin_seq <- seq(lowrange, toprange, length.out = nbins)
     
     timestepvector <- numeric()
     histcounts <- numeric()
     histmids <- numeric()
     for (i in 1:length(unique_timesteps)){
       select_indices <- which(all_timesteps %in% c(unique_timesteps[i]))
       current_hist <- hist(velocity_data[select_indices], breaks = bin_seq, plot = FALSE)
       histcounts <- c(histcounts, current_hist$counts)
       histmids <- c(histmids, current_hist$mids)
       timestepvector <- c(timestepvector, rep(unique_timesteps[i], length(current_hist$mids)))
     }
     
     histogramdata <- data.frame(time_steps = timestepvector, 
                                 hist_counts = histcounts, hist_mids = histmids)
     
     plot_ly(histogramdata, x=~hist_mids, y=~hist_counts,frame = ~time_steps, type ="scatter",
             mode = "markers+lines",fill = 'tozeroy',showlegend=F)%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   
   })
   
   output$Volume <- renderPlotly({
     plot_ly(data = results()[[4]], x=~time_step,y=~volume,
             frame=~frame, type="scatter",mode="markers+lines")%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   })
   
   output$Pressure <- renderPlotly({
     plot_ly(data = results()[[4]], x=~time_step,y=~pressure,
             frame=~frame, type="scatter",mode="markers+lines")%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   })
   
   output$Enthalpy <- renderPlotly({
     plot_ly(data = results()[[4]], x=~time_step,y=~enthalpy,
             frame=~frame, type="scatter",mode="markers+lines")%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   })
   
   output$Density <- renderPlotly({
     plot_ly(data = results()[[4]], x=~time_step,y=~density,
             frame=~frame, type="scatter",mode="markers+lines")%>% 
       animation_opts(frame = 1000,
                      transition=500,redraw=F,easing="linear")%>% 
       animation_slider(
         currentvalue = list(
           prefix = "Time Step "))
   })
   
   # downloads
   
   output$downloadcoords <- downloadHandler(
     filename = function() {
       paste("AtomPropertyHistory", ".csv", sep = "")
     },
     content = function(file) {
       write.csv(results()[[1]], file, row.names = FALSE)
     }
   )
   output$downloadRDF <- downloadHandler(
     filename = function() {
       paste("RDFHistory", ".csv", sep = "")
     },
     content = function(file) {
       write.csv(results()[[2]], file, row.names = FALSE)
     }
   )
   output$downloadMSD <- downloadHandler(
     filename = function() {
       paste("MSDHistory", ".csv", sep = "")
     },
     content = function(file) {
       write.csv(results()[[3]], file, row.names = FALSE)
     }
   )
   output$downloadthermo <- downloadHandler(
     filename = function() {
       paste("ThermoHistory", ".csv", sep = "")
     },
     content = function(file) {
       write.csv(results()[[4]], file, row.names = FALSE)
     }
   )
}

# Run the application 
shinyApp(ui = ui, server = server)
