library(shiny)
library(ggplot2)
library(gridExtra)


shinyUI(
    navbarPage("Antigenic Drift Simulation",
               #' Parameter exploration panel
               tabPanel("Parameters",
                        sidebarPanel(
                             h3(strong("Input Parameters")),
                            numericInput("p",label="p - change in binding avidity affecting probability of immune escape:",value=3),
                            numericInput("r",label="r - impact of previous exposure on immunity:",value=50),
                            numericInput("b",label="b - shape of V/survival relationship:",value=3),
                            numericInput("a",label="a - rate of change of V/survival relationship:",value=0.7),
                            numericInput("c",label="c - per day contact rate:",value=1),
                            numericInput("n",label="n - virus offspring per replication:",value=2),
                            numericInput("v",label="v - number of virions transmitted:",value=2),
                            numericInput("q",label="q - shift on x-axis:",value=1),
                            numericInput("N_reinfect",label="Number of Reinfections:",value=20),
                            numericInput("delta",label="Delta - distance between virus and host immunity:",value=0.2)
                        ),
                        
                        mainPanel(
                            plotOutput("Main")
                            
                        )
                        ),
               #' Simulation exploration panel
               tabPanel("Simulation",
                        sidebarPanel(
                            checkboxGroupInput("flags",
                                               label=h3(strong("Flags")),
                                               choices=list(
                                                   "Save SIR output"=1,
                                                   "Save virus characteristics"=2,
                                                   "Save virus pairwise distances"=3,
                                                   "Record run time"=4,
                                                   "Verbose"=5
                                                   ),
                                               selected=1),
                            h3(strong("Host Population Parameters")),
                            numericInput("dur","Duration",365),
                            numericInput("s0", "S0", 900000),
                            numericInput("i0", "I0", 100),
                            numericInput("r0", "R0", 50000),

                            numericInput("c", "Contact Rate", 1.5),
                            numericInput("mu", "Birth and Death Rate (ie. life expectancy in years)", 40),
                            numericInput("wane", "Waning Rate (waning immunity in days)", 25),
                            numericInput("gamma", "Recovery Time (in days)", 3),
                            numericInput("iniBinding", "Initial Binding Avidity ", 0.8),
                            checkboxGroupInput(
                                "scenarios",
                                h3(strong("Scenarios")),
                                choices=list(
                                    "Random drift; fixed binding avidity" = 1,
                                    "No drift; adaptive binding avidity; adaptive antigenic change" = 2,
                                    "Random drift; adaptive binding avidity; adaptive antigenic change" = 3,
                                    "Random drift; adaptive binding avidity; no adaptive antigenic change" = 4
                                    ),
                                selected=1
                            ),
                            h3(strong("Actions")),
                            actionButton("run","Run Simulation"),
                            actionButton("download","Download Results"),
                            actionButton("download_plots","Download Plots")                           
                            
                            
                            ),
                        mainPanel(
                            h3(strong("Scenario 1: Random drift; fixed binding avidity")),
                            plotOutput("sim_main_1"),
                            h3(strong("Scenario 2: No drift; adaptive binding avidity; adaptive antigenic change")),
                            plotOutput("sim_main_2"),
                            h3(strong("Scenario 3: Random drift; adaptive binding avidity; adaptive antigenic change")),
                            plotOutput("sim_main_3"),
                            h3(strong("Scenario 4: Random drift; adaptive binding avidity; no adaptive antigenic change")),
                            plotOutput("sim_main_4")
                        )
                        ),
               tabPanel("Phylogeny")
               )
)



