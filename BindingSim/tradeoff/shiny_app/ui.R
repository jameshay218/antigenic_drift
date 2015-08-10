library(shiny)
library(ggplot2)
library(gridExtra)

numericInputRow<-function (inputId, label, value = "") 
{
    div(style="display:inline-block",
        tags$label(label, `for` = inputId), 
        tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(fluidPage(
    sidebarLayout(
        sidebarPanel(
                numericInputRow(inputId="p","p:",value=1),br(),
                numericInputRow(inputId="r","r:",value=1),br(),
                numericInputRow("b","b:",value=3),br(),
                numericInputRow("a","a:",value=0.7),br(),
                numericInputRow("c","c:",value=1),br(),
                numericInputRow("nv","nv:",value=2),br(),
                numericInputRow("q","q:",value=1),br(),
            numericInputRow("lifespan","Life span:",value=70),br(),
            numericInputRow("infectious_period","Infectious Period:",value=5),br(),
            numericInputRow("N_reinfect","Number of Reinfections:",value=20),br(),
            numericInputRow("D","D:",value=0.5),br(),
            numericInputRow("delta","delta:",value=0.2),br()
            ),
        
        mainPanel(
            plotOutput("Main")
            
            )
        
        )))
