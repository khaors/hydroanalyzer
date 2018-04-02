#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("HydroAnalyzer-GUI: Shiny Interface (v0.1)"),

  #### Panel 'About' (right hand side)
  ##############################################################################
  sidebarPanel(
    imageOutput("uptc.logo", inline=TRUE),
    #
    p(HTML("<h5>This is HydroAnalyzer-GUI, the Shiny interface for analysis and
           evaluation of hydrological data in <strong>R</strong>.</h5>
           This application can be used for the Exploratory Data Analysis of
           hydrological variables (precipitation, discharge, temperature), consistency
           analysis, watershed analysis (watershed delineation, river network extraction),
           spatial analysis of hydrological variables (spatial correlation), water budget
           calculation (using direct, abdc model and long term approach), frequency
           analysis (return period estimation using maximum and minimum values), and
           hydrologic regionalization (using the regression and L-moments approaches).")),
    p(HTML('This package was developed by Oscar Garcia-Cabrejo, School of Geological
           Engineering, Universidad Pedagogica y Tecnologica de Colombia, Sogamoso,
           Boyaca, Colombia. Its source code is freely available on
           <a href="http://www.github.com/khaors/hydroanalizer">github</a>.')),
    br(),

    h3('References:'),
    p(HTML('<li> <span style="font-variant: small-caps;">V. T. Chow, Maidment, D. &
            Mays, L.</span> (1988).<I>Applied Hydrology</I>.
             McGraw-Hill Publishing Company; International Ed edition .</li>
           <li> <span style="font-variant: small-caps;">Maidment, D.</span>(1993).
           <I>Handbook of Hydrology</I>. McGraw-Hill Education. </li>
           <li> <span style="font-variant: small-caps;">Davie, T.</span> (2002).
           <I> Fundamentals of Hydrology</I> Routledge Fundamentals of Physical Geography,
           Routledge.</li>'))
    ),

  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      #########################################################################
      #                    Panel 'Import Data'
      #########################################################################
      tabPanel("Import Data",
               icon = icon("file"),
               h3("First step: import data"),

               p(HTML("To run the application, import your data set using the
                      import button below. Your data must be supplied in the form of a text/csv file.
                      If the importation is done properly, a preview of the data is displayed below.
                      When this is done, go to the next step: Exploratory Data Analysis.")),
               #
               br(),
               checkboxInput('header', ' Header?', TRUE),
               checkboxInput('rownames', ' Row names?', FALSE),
               selectInput('sep', 'Separator:',
                           c("Comma","Semicolon","Tab","Space"), 'Comma'),
               selectInput('quote', 'Quote:',
                           c("None","Double Quote","Single Quote"),
                           'Double Quote'),
               selectInput('dec', 'Decimal mark', c("Period", "Comma"),
                           'Period'),
               numericInput('nrow.preview','Number of rows in the preview:',20),
               numericInput('ncol.preview', 'Number of columns in the preview:',
                            10),
               fileInput('file1', 'Choose CSV/TXT File'),
               helpText("Note: The preview only shows a given number of observations, but
                         the analysis will consider the whole dataset."),
               tableOutput("view")
               ),
      #########################################################################
      #                    Panel 'Import Data'
      #########################################################################
      tabPanel("Exploratory Data Analysis",
               icon = icon("bar-chart-o"),
               h3("Second Step: Start to look at our data"),
               br(),
               p(HTML("In this step, a set of tools is used to gain insight into the data,
                      uncover the underlying structure, define important variables, detect
                      outliers and anomalies, test underlying assumptions, develop
                      parsimonious models.")),
               br(),
               h4("Variable"),
               uiOutput("EDAvarnames"),
               br(),
               h4("Histogram"),
               textInput(inputId = "EDAnbins", label = "Number Bins", value = "30"),
               radioButtons(inputId = "EDAloghist", label ='Scale',
                            choices = c("Arithmetic", "Log"), selected = "Arithmetic"),
               br(),
               h4("Autocorrelation Function"),
               textInput(inputId = "EDAmaxlag", label = "Maximum Lag",
                         value = '24'),
               br(),
               h4("Periodogram"),
               textInput(inputId = 'EDAfilter', label = 'Filter', value = "3,5"),
               textInput(inputId = 'EDAtaper', label = 'Taper', value = '0.1'),
               radioButtons(inputId = 'EDAlogspec', label = 'Scale',
                            choices = c("Arithmetic", "Log"), selected = "Arithmetic"),
               plotOutput("EDA.plot")
      ),
      #########################################################################
      #                    Panel 'Consistency'
      #########################################################################
      tabPanel("Consistency Analysis",
               icon = icon("newspaper-o"),
               h3("Consistency Analysis"),
               br(),
               h5("The tests included in this tab are used to determine if a
                  time series is homogeneous or not, or if two given time series
                  are consisten to each other. This type of analysis is helpful
                  in determining if corrections to the hydrological measurements
                  are requiered in the time series."),
               br(),
               selectInput(inputId = "consisttype", label = "Type",
                           choices = c(None= "None", Homogeneous="Homogeneous",
                                       Consistency = "Consistency"),
                           selected = "None"),
               br(),
               conditionalPanel(condition = 'input.consisttype == "Homogeneous"',
                                selectInput(inputId = "homogeneousmethod", 
                                            label = "Method",
                                            choices = c(None = "None", 
                                                        VonNeumannTest = "VonNeumannTest", 
                                                        CumulativeResiduals = "CumulativeResiduals"), 
                                            selected = "None")
                                ),
               br(),
               conditionalPanel(condition = 'input.consisttype == "Consistency"',
                                selectInput(inputId = "consistmethod", label = "Method",
                                            choices = c(None="None", DoubleMass="DoubleMass", Bois="Bois"),
                                            selected = "None")),
               uiOutput('consist1'),
               uiOutput('consist2'),
               uiOutput('consist3'),
               br(),
               plotOutput("consistency"),
               br(),
               tableOutput("homogeneity")
      ),
      #########################################################################
      #                    Panel 'Filling Missing Data'
      #########################################################################
      tabPanel("Filling Missing Data",
               icon = icon("paint-brush"),
               h3("Filling Missing Observations"),
               br(),
               h5("Sometimes the hydrologic time series are not complete due to
                  different reasons (equipment failure, extreme events, human
                  disturbances, mishandling of data records, accidental losses,
                  etc)"),
               br(),
               selectInput(inputId = "fillingtype", label = "Method",
                           choices = c(StationAverage="StationAverage",
                                       MonthAverage="MonthAverage",
                                       NormalRatio="NormalRatio",
                                       IDW = "IDW", Regression = "Regression")),
               br(),
               uiOutput("filling1"),
               uiOutput("filling2"),
               uiOutput("filling3")

      ),
      #########################################################################
      #                    Panel 'Watershed Analysis'
      #########################################################################
      tabPanel("Watershed Analysis",
               icon = icon("wrench")
      ),
      #########################################################################
      #                    Panel 'Spatial Analyisis'
      #########################################################################
      tabPanel("Spatial Analysis",
               icon = icon("map-o"),
               h3("Estimating the spatial distribution of hydrological variables"),
               br(),
               h5("In most cases the hydrological variables are sampled at specific
                   locations in space and using this information it is required to
                   know the values at unsampled locations. In this case geostatistics
                   comes to the rescue offering us a set of tools to solve this
                   challenging problem"),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   h4("Input Variables"),
                   br(),
                   fileInput('watershed.limit.fl', 'Choose a SHP File'),
                   fileInput('DEM.fl', 'Choose a DEM File'),
                   fileInput('rainfall.fl', 'Choose a CSV file')
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Spatial Correlation",
                              br()
                     ),
                     tabPanel("Hydrological Maps",
                       br(),
                       plotOutput(outputId = "hydrologic.maps")
                     )
                   )
                 )
               )

      ),
      #########################################################################
      #                    Panel 'Water Budget'
      #########################################################################
      tabPanel("Water Budget",
               icon = icon("money"),
               h3("Water Budget: How much water there is in a region during a period of
                  time"),
               br(),
               p(HTML("Using information of precipitation, temperature and discharge,
                       the different component of the hydrologic cycle are
                      determined.")),
               br(),
               selectInput(inputId = "budgetmethod", label = "Water Budget Method",
                           choices = c(None = "None",
                                       Direct="Direct",
                                       LongTerm = "LongTerm",
                                       ABCD="ABCD"), selected = "None"),
               br(),
               conditionalPanel(
                 condition = "input.budgetmethod == 'Direct'",
                 br()),
               uiOutput('budget1'),
               uiOutput('budget2'),
               uiOutput('budget3'),
               uiOutput('budget4'),
               uiOutput('budget5'),
               uiOutput('budget6'),
               uiOutput('budget7'),
               uiOutput('budget8'),
               plotOutput('water.budget'),
               br(),
               h4('Water Budget Results'),
               br(),
               #uiOutput('view.budget')
               dataTableOutput("view.budget")
      ),
      #########################################################################
      #                    Panel 'Base Flow Analysis'
      #########################################################################
      tabPanel("Base Flow Analysis",
               icon = icon("bath"),
               br(),
               h4("Discharge"),
               uiOutput("BFvarnames"),
               selectInput('time.base', "Time Base:", c('day','month','year')),
               selectInput("method", "Method:", c('None','Graphical','Nathan', 'Chapman', 'Eckhardt')),
               conditionalPanel(
                 condition = "input.method == 'Nathan'",
                 textInput("nathan.alpha", label = h5("alpha"), value = "0.925")),
               conditionalPanel(
                 condition = "input.method == 'Chapman'",
                 textInput("chapman.alpha", label = h5("alpha"), value = "0.925")),
               conditionalPanel(
                 condition = "input.method == 'Eckhardt'",
                 textInput("eckhardt.alpha", label = h5("alpha"), value = "0.925"),
                 textInput("eckhardt.bfi", label = h5('BFImax'), value = "0.8")),
               #sliderInput("plot.range","Time range= ", min = -1, max = 0, value = c(-.6,-.5)),
               plotOutput("baseflow")
      ),
      #########################################################################
      #                    Panel 'Frequency Analysis'
      #########################################################################
      tabPanel("Frequency Analysis",
               icon = icon("repeat"),
               h3("Frequency Analysis: Build a frequency model of your data"),
               br(),
               selectInput(inputId="freqselect", label = "Step",
                           choices = c(None = "None", SelectModel = "SelectModel",
                                       ParameterEstimation = "ParameterEstimation",
                                       ModelValidation = "ModelValidation",
                                       UncertaintyAnalysis = "UncertaintyAnalysis")),
               br(),
               uiOutput("freq1"),
               uiOutput("freq2"),
               uiOutput("freq3"),
               plotOutput("frequency"),
               uiOutput("parameter.estimates")
      ),
      #########################################################################
      #                    Panel 'Regionalization'
      #########################################################################
      tabPanel("Regionalization",
               icon = icon("globe")
      )

    )
)))
