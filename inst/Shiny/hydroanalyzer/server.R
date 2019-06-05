#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(moments)
library(lmom)
library(sp)
library(gstat)
library(rgdal)
library(DT)
library(hydroanalyzer)
#
# Define Global Variables
#
options(shiny.maxRequestSize=30*1024^2)
#
model.types <- c("None", "normal", "lognormal", "exponential", "gumbel", "weibull",
                  "gev", "pearson3", "logpearson3", "uniform")


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  ## Server variables
  server.env <- environment() # used to allocate in functions
  current.table <- NULL
  start.date <- NULL
  start.date.vec <- NULL
  current.ts <- NULL # this variable will contain the current time series
  current.sp <- NULL # This variable will contain the current spatial info
  current.varnames <- NULL
  cfreq <- NULL
  first <- TRUE
  water.budget.mt <- NULL
  watershed.limit.shp <- NULL
  watershed.dem.raster <- NULL
  watershed.rainfall.sp <- NULL
  # Output the uptc logo :
  output$uptc.logo <- renderImage(list(src="uptc_jpg.jpg"),
                                  deleteFile=FALSE)
  ## Panel 'Import data'
  dInput <- reactive({
    in.file <- input$file1
    if (is.null(in.file))
      return(NULL)

    the.sep <- switch(input$sep, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec, "Period"=".", "Comma"=",")
    if (input$rownames) {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, dec=the.dec)
    }
    if(!first){
      # Reset EDA tab
      varnames <- c("None", names(the.table))
      updateSelectInput(session, inputId = "EDAvarnames1", choices = varnames,
                        selected = "None")
      # Reset consistency tab
      updateSelectInput(session, inputId = "consisttype", selected = "None")
      # Reset water budget tab
      updateSelectInput(session, inputId = "Budget.prec", choices = varnames,
                        selected = "None")
      updateSelectInput(session, inputId = "Budget.evtcol", choices = varnames,
                        selected = "None")
      updateSelectInput(session, inputId = "budgetmethod", selected = "None")
      output$view.budget <- renderTable({ NULL })
      # Reset baseflow tab
      updateSelectInput(session, inputId = "BFvarnames1", choices = varnames,
                        selected = "None")
      updateSelectInput(session, inputId = "time.base", selected = "day")
      updateSelectInput(session, inputId = "method", selected = "None")
    }
    if(first)
      first <- FALSE
    server.env$first <- first
    # return the table
    nvar <- ncol(the.table)
    ndat <- nrow(the.table)
    server.env$current.table <- the.table
    server.env$current.varnames <- names(the.table)
    print(the.table)
    the.table
  })
  # data preview table
  output$view <- renderTable({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
    if (ncol(d.input)>input$ncol.preview)
      d.input <- d.input[,1:input$ncol.preview]
    head(d.input, n=input$nrow.preview)
  })
  #######################################################################################
  #                                  EDA Tab
  #######################################################################################
  output$EDAvarnames <- renderUI({
    if(!is.null(server.env$current.table)){
      current.names <-  c("None", server.env$current.varnames)
      selectInput(inputId= "EDAvarnames1", label = "Current Variable",
                  choices = current.names)
    }
  })
  #
  output$EDA.plot <- renderPlot({
    current.table <- server.env$current.table
    current.var <- input$EDAvarnames1
    if(is.null(current.table))
      return(NULL)
    if(is.null(current.var) || current.var == "None")
      return(NULL)
    # Read input parameters
    nd <- nrow(current.table)
    nbins <- as.numeric(input$EDAnbins)
    log.hist <- input$EDAloghist
    max.lag <- as.numeric(input$EDAmaxlag)
    if(max.lag >= nd){
      max.lag <- round(nd/2)
    }
    span.filter <-  as.numeric(unlist(strsplit(input$EDAfilter,",")))
    v.taper <- as.numeric(input$EDAtaper)
    log.spec <- input$EDAlogspec
    tt <- seq(1, nd, 1)
    V <- as.matrix(unname(current.table[current.var]))
    data.df <- data.frame(t = tt, V = V)
    current.ts <- ts(V)
    fn <- frequency(current.ts)
    # Calculate the autocorrelation function
    data.acf <- acf(current.ts, lag.max = max.lag, plot = F)
    data.acf.df <- data.frame(lag = data.acf$lag, ymax = data.acf$acf,
                                ymin = matrix(0,max.lag+1))
    # Calculate the spectral density
    current_spectrum <- spec.pgram(current.ts, spans = span.filter, taper = v.taper,
                                   detrend = T, plot = F)
    spectrum.df <- data.frame(freq = current_spectrum$freq/fn, spec = current_spectrum$spec)
    #
    # Create time series plot
    #
    tsbase <- ggplot(data = data.df) + geom_line(aes(x = t, y = V), color = "#0066CC") +
              xlab("Date") + ylab("Variable") + ggtitle("a) Time Series")
    #
    # Create histogram
    #
    pbaseh <- ggplot(data = data.df, aes(x = V)) + geom_histogram(aes(y = ..density..),
                                                                colour="black",
                                                                fill="white",
                                                                bins = nbins) +
      geom_density(alpha=0.2, fill = "#FF6666") +
      geom_vline(aes(xintercept=mean(V, na.rm=T)),   # Ignore NA values for mean
                 color = "red", linetype = "dashed", size = 1)+
      xlab("Variable") + ggtitle("b) PDF")

    if(log.hist == "Log"){
      pbaseh <- pbaseh + scale_x_log10()
    }
    #
    # Autocorrelation function
    #
    pacf <- ggplot() + geom_linerange(aes(x = lag, ymin = ymin, ymax = ymax),
                                      data = data.acf.df,
                                      color = "#6666CC" ) +
      xlab("Lag") + ylab("Correlation Coeff.") + ggtitle("c) Autocorrelation Function")
    #
    # Periodogram function
    #
    pspec <- ggplot() + geom_line(aes(x = freq, y = spec), data = spectrum.df, color = "#6666CC") +
      xlab("Frequency") + ylab("Spectral Energy") + ggtitle("d) Smoothed Periodogram")
    if(log.spec == "Log"){
      pspec <- pspec + scale_y_log10()
    }
    layout.matrix <- rbind(c(1, 2, 5), c(3, 4, 5))
    grid.arrange(tsbase, pbaseh, pacf, pspec,
                 nrow=2, ncol = 2,
                 as.table=TRUE,
                 heights=c(1,1))

    val <- c(format(mean(data.df$V), digits = 5),
             format(sqrt(var(data.df$V)), digits = 5),
             format(skewness(data.df$V), digits = 5),
             format(kurtosis(data.df$V), digits = 5),
             format(mean(data.df$V)/sqrt(var(data.df$V)), digits = 5),
             length(data.df$V),
             format(max(data.df$V), digits = 5),
             format(quantile(data.df$V, 0.75), digits = 5),
             format(quantile(data.df$V,.5), digits = 5),
             format(quantile(data.df$V,.25), digits = 5),
             format(min(data.df$V), digits = 5))
    val <- matrix(val, 11)
    #SummaryTable <- data.frame( mean = )
    tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
    tbl <- tableGrob(val,
                     rows = c("Mean", "Std.Deviation", "Skewness", "Kurtosis",
                              "Coef.variation", "Number Data", "Maximum",
                              "Q3", "Q2", "Q1", "Minimum"),
                     cols = 'Value',
                     theme=tt)
    #
    grid.arrange(tsbase, pbaseh, pacf, pspec, tbl,
                 layout_matrix = layout.matrix, widths = c(1, 1, 0.5))
  })
  ##########################################################################################
  #                                Consistency Tab
  ##########################################################################################
  # #
  output$consist1 <- renderUI({
    tmp <- NULL
    if( is.null(input$consisttype) ||  input$consisttype == "None")
      return(NULL)
    #if(input$consistmethod == "None")
    #  return(NULL)
    #print("InsideUI1")
    #print(input$consisttype)
    if(input$consisttype == "Homogeneous"){
      if(input$homogeneousmethod == "VonNeumannTest"){
        if(!is.null(server.env$current.table)){
          current.names <- c("None", server.env$current.varnames)
          tmp <- selectInput(inputId = "vonnewmann.ref", label = "Test Station", 
                             choices = current.names)
        }
      }
      else if(input$homogeneousmethod == "CumulativeResiduals"){
        if(!is.null(server.env$current.table)){
          current.names <- c("None", server.env$current.varnames)
          tmp <- selectInput(inputId = "cumresiduals.ref", label = "Test Station", 
                             choices = current.names)
        }
      }
    }
    if(input$consisttype == "Consistency"){
      if(input$consistmethod == "DoubleMass"){
        if(!is.null(server.env$current.table)){
          current.names <-  c("None", server.env$current.varnames)
          tmp <- selectInput(inputId= "doublemass.ref", label = "Reference Station",
                             choices = current.names)
        }
      }
      else if(input$consistmethod == "Bois"){
        if(!is.null(server.env$current.table)){
          current.names <-  c("None", server.env$current.varnames)
          tmp <- selectInput(inputId= "bois.ref", label = "Reference Station",
                             choices = current.names)
        }
      }
    }
    return(tmp)
  })
  #
  output$consist2 <- renderUI({
    tmp <- NULL
    if(input$consisttype == "None")
      return(NULL)
    #if(input$consistmethod == "None")
    #  return(NULL)
    if(input$consisttype == "Consistency"){
      if(input$consistmethod == "DoubleMass"){
        if(!is.null(server.env$current.table)){
          current.names <-  c("None", server.env$current.varnames)
          tmp <- selectInput(inputId= "doublemass.test", label = "Test Station",
                             choices = current.names)
        }
      }
      else if(input$consistmethod == "Bois"){
        if(!is.null(server.env$current.table)){
          current.names <-  c("None", "Null", server.env$current.varnames)
          tmp <- selectInput(inputId= "bois.test", label = "Test Station",
                             choices = current.names)
        }
      }
    }
    return(tmp)
  })
  #
  output$consist3 <- renderUI({
    tmp <- NULL
    if(input$consisttype == "None")
      return(NULL)
    if(input$consistmethod == "None")
      return(NULL)
    if(input$consisttype == "Consistency"){
      if(input$consistmethod == "Bois"){
        tmp <- textInput(inputId = "bois.alpha", label = "Significance Level alpha",
                         value = "0.025")
      }
    }
    return(tmp)
  })
  #
  output$consistency <- renderPlot({
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    if(input$consisttype == "None")
      return(NULL)
    if(input$consistmethod == "None"){
      return(NULL)
    }
    pres <- NULL
    if(input$consisttype == "Consistency"){
      if(input$consistmethod == "DoubleMass"){
        ref.var <- input$doublemass.ref
        if(is.null(ref.var) || ref.var == "None")
          return(NULL)
        test.var <- input$doublemass.test
        if(is.null(test.var) || test.var == "None")
          return(NULL)
        ref <- as.matrix(unname(current.table[ref.var]))
        test <- as.matrix(unname(current.table[test.var]))
        doublemass <- double_mass_curve(ref,test)
        dm.df <- data.frame(ref = doublemass$S1ref, test = doublemass$S1)
        #
        # Create plot
        #
        mx.val <- max(c(max(doublemass$S1ref),max(doublemass$S1)))
        line.df <- data.frame(x=c(0,mx.val),y=c(0,mx.val))
        pres <- ggplot() + geom_line(aes(x = ref, y = test), data = dm.df) +
          geom_line(aes(x = x, y = y), data =line.df, col ="red") +
          coord_fixed() + xlab(paste0("Reference Station:",ref.var)) +
          ylab(paste0("Test Station:",test.var))
        #
      }
      else if(input$consistmethod == "Bois"){
        ref.var <- input$bois.ref
        if(is.null(ref.var) ||
           ref.var == "None")
          return(NULL)
        test.var <- input$bois.test
        if(is.null(test.var) || test.var == "None")
          return(NULL)
        ref <- as.matrix(unname(current.table[ref.var]))
        ndat <- nrow(ref)
        alpha <- as.numeric(input$bois.alpha)
        bois.results <- NULL
        if(test.var != "NUll"){
          test <- as.matrix(unname(current.table[test.var]))
          bois.results <- bois_test(Serie1 = ref, Serie2 = test, alpha = alpha)
        }
        else{
          bois.results <- bois_test(Serie1 = ref, alpha = alpha)
        }
        bois.res.df <- data.frame(tt = 0:ndat, residuals = bois.results$residuals,
                                  ellipse.inf = bois.results$ellipse.inf,
                                  ellipse.sup = bois.results$ellipse.sup)
        #
        # Check residuals outside the ellipse
        #
        pos_pos <- bois.results$residuals > bois.results$ellipse.sup
        pos_neg <- bois.results$residuals < bois.results$ellipse.inf
        #
        # Create plot
        #
        pres <- ggplot() + geom_line(aes(x = tt, y = residuals), data = bois.res.df, col = "red") +
          geom_line(aes(x = tt, y = ellipse.inf), data = bois.res.df) +
          geom_line(aes(x = tt, y = ellipse.sup), data = bois.res.df) +
          xlab("Time") + ylab("Cumulative Residual")
        if(test.var != "Null"){
          pres <- pres + ggtitle(paste0("Results of Bois Test: Stations ", ref.var,"(Ref.), ",
                                        test.var,"(Test)"))
        }
        else {
          pres <- pres + ggtitle(paste0("Results of Bois Test: Stations ", ref.var,"(Ref.)"))
        }
        if(!is.null(pos_pos) & length(pos_pos) > 0){
          tmp <- cbind(bois.res.df$ellipse.sup, bois.res.df$residuals)
          bois.res.df$area.sup <- apply(tmp, 1, max)
          pres <- pres + geom_ribbon(aes(x = tt, ymin = ellipse.sup, ymax = area.sup),
                                     data = bois.res.df, color = "black", alpha = 0.6, fill = "red")
        }
        if(!is.null(pos_neg) & length(pos_neg) > 0){
          tmp <- cbind(bois.res.df$ellipse.inf, bois.res.df$residuals)
          bois.res.df$area.inf <- apply(tmp, 1, min)
          pres <- pres + geom_ribbon(aes(x = tt, ymin = ellipse.inf, ymax = area.inf),
                                     data = bois.res.df, color = "black", alpha = 0.6, fill = "red")
        }
      }
    }
    return(pres)
  })
  #
  output$homogeneity <- renderTable({
    #print("inside")
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    if(input$consisttype == "None")
      return(NULL)
    #if(input$consistmethod == "None"){
    #  return(NULL)
    #}
    tab <- NULL
    if(input$consisttype == "Homogeneous"){
      if(input$homogeneousmethod == "VonNeumannTest"){
        test.var <- input$vonnewmann.ref
        if(is.null(test.var) || test.var == "None")
          return(NULL)
        test <- as.matrix(unname(current.table[test.var]))
        res <- von_neumann_test(test)
        print(res)
        tab <- as.data.frame(as.matrix(summary(test)))
        names(tab) <- test.var
        print(tab)
      }
      else if(input$homogeneousmethod == "CumulativeResiduals"){
        test.var <- input$cumresiduals.ref
        if(is.null(test.var) || test.var == "None")
          return(NULL)
        test <- as.matrix(unname(current.table[test.var]))
        res <- cumulative_deviation_test(test, 0.025) 
      }
    }
    return(tab)
  })
  ########################################################################################
  #                             Spatial Analysis TAB
  ########################################################################################
  read_watershed_limit <- reactive({
    in.file <- input$watershed.limit.fl
    if (is.null(in.file))
      return(NULL)
    watershed.limit.shp <- read.OGR(dsn=in.file, layer)
    server.env$watershed.limit.shp <- watershed.limit.shp
  })
  #
  read_dem_file <- reactive({
    in.file <- input$DEM.fl
    if(is.null(in.file))
      return(NULL)
    watershed.dem.raster <- readGDAL(in.file)
    server.env$watershed.dem.raster <- watershed.dem.raster
  })
  #
  read_rainfall_file <- reactive({
    in.file <- input$rainfall.fl
    if(is.null(in.file))
      return(NULL)
    rainfall
  })
  ##########################################################################################
  #                                Water Budget Tab
  ##########################################################################################
  output$budget1 <- renderUI({
    if(input$budgetmethod == "Direct" || input$budgetmethod == "ABCD"){
      if(!is.null(server.env$current.table)){
        current.names <-  c("None", server.env$current.varnames)
        selectInput(inputId= "Budget.prec", label = "Precipitation",
                    choices = current.names)}
    }
  })
  #
  output$budget2 <- renderUI({
    if(input$budgetmethod == "Direct" || input$budgetmethod == "ABCD"){
      if(!is.null(server.env$current.table)){
        radioButtons(inputId = "EVTcalculation", label = "EVT Calculation",
                     choices = c(EVT = "EVT", Temperature = "Temperature"))
        }
    }
  })
  #
  output$budget3 <- renderUI({
    if(input$budgetmethod == "Direct" || input$budgetmethod == "ABCD"){
      if(!is.null(server.env$current.table)){
        current.names <-  c("None", server.env$current.varnames)
        selectInput(inputId= "Budget.evtcol", label = "Variable",
                    choices = current.names)}
    }
  })
  #
  output$budget4 <- renderUI({
    if(input$budgetmethod == "Direct"){
      if(!is.null(server.env$current.table)){
        current.names <-  c("None", server.env$current.varnames)
        textInput(inputId= "Budget.rmax", label = "Max. Soil Retention",
                    value = "100")}
    }
    else if(input$budgetmethod == "ABCD"){
      selectInput(inputId = "ABCD.scale", label = "Temporal Scale ", 
                choices = c("None", "Month", "Year"), selected = "None")
    }
  })
  #
  output$budget5 <- renderUI({
    if(input$budgetmethod == "Direct"){
      if(!is.null(server.env$current.table)){
        checkboxInput(inputId = "Budget.table", label = "Water Budget Table",
                      value = FALSE)
      }
    }
    else if(input$budgetmethod == "ABCD"){
      if(!is.null(server.env$current.table)){
        current.names <-  c("None", server.env$current.varnames)
        selectInput(inputId = "Budget.Q", label = "Discharge", 
                    choices = current.names, selected = "None")
      }
    }
  })
  #
  output$budget6 <- renderUI({
    tmp <- NULL
    if(input$budgetmethod == "Direct" || input$budgetmethod == "ABCD"){
      if(!is.null(server.env$current.table)){
        tmp <- actionButton(inputId = "run.budget", label = "Calculate Water Budget")
      }
    }
    return(tmp)
  })
  #
  calculate_water_budget_direct <- function(){
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    var.prec <- isolate(input$Budget.prec)
    if(is.null(var.prec) || var.prec == "None")
      return(NULL)
    var.evt <- isolate(input$Budget.evtcol)
    if(is.null(var.evt) || var.evt == "None")
      return(NULL)
    EVT <- NULL
    Prec <- as.matrix(current.table[var.prec])
    print(Prec)
    if(var.evt == "EVT"){
      EVT <- as.matrix(current.table[var.evt])
    }
    else{
      EVT <- evt_thornthwaite(as.matrix(current.table[var.evt]))
    }
    Rmax <- isolate(as.numeric(input$Budget.rmax))
    water.budget <- NULL
    budget.method <- isolate(input$budgetmethod)
    if(budget.method == "Direct"){
      water.budget <- water_budget_direct(Prec, EVT, Rmax)
    }
    if(is.null(water.budget))
      return(NULL)
    precip <- cbind(rep("Precipitation", 12))
    evtt <- cbind(rep("PotentialEVT", 12))
    storage <- cbind(rep("Storage",12))
    deficit <- cbind(rep("Deficit",12))
    excess <- cbind(rep("Excess",12))
    runoff <- cbind(rep("Runoff",12))
    perc <- cbind(rep("Recharge",12))
    evtr <- cbind(rep("RealEVT",12))
    process <- rbind(storage, deficit, excess, runoff, perc, evtr)
    quantities <- rbind(cbind(water.budget$R), cbind(water.budget$WD),
                        cbind(water.budget$WE), cbind(water.budget$RN),
                        cbind(water.budget$PC), cbind(water.budget$ETR))
    months <- cbind(rep(seq(1,12,1),6))
    water.budget.df <- data.frame(month = months, quantities = quantities,
                                  process = process)
    process1 <- rbind(cbind(rep("Precipitation", 12)), cbind(rep("EVT", 12)))
    data.df <- data.frame(month = months[1:24], quantities = rbind(unname(Prec), 
                                                                   unname(EVT)),
                          process = process1)
    names(data.df) <- c("month", "quantities", "process")
    #
    # Create plots
    #
    vmax <- max(c(max(Prec),max(EVT)))
    phydro <- ggplot() + geom_line(aes(x = month, y = quantities, color = process),
                                   data = data.df) +
      ylim(0, vmax) +
      theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "Storage"
    preserve <- ggplot() +  geom_area(aes(x = month, y = quantities, fill = process),
                                      data = water.budget.df[pos,]) +
      ylim(0, vmax) + xlim(1,12) +
      scale_fill_manual(values=c("#3399FF"))  + theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "RealEVT"
    pevt <- ggplot() +  geom_area(aes(x = month, y = quantities, fill = process),
                                  data = water.budget.df[pos,]) + ylim(0, vmax) +
      scale_fill_manual(values=c("#00CC33")) + theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "Deficit"
    pdeficit <- ggplot() +  geom_area(aes(x = month, y = quantities, fill = process),
                                      data = water.budget.df[pos,]) + ylim(0, vmax) +
      scale_fill_manual(values=c("#CC6666")) + theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "Excess"
    pexcess <- ggplot() +  geom_area(aes(x = month, y = quantities, fill = process),
                                     data = water.budget.df[pos,]) + ylim(0, vmax) +
      scale_fill_manual(values=c("#3366FF")) + theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "Runoff"
    prunoff <- ggplot() +  geom_area(aes(x = month, y = quantities, fill = process),
                                     data = water.budget.df[pos,]) + ylim(0, vmax) +
      scale_fill_manual(values=c("#FF0033")) + theme(legend.title=element_blank())
    #
    pos <- water.budget.df$process == "Recharge"
    precharge <- ggplot() + geom_area(aes(x = month, y = quantities, fill = process),
                                      data = water.budget.df[pos,]) + 
      ylim(-vmax, vmax) +
      scale_fill_manual(values=c("#3366FF")) + 
      theme(legend.title=element_blank())
    #
    # Create dataframe for Table
    #
    if(input$Budget.table){
      water.budget.mt <- matrix(0.0, nrow = 12, ncol = 8)
      water.budget.mt[,1] <- Prec
      water.budget.mt[,2] <- EVT
      water.budget.mt[,3] <- as.numeric(water.budget$R)
      water.budget.mt[,4] <- as.numeric(water.budget$ETR)
      water.budget.mt[,5] <- as.numeric(water.budget$WE)
      water.budget.mt[,6] <- as.numeric(water.budget$WD)
      water.budget.mt[,7] <- as.numeric(water.budget$RN)
      water.budget.mt[,8] <- as.numeric(water.budget$PC)
      #
      water.budget.mt <- t(water.budget.mt)
      rownames(water.budget.mt) <- c("Prec","EVT","Reserve","ETR","Excess","Deficit","Runoff","Recharge")
      colnames(water.budget.mt) <- c("Jan","Feb","March","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      server.env$water.budget.mt <- water.budget.mt
      #output$view.budget <- renderTable({ water.budget.mt }, rownames = TRUE)
      output$view.budget <- renderDataTable({
        datatable(water.budget.mt,
                  extensions = c('Buttons'),
                  options = list(
                    pageLength = 10,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                    text = 'Download',
                    scrollY = 200)) %>% formatRound(1:12, 3)
      })
    }
    #
    res <- grid.arrange(phydro, preserve, pevt, pdeficit, precharge, prunoff,
                        nrow = 3, ncol = 2,
                        as.table=TRUE,
                        heights=c(3,3,3))
    return(res)
  }
  #
  calculate_water_budget_abcd <- function(){
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    var.prec <- isolate(input$Budget.prec)
    if(is.null(var.prec) || var.prec == "None")
      return(NULL)
    var.evt <- isolate(input$Budget.evtcol)
    if(is.null(var.evt) || var.evt == "None")
      return(NULL)
    EVT <- NULL
    Prec <- as.matrix(current.table[var.prec])
    if(var.evt == "EVT"){
      EVT <- as.matrix(current.table[var.evt])
    }
    else{
      EVT <- evt_thornthwaite(as.matrix(current.table[var.evt]))
    }
    var.Q <- isolate(input$Budget.Q)
    if(is.null(var.Q) || var.Q == "None")
      return(NULL)
    Q <- as.matrix(current.table[var.Q])
    # Calibration parameters
    obj.fn <- isolate(input$objfn) #c("rss", "mnad", "mxad")
    opt.method <- isolate(input$optmethod) # c("lbfgs", "ga", "sa")
    # Search region
    if(input$budgetmethod == "ABCD"){
      scale <- isolate(input$ABCD.scale)
      if(is.null(scale) || scale == "None")
        return(NULL)
      if(scale == "Month"){
        args <- list(Prec = P, PET = EVT)
        fn.model <- "abcd.month.model"
        res <- calibrate(Q, fn.model, obj.fn = obj.fn,
                         opt.method = opt.method, args)
      }
      else if(scale == "Year"){
        args <- list(Prec = P, PET = EVT)
        fn.model <- "abcd.year.model"
        res <- calibrate(Q, fn.model, obj.fn = objn.fn,
                         opt.method = opt.method, args)
      }
    }

    return(res)
  }
  #
  observeEvent(input$run.budget,{
    print("run.budget")
    output$water.budget <- renderPlot({
      res.plot <- NULL
      budget.method <- isolate(input$budgetmethod)
      if(budget.method == "Direct"){
        res.plot <- calculate_water_budget_direct()
      }
      else if(budget.method == "ABCD"){
        res.plot <- calculate_water_budget_abcd()
      }
      return(res.plot)
    })
  })

  # observeEvent(input$transformationRun, {
  #   output$transform_results <- renderDataTable({
  #     current.ves <- server.env$current.ves
  #     if(input$transform_results_plot){
  #       validate(
  #         need(!is.null(current.ves), "The VES is not defined")
  #       )
  #       #
  #       if(is.null(current.ves))
  #         return(NULL)
  #       #
  #       current.transformation <- isolate(input$transformation.type)
  #       validate(
  #         need(current.transformation != "None", "Please choose a Transformation")
  #       )
  #       current.transf <- NULL
  #       if(current.transformation == "Direct"){
  #         current.transf <- "direct"
  #       }
  #       else if(current.transformation == "Scaling"){
  #         current.transf <- "scaling"
  #       }
  #       else if(current.transformation == "Zohdy"){
  #         current.transf <- "zohdy"
  #       }
  #       else if(current.transformation == "Smoothed.Zohdy"){
  #         current.transf <- "smoothed_zohdy"
  #       }
  #       #
  #       base_transform <- "transform_"
  #       def_transform <- paste0(base_transform, current.transf)
  #       args <- list(ves = current.ves)
  #       res_transform <- do.call(def_transform, args)
  #       res <- data.frame(depth = res_transform$depth,
  #                         real.resistivity = res_transform$real.res)
  #       return(res)
  #     }
  #     
  #   }, extensions = c('Buttons'),
  #   options = list(
  #     pageLength = length(server.env$current.ves$ab2),
  #     dom = 'Bfrtip',
  #     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
  #     text = 'Download',
  #     scrollY = 200,
  #     scroller = TRUE))
  #
  ##########################################################################################
  #                             Baseflow Tab
  ##########################################################################################
  output$BFvarnames <- renderUI({
    if(!is.null(server.env$current.table)){
      current.names <-  c("None", server.env$current.varnames)
      selectInput(inputId= "BFvarnames1", label = "Discharge",
                  choices = current.names)
    }
  })
  #
  output$baseflow <- renderPlot({
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    current.var <- input$BFvarnames1
    if(is.null(current.var) || current.var == "None")
      return(NULL)
    Q <- as.numeric(as.matrix(current.table[current.var]))
    nd <- length(Q)
    #print(Q)
    alpha <- 0.0
    BSImax <- 0.0
    baseflow <- vector('numeric', length = nd)
    if(input$method == 'None'){
      return(NULL)
    }
    if(input$method == 'Graphical'){
      res_baseflow <- baseflow_graphical(Q)
      t <- res_baseflow$t
      baseflow <- c(res_baseflow$Qbf[2:nd],res_baseflow$Qbf[nd])
    }
    if(input$method == 'Nathan'){
      alpha1 <- as.numeric(input$nathan.alpha)
      baseflow <- nathan_filter(Q, alpha1)
    }
    if(input$method == 'Chapman'){
      alpha2 <- as.numeric(input$chapman.alpha)
      baseflow <- chapman_filter(Q, alpha2)
    }
    if(input$method == 'Eckhardt'){
      alpha <- as.numeric(input$eckhardt.alpha)
      BFImax <- as.numeric(input$eckhardt.bfi)
      baseflow <- eckhardt_filter(Q, alpha, BFImax)
    }
    #
    rtmn <- input$plot.range
    #
    # Create plot
    #
    nd <- length(Q)
    tt <- rep(seq(1, nd, 1),2)
    Q1 <- as.matrix(Q,nd,1)
    bf1 <- as.matrix(baseflow,nd,1)
    discharge <- rbind(Q1,bf1)
    Qlabel <- rep("Q", nd)
    bflabel <- rep("Baseflow", nd)
    process <- c(Qlabel, bflabel)
    baseflow.df <- data.frame(t = tt, discharge = discharge, process = process)
    p <- ggplot() + geom_line(aes(x = t, y = discharge, color = process), data = baseflow.df)
    return(p)
  })
  ##########################################################################################
  #                             Frequency Analysis Tab
  ##########################################################################################
  output$freq1 <- renderUI({
    if(!is.null(server.env$current.table)){
      current.names <-  c("None", server.env$current.varnames)
      selectInput(inputId= "freqvarnames", label = "Variable",
                  choices = current.names, selected = "None")
    }
  })
  #
  output$freq2 <- renderUI({
    current.table <- server.env$current.table
    tmp <- NULL
    if(!is.null(current.table)){
      if(input$freqselect == "SelectModel")
        tmp <- selectInput(inputId = "freqmethod", label = "Method",
                          choices = c("None", "Frequency Plot","Quantile Plot",
                                      "Moment Diagram", "L-Moment Diagram"))
      else if(input$freqselect == "ParameterEstimation"){
        tmp <- selectInput(inputId = "Parestmethod", label = "Method", choices =
                             c( None = "None", Moments="Moments", MLE = "MLE",
                                LMoments = "LMoments"))
      }
    }
    return(tmp)
  })
  #
  output$freq3 <- renderUI({
    current.table <- server.env$current.table
    tmp <- NULL
    if(!is.null(current.table)){
      if(input$freqselect == "SelectModel"){
        if(is.null(input$freqmethod) || input$freqmethod == "None"){
          return(NULL)
          #tmp <- NULL
        }
        if(input$freqmethod == "Frequency Plot"){
          tmp <- selectInput(inputId = "freqmodel", label = "Model", choices = model.types )
        }
        if(input$freqmethod == "Quantile Plot"){
          tmp <- selectInput(inputId = "freqmodel", label = "Model", choices = model.types )
        }
        else if(input$freqmethod == "Moment Diagram"){
          tmp <- NULL
        }
      }
      #
      if(input$freqselect == "ParameterEstimation"){
        tmp <- radioButtons(inputId = "freqmodel1", label = "PDF model",
                            choices = model.types)
      }
    }
    return(tmp)
  })
  #
  #
  output$frequency <- renderPlot({
    pfreq <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table))
      return(NULL)
    if(input$freqselect == "None")
      return(NULL)
    if(input$freqselect == "SelectModel"){
      if(is.null(input$freqmethod))
        return(NULL)
      if(input$freqmethod == "Frequency Plot"){
        current.model <- input$freqmodel
        if(is.null(current.model))
          return(NULL)
        if(current.model == "None" )
          return(NULL)
        current.var <- input$freqvarnames
        if(current.var == "None")
          return(NULL)
        var <- as.matrix(unname(current.table[current.var]))
        freq.results <- empirical_frequency(var, current.model)
        prob.results <- probability_plot(var, current.model)
        current.title <- paste0("Empirical Frequency Diagram: Model ", current.model)
        current.ylabel <- paste0("Standard ", current.model, " variable")
        Prob.df <- data.frame(var = prob.results$Var, model.var = prob.results$z)
        pfreq <- ggplot() + geom_point(aes(x = var, y = model.var), data = Prob.df) +
          geom_smooth(aes(x = var, y = model.var), data = Prob.df, method = "lm",
                      formula = y ~ x, se = FALSE)+
          xlab("Variable") + ylab(current.ylabel) +
          ggtitle(current.title)
        if(input$freqmodel == "lognormal" | input$freqmodel == "logpearson3")
          pfreq <- pfreq + scale_x_log10()
      }
      else if(input$freqmethod == "Quantile Plot"){
        current.model <- input$freqmodel
        if(is.null(current.model))
          return(NULL)
        if(current.model == "None" )
          return(NULL)
        current.var <- input$freqvarnames
        if(current.var == "None")
          return(NULL)
        var <- as.matrix(unname(current.table[current.var]))
        res.plplot <- probability_plot(var, current.model)
        q.plt.df <- data.frame(x = res.plplot$Var, y = res.plplot$z)
        pfreq <- ggplot() + geom_point(aes(x = x, y = y), data = q.plt.df)
      }
      else if(input$freqmethod == "Moment Diagram"){
        current.var <- input$freqvarnames
        if(current.var == "None")
          return(NULL)
        cs <- skewness(unname(current.table[current.var]))
        ck <- kurtosis(unname(current.table[current.var]))
        var <- data.frame(x=cs, y=ck)
        moment.diagram <- calculate_diagram_moments()
        LN3 <- as.data.frame(moment.diagram$LN3)
        P3 <- as.data.frame(moment.diagram$P3)
        GEV <- as.data.frame(moment.diagram$GEV)
        N <- data.frame(x=0,y=0)
        pfreq <- ggplot() + geom_line(aes(x = Cs, y = Ck), data = LN3, color = "red") +
          geom_line(aes(x = Cs, y = Ck), data = P3, color = "blue") +
          geom_line(aes(x = Cs, y = Ck), data = GEV, color = "green") +
          geom_point(aes(x = x, y = y), data = N, color = "black") +
          geom_point(aes(x = x, y = y), data = var, color = "yellow") +
          xlim(-2.5*cs, 2.5*cs) +
          ylim(0, 2.5*ck)
      }
      else if(input$freqmethod == "L-Moment Diagram"){
        current.var <- input$freqvarnames
        if(current.var == "None")
          return(NULL)
        tx <- as.matrix(current.table[current.var])
        sample.lmom <- samlmu(tx)
        pfreq <- lmrd(sample.lmom)
      }
    }
    return(pfreq)
  })
})
#
# output$parameter.estimates <- renderUI({
#   current.table <- server.env$current.table
#   if(is.null(current.table))
#     return(null)
#   if(input$freqselect == "ParameterEstimation"){
#
#   }
#
# })
