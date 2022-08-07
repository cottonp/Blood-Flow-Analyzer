#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
library(tidyverse)
library(plotly)
library(gridExtra)
library(segmented)
library(DT)
library(nls.multstart)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Blood Flow Analyzer"),

    # ---- Sidebar panel for inputs ----
    sidebarPanel(width = 3,
                 
                 # ---- Input: Select a file ----
                 fileInput("file", "Choose File:",
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv"
                           )),
                 bsTooltip("file", "Browse for time-intensity .csv file to process",
                           "right", options = list(container = "body")),
                 
                 # --- Input: Select region of interest ---
                 selectInput("roi_selection", "Region of interest:",
                             choices = c("Chest Wall", "Diaphragm", "Liver")),
                 bsTooltip("roi_selection", "Selects a region of interest to process",
                           "right", options = list(container = "body")),
                 
                 # --- Input: Slider to define baseline range ---
                 sliderInput("baseline_slider", "Baseline range:",
                             value = c(1,10),
                             min = 0,
                             max = 30,
                             round = -2,
                             step = 0.5
                             ),
                 bsTooltip("baseline_slider", "Use this slider to select a time range for the calculation of the baseline
                             acoustic intensity",
                           "right", options = list(container = "body")),
                 
                 # --- Input: Slider for max intensity ---
                 sliderInput("max_slider", "BFI anlysis range:",
                             value = 30,
                             min = 0,
                             max = 30,
                             round = -2,
                             step = 1
                 ),
                 bsTooltip("max_slider", "Use this slider to manually change the location of peak 
                             acoustic intensity"
                           ),
                 
                 # --- Input: Store BFI parameters ---
                 actionButton("store_bfi", label = "Store BFI Parameters"),
                 ),
    # ---- Main panel ----
    mainPanel(
        
        tabsetPanel(type = "tabs",
                    tabPanel("Blood Flow Analysis",
                             selectInput("regression_selection",
                                         "Regression method:",
                                         choices = c("Segmented", "Standard"),
                                         selected = "Standard"),
                             h2("BFI Analysis"),
                             plotlyOutput(outputId = "bfi_plot"),
                             textOutput("exclusion_comment"),
                             h3("BFI Parameters"),
                             DTOutput('bfi_table')
                             ),
                    tabPanel("Indicator-dilution Modelling",
                             br(),
                             actionButton("fit_model", label = "Fit models"),
                             selectInput("model_selection", "Select Models to Plot:",
                                         choices = c("Lognormal", 
                                                     "Gamma Variate", 
                                                     "LDRW",
                                                     "FPT",
                                                     "All"),
                                         selected = 'All'),
                             h3("Chest Wall"),
                             plotOutput(outputId = "cw_plot"),
                             h3("Diaphragm"),
                             plotOutput(outputId = "dia_plot"),
                             h3("Liver"),
                             plotOutput(outputId = "liver_plot"),
                             br(), DTOutput('modelling_table'),
                             DTOutput('cumulative_error'),
                             p('This table presents the cumulative and average residual standard error (RSE)
                               of each model across all three different regions of interest 
                               (chest wall, diaphragm and liver).'),
                             h2("Goodness of Fit"),
                             p("The purpose of this section is to visually assess whether or not the 
                               fitted model is violating any of the assumptions of non-linear 
                               least-squares regression."),
                             br(), plotlyOutput("standardized_residuals"),
                             p("These assumptions are:"),
                             strong("1) The variance of the error term is constant 
                                    (i.e., homoscedasticity)"),
                             br(),
                             strong("2) The expected value (i.e., average or mean) of the error
                                    term is equal to zero")),
                    tabPanel("Summary",
                             selectInput("model_final", "Select the Model you want to use:",
                                         choices = c("Lognormal", 
                                                     "Gamma Variate", 
                                                     "LDRW",
                                                     "FPT")),
                             h2("Summary Plots"),
                             plotlyOutput(outputId = "summary_plot"),
                             plotlyOutput(outputId = "summary_plot1"),
                             p('Note that all the traces plotted above (raw data, filtered data, and modelled data)
                               have been baseline and time zero (t0) corrected (utilizing the baseline
                               and t0 values obtained from the BFI analysis) to permit relative comparisons
                               between the different regions of interest. Similarly, the "Adjusted MTT" 
                               is the mean transit time (MTT) corrected to aforementioned t0 value.
                               Thus, this "Adjusted MTT" value should only be interpreted in the context
                               of the plots above.'),
                             h2("Summary Table"),
                             DTOutput(outputId = "summary_table"),
                             textInput('id', "ID:", value = " "),
                             textInput('condition', "Condition:", value = " "),
                             p("The ID and condition values are inferred from the file that is loaded, and
                             these values will be used to name the files that are saved from this 
                             analysis. However, this process is subject to error. Thus, you can 
                             change the value of these tags to anything you desire."),
                             actionButton('save_data', "Save Data"),
                             downloadButton('write_bf_table', "Write Blood Flow Data"),
                             downloadButton('write_TIC_data', "Write TIC Data"),
                             downloadButton('write_TIC_data_long', "Write long-form TIC Data"),
                             downloadButton('write_model_params', "Write Model Parameters"),
                             downloadButton('save_summary', 'Write Summary Report'),
                             br()
                             )
        ))
        # ---- Output: Tabset w/ parsing table, selected data, and plot ----
      
        )
    
    

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    values <- reactiveValues(
        df_in = NULL,
        df_out = NULL,
        region_dict = list("dia" = "Diaphragm",
                       "cw" = "Chest Wall",
                       'liver' = "Liver"),
        modelling_df = NULL,
        modelling_df_temp = NULL,
        modelling_table = NULL,
        summary_plot = NULL,
        summary_plot= NULL,
        roi = NULL,
        bfi_df = NULL,
        bfi_df_temp = NULL,
        fit_list = NULL,
        max_time = NULL,
        max_intensity = NULL,
        max_intensity_time = NULL,
        raw_data = NULL,
        filtered_data = NULL,
        fit_data = NULL,
        smoothed_data = NULL,
        regression_data = NULL,
        regression_t_zero = NULL,
        seg_slope_extrapolation = NULL,
        baseline = NULL,
        baseline_text = NULL,
        baseline_rect = NULL,
        start_rise = NULL,
        end_rise = NULL,
        rise_time = NULL,
        rise_time_text = NULL,
        rise_time_rect = NULL,
        peak_marker = NULL,
        peak_text = NULL,
        df_reg = NULL,
        #df_reg_plot = NULL,
        bfi_plot_shapes = NULL,
        bfi_plot_annotations = NULL,
        modelling_plot = NULL,
        modelling_plot1 = NULL,
        bfi = NULL,
        bf_table_out = NULL,
        TIC_data = NULL,
        TIC_data_long = NULL,
        model_fit_params = NULL
    )
    # ---- Observe file input and read in csv ----
    observeEvent(input$file$datapath, {
        req(input$file$datapath)
        
        tryCatch( {
            df <- read.csv(input$file$datapath)
            values$df <- df
            
            
        }, error = function(e) {
            stop(safeError(e))
        }
      
        
        )
      
      tryCatch({
        # split file name into ID and condition
        temp <- strsplit(input$file$name,"_")
        values$ID <- temp[[1]][1]
        values$condition <- temp[[1]][2]
        
        # update default ID and Condition inputs
        updateTextInput(session, 'id', 
                        value = values$ID)
        
        updateTextInput(session, 'condition', 
                        value = values$condition)
      },error = function(e) {
        stop(safeError(e))
      })
    })
    
    # ---- Observe changes to loaded data or region selection ----
    dataListen <- reactive({list(input$file, input$roi_selection)})
    observeEvent(dataListen(), {
        if(input$roi_selection == "Chest Wall"){
            roi <- 'cw'
        }else if(input$roi_selection == "Diaphragm"){
            roi <- 'dia'
        }else{
            roi <- 'liver'
        }
        values$raw_data <- roi; values$filtered_data <- paste(roi,"filt",sep = "_")
        values$fit_data <- paste(roi,"fit",sep = "_")
        values$smoothed_data <- paste(roi,"smooth",sep = "_")
        values$regression_data <- paste(roi, "regression",sep = "_")
        df <- values$df
        if(is.null(df)){
          
        }else{
          values$max_intensity <- max(df[values$raw_data])
          values$max_intensity_time <- df$time[which(df[values$raw_data] == values$max_intensity)]
          values$max_time <- df$time[length(df$time)]
        }
        
        values$roi <- roi
    })
    
    # ---- Observe changes in max intensity and max time: update slider values accordingly ----
    observe({
        req(values$max_time, values$max_intensity_time)
        updateSliderInput(session, 'baseline_slider', 
                          max = round(values$max_time))
        updateSliderInput(session, 'max_slider', value = round(values$max_time),
                          max = round(values$max_time))
    })
    
    
    # ---- Observe changes in regression selection ----
    observeEvent(input$regression_selection, {
        req(input$file)
        df <- values$df
        x1 <- input$baseline_slider[1]
        x2 <- input$baseline_slider[2]
        x3 <- input$max_slider
        start_rise <- values$start_rise
        end_rise <- values$end_rise
        baseline <- values$baseline
        max_intensity <- values$max_intensity
        # run regression for extrapolation of time zero
        if(input$regression_selection == "Segmented"){
            df_reg <- df %>% filter(time >= x1, time <=end_rise)
            lm_formula <- as.formula(paste0(values$raw_data,"~time"))
            my.lm <- lm(lm_formula, data = df_reg)
            seg <- segmented(my.lm,
                             seg.Z = ~time,
                             npsi = 1)
            b0 <- coef(seg)[[1]]
            b1 <- coef(seg)[[2]]
            c1 <- b1+coef(seg)[[3]]
            break_1 <- seg$psi[[2]]
            c0 <- b0+break_1*(b1-c1)
            df_reg[values$regression_data] <- predict(seg)
            t_zero <- (baseline-c0)/c1
            seg_slope_extrapolation <- df_reg %>% filter(time >= t_zero-2) %>% 
                mutate(y = c0+c1*time) %>% select(time, y)
            df_reg_plot <- df_reg
            values$seg_slope_extrapolation <- seg_slope_extrapolation
            
        }else{
            df_reg <- df %>% filter(time >= start_rise, time <= end_rise)
            my.lm <- lm(get(values$raw_data)~time, data = df_reg)
            b0 <- coef(my.lm)[1]; b1 <- coef(my.lm)[2]
            t_zero <- (baseline-b0)/b1
            df_reg_plot <- df %>% filter(time >= (t_zero-2), time <= end_rise)
            df_reg_plot[values$regression_data] <- (b0+b1*df_reg_plot$time)
        }
        
        # update t-zero annotation
        t_zero_text <- paste0('Time zero: ', round(t_zero,2),' s')
        t_zero_annotation <- list(text = t_zero_text, 
                                  xref = 'x', 
                                  yref = 'y',
                                  x = t_zero, 
                                  y = baseline,
                                  axref = 'x', 
                                  ayref = 'y',
                                  ax = t_zero, 
                                  ay = (max_intensity-baseline)/2+baseline,
                                  borderpad = 4,
                                  borderwidth = 2,
                                  bordercolor = 'rgb(0,0,0)',
                                  bgcolor = 'rgb(255,255,255)',
                                  arrowhead = 2,
                                  arrowwidth = 1.5,
                                  showarrow = T)
        
        # update reactive values
        values$bfi_plot_annotations[[1]] <- t_zero_annotation
        values$df_reg <- df_reg
        values$df_reg_plot <- df_reg_plot
        values$regression_t_zero <- t_zero
    })
 
    # ---- Observe changes in sliders and update graph accordingly ----
    sliderListen <- reactive({
        list(input$baseline_slider, input$max_slider)
    })
    observeEvent(sliderListen(),{
        req(input$file)
        df <- values$df
        x1 <- input$baseline_slider[1]
        x2 <- input$baseline_slider[2]
        x3 <- input$max_slider
        
        # calculate baseline
        
        baseline <- mean(df[values$raw_data][which((df$time >= x1) & (df$time <= x2)),])
        
        # calculate max intensity
        temp <- df[which((df$time <= x3) & (df$time >= x2)),]
        
        max_intensity <- max(temp[values$raw_data])
        max_intensity_time <- temp$time[which(temp[values$raw_data] == max_intensity)]
        
        # calculate intensity range
        intensity_range <- max_intensity - baseline
        
        # calculate 10 and 90% of the intensity range for the calculation of bfi
        ten_percent <- 0.1*intensity_range+baseline
        ninety_percent <- 0.9*intensity_range+baseline
        
        # determine the times associated with 10 and 90% of the intensity range
        start_rise <- temp$time[which(temp[values$raw_data]>=ten_percent)][1]
        end_rise <- temp$time[which(temp[values$raw_data]>=ninety_percent)][1]
        
        # calculate rise time and bfi
        rise_time <- end_rise-start_rise
        bfi <- intensity_range/rise_time
        
        # define colors to use for plotting
        blue_rgb <- 'rgb(41,110,175)'
        grey_rgb <- 'rgb(135,135,135)'
        red_rgb <- 'rgb(251,0,6)'
        salmon_rgb <- 'rgb(250, 147, 123)'
        black_rgb <- 'rgb(0,0,0)'
        white_rgb <- 'rgb(255,255,255)'
        
        # run regression for extrapolation of time zero
        if(input$regression_selection == "Segmented"){
         
            df_reg <- df %>% filter(time >= x1, time <=end_rise)
            lm_formula <- as.formula(paste0(values$raw_data,"~time"))
            my.lm <- lm(lm_formula, data = df_reg)
            seg <- segmented(my.lm,
                             seg.Z = ~time,
                             npsi = 1)
            b0 <- coef(seg)[[1]]
            b1 <- coef(seg)[[2]]
            c1 <- b1+coef(seg)[[3]]
            break_1 <- seg$psi[[2]]
            c0 <- b0+break_1*(b1-c1)
            df_reg[values$regression_data] <- predict(seg)
            t_zero <- (baseline-c0)/c1
            seg_slope_extrapolation <- df_reg %>%filter(time >=t_zero-2) %>% 
                mutate(y = c0+c1*time) %>% select(time, y)
            df_reg_plot <- df_reg
            values$seg_slope_extrapolation <- seg_slope_extrapolation
        }else{
            df_reg <- df %>% filter(time >= start_rise, time <= end_rise)
            my.lm <- lm(get(values$raw_data)~time, data = df_reg)
            b0 <- coef(my.lm)[1]; b1 <- coef(my.lm)[2]
            t_zero <- (baseline-b0)/b1
            df_reg_plot <- df %>% filter(time >= (t_zero-2), time <= end_rise)
            df_reg_plot[values$regression_data] <- (b0+b1*df_reg_plot$time)
        }
        
        # update reactive values
        values$df_reg <- df_reg
        values$df_reg_plot <- df_reg_plot
        values$regression_t_zero <- t_zero
       
        
        # create annotations for bfi plot
        # create annotation for marker 3 or "peak" marker
        arrow_x <-  ((values$max_time - values$max_intensity_time)/2)+values$max_intensity_time
        peak_text <- list(
                          y = max_intensity, 
                          yref = 'y',
                          xref = 'x',
                          x = arrow_x,
                          text = paste0('Peak = ', round(max_intensity,2),' dB'),
                          font = list(color = black_rgb),
                          showarrow = F,
                          yshift = 10,
                          opacity = 0.6
                          )
        
        baseline_annotation <- list(
            y = baseline, 
            yref = 'y',
            xref = 'x',
            x = arrow_x,
            text = paste0('Baseline = ', round(baseline,2),' dB'),
            font = list(color = black_rgb),
            showarrow = F,
            yshift = -10,
            opacity = 0.6
        )
        
        # create baseline text annotation
        baseline_text <- list(x = ((x1+x2)/2), xref = 'x',
                              y = 0, yref = 'paper',
                              text = 'Baseline range',
                              font = list(color = blue_rgb),
                              showarrow = F)
        
        # create rise time text annotation
        rise_time_text <- paste0("Rise time = ", as.character(round(rise_time,2))," s")
        rise_time_text <- list(x = ((rise_time/2)+start_rise),
                               y = 0,
                               yref = 'paper',
                               xref = 'x',
                               text = rise_time_text,
                               font = list(color = salmon_rgb),
                               showarrow = F)
        
        # create intensity range annotation
        
        intensity_range_arrow <- list(xref = 'x', 
                                      yref = 'y', 
                                      axref = 'x', 
                                      ayref = 'y',
                                      x = arrow_x, 
                                      y = baseline, 
                                      ax = arrow_x, 
                                      ay = max_intensity, 
                                      showarrow = T,
                                      startarrowhead = 2, 
                                      arrowhead = 2, 
                                      arrowside = 'end+start',
                                      arrowwidth = 1.5,
                                      arrowcolor = black_rgb
                                      )
        
        intensity_range_text <- paste0('Intensity range = ', round(intensity_range,2), ' dB')
        intensity_range_text <- list(x = arrow_x, xref = 'x',
                                     y = (max_intensity-baseline)/2+baseline,
                                     yref = 'y',
                                     borderpad = 4,
                                     borderwidth = 2,
                                     showarrow = F,
                                     bgcolor = white_rgb,
                                     bordercolor = black_rgb,
                                     text = intensity_range_text
                                     )
        # create time_zero annotation
        t_zero_text <- paste0('Time zero: ', round(t_zero,2),' s')
        t_zero_annotation <- list(text = t_zero_text, 
                                  xref = 'x', 
                                  yref = 'y',
                                  x = t_zero, 
                                  y = baseline,
                                  axref = 'x', 
                                  ayref = 'y',
                                  ax = t_zero, 
                                  ay = (max_intensity-baseline)/2+baseline,
                                  borderpad = 4,
                                  borderwidth = 2,
                                  bordercolor = black_rgb,
                                  bgcolor = white_rgb,
                                  arrowhead = 2,
                                  arrowwidth = 1.5,
                                  showarrow = T)
        
        # create bfi text annotation
        bfi_text <- paste0('BFI = ', round(bfi,2), ' dB/s')
        bfi_text <- list(x = arrow_x, xref = 'x',
                         y = 0,
                         yref = 'paper',
                         text = bfi_text,
                         showarrow = F
                         )
                         
        # create shapes for plot
        # peak marker
        peak_filter <- list(type = "line",
                            x0 = x3, x1 = x3, xref = 'x',
                            y0 = 0, y1 = 1, yref = 'paper',
                            opacity = 0.6,
                            layer = 'below',
                            line = list(color = black_rgb)
                            )
        peak_marker <- list(type = "line",
                            x0 = max_intensity_time, 
                            x1 = max_intensity_time, 
                            xref = 'x',
                            y0 = 0, y1 = 1, yref = 'paper',
                            line = list(color = salmon_rgb,
                                        dash = 'dash'),
                            name = "Peak acoustic intensity")
        
        # rise time span
        rise_time_rect <- list(type = "rect",
                               fillcolor = salmon_rgb,
                               line = list(color = salmon_rgb),
                               opacity = 0.2,
                               x0 = start_rise, x1 = end_rise, xref = "x",        
                               y0 = 0, y1 = 1, yref = "paper"
                               )
        
        # baseline span
        baseline_rect <- list(type = "rect",
                              fillcolor = blue_rgb,
                              line = list(color = blue_rgb),
                              opacity = 0.2,
                              x0 = x1, x1 = x2, xref = "x",
                              y0 = 0, y1 = 1, yref = "paper"
                              )
        
        # exclusion span
        exclusion_rect <- list(type = "rect",
                               fillcolor = grey_rgb,
                               line = list(color = grey_rgb),
                               opacity = 0.1,
                               x0 = x3, 
                               x1 = values$max_time, 
                               xref = "x",
                               y0 = 0, y1 = 1, yref = "paper"
        )
        # baseline reference line
        baseline_trace <- list(type = "line",
                             x0 = 0, x1 = 1, xref = 'paper',
                             y0 = baseline, y1 = baseline, yref = 'y',
                             line = list(dash = 'longdash'),
                             opacity = 0.8
                             )
        
        # max intensity reference line
        max_intensity_trace <- list(type = "line",
                                    x0 = max_intensity_time, 
                                    x1 = values$max_time, 
                                    xref = 'x',
                                    y0 = max_intensity, 
                                    y1 = max_intensity, 
                                    yref = 'y',
                                    line = list(dash = 'longdash'),
                                    opacity = 0.8
                                    )

        # set reactive values
        values$start_rise <- start_rise
        values$end_rise <- end_rise
        values$baseline <- baseline
        values$max_intensity <- max_intensity
        values$rise_time <- rise_time
        values$bfi <- bfi
        values$bfi_plot_annotations <- list(t_zero_annotation,baseline_text, rise_time_text, bfi_text, 
                                             intensity_range_arrow,intensity_range_text,
                                            peak_text,baseline_annotation)
        values$bfi_plot_shapes <- list(baseline_rect, rise_time_rect, baseline_trace,
                                       max_intensity_trace, peak_marker, exclusion_rect,
                                       peak_filter)
        
    })
    
    
    
    # ---- Fit model button ----
    observeEvent(input$fit_model, {
      req(input$file)
      withProgress(message = "Fitting in Progress", {
      incProgress(amount = 0)
      df <- values$df
      regions <- c('cw','dia','liver')
      
      for (i in 1:length(regions)) {
   
      # define region of interest
      roi <- regions[i]
      
      # filter bfi df to current roi and set starting parameter limits
      bfi_df <- values$bfi_df %>% filter(region == values$region_dict[roi])
      upper_t <- bfi_df$t0*1.5
      upper_C <- 2*bfi_df$baseline
      lower_C <- 0
      
      
      # define different modelling functions
      log_norm <- function(AUC, t, u, s, t0, C) ifelse (t0 >= t, C, AUC/(sqrt(2*pi)*s*(t-t0))*exp(-((log(t-t0)-u)^2)/(2*s^2)) + C)
      gamma_variate <- function(AUC, t, a, b, t0, C) ifelse (t0 >= t, C, AUC/(b^(a+1)*gamma(a+1))*((t-t0)^a)*exp(-(t-t0)/b) + C)
      ldrw <- function(AUC, t, u, lambda, t0, C) ifelse (t0 >= t, C, AUC*exp(lambda)/u*sqrt(u*lambda/((t-t0)*2*pi))*exp(-1/2*lambda*((u/(t-t0))+(t-t0)/u)) + C)
      fpt<- function(AUC, t, u, lambda, t0, C) ifelse (t0 >= t, C, AUC*exp(lambda)/u*sqrt(lambda/2*pi)*(u/(t-t0))^(3/2)*exp(-1/2*lambda*((u/(t-t0))+(t-t0)/u)) + C)
      
     
      
      lognormal_plot <- paste(roi,"lognormal", sep = '_')
      gamma_plot <- paste(roi, 'gamma',sep = '_')
      ldrw_plot <- paste(roi,'ldrw',sep = '_')
      fpt_plot <- paste(roi, 'fpt', sep = '_')
      
      
      # fit lognormal model
      nls_formula <- as.formula(substitute(y ~ log_norm(AUC, t = time, u, s, t0, C), list(y = as.name(roi))))
      lognormal_fit <- nls_multstart(nls_formula, 
                            data = df,
                            iter = 200,
                            start_lower = c(AUC = 0, u = 0, s = 0, t0 = 0, C = lower_C),
                            start_upper = c(AUC = 10, u = 10, s = 1, t0 = upper_t, C = upper_C),
                            supp_errors = "Y",
                            lower = c(AUC = 0, u = 0, s = 0, t0 = 0, C = 0)
      )
      
    
      df[lognormal_plot] <- predict(lognormal_fit)
      u <- coef(lognormal_fit)[2]
      s <- coef(lognormal_fit)[3]
      t_zero <- coef(lognormal_fit)[4]
      baseline <- coef(lognormal_fit)[5]
      MTT <- exp(u + s^2 / 2)
      Tp <- exp(u - s^2)
      AUC <- coef(lognormal_fit)[1]
      AUC_inverse <- 1/AUC
      volume <- MTT/AUC
      modelling_table <- data.frame(Model = "Lognormal",
                               RSE = sigma(lognormal_fit),
                               t0 = t_zero,
                               AUC = AUC,
                               MTT = MTT,
                               Tp = Tp,
                               Baseline = baseline,
                               AUC_inverse = AUC_inverse,
                               Volume = volume
                               )
      
      
      # fit gamma variate model
      nls_formula <- as.formula(substitute(y ~ gamma_variate(AUC, t = time, a, b, t0, C), list(y = as.name(roi))))
      gammavariate_fit <- nls_multstart(nls_formula, 
                                        data = df,
                                        iter = 200,
                                        start_lower = c(AUC = 0, a = 0, b = 0, t0 = 0, C = lower_C),
                                        start_upper = c(AUC = 10, a = 10, b = 10, t0 = upper_t, C = upper_C),
                                        supp_errors = "Y",
                                        lower = c(AUC = 0, a = 1, b = 1, t0 = 0, C = 0)
      )
      
      
      df[gamma_plot] <- predict(gammavariate_fit)
      alpha <- coef(gammavariate_fit)[2]
      beta <- coef(gammavariate_fit)[3]
      t0 <- coef(gammavariate_fit)[4]
      baseline <- coef(gammavariate_fit)[5]
      RSE <- sigma(gammavariate_fit)
      MTT <- beta*(alpha+1)
      Tp <- alpha*beta
      AUC <- coef(gammavariate_fit)[1]
      AUC_inverse <- 1/AUC
      volume <- MTT/AUC
      modelling_table <- rbind(modelling_table, c("Gamma Variate", RSE,t0, AUC, MTT,Tp,baseline, AUC_inverse,
                                                  volume))
      
      
      # fit ldrw model
      nls_formula <- as.formula(substitute(y ~ ldrw(AUC, t = time, u, lambda, t0, C), list(y = as.name(roi))))
      ldrw_fit <- nls_multstart(nls_formula, 
                                data = df,
                                iter = 200,
                                start_lower = c(AUC = 0, u = 0, lambda = 0, t0 = 0, C = lower_C),
                                start_upper = c(AUC = 10, u = 10, lambda = 10, t0 = upper_t, C = upper_C),
                                supp_errors = "Y",
                                lower = c(AUC = 0, u = 1, lambda = 1, t0 = 0, C = 0)
      )
      
      df[ldrw_plot] <- predict(ldrw_fit)
      u <- coef(ldrw_fit)[2]
      lambda <- coef(ldrw_fit)[3]
      t0 <- coef(ldrw_fit)[4]
      RSE <- sigma(ldrw_fit)
      MTT <- u
      Tp <- (u/2*lambda)*(sqrt(1+4*lambda^2)-1)
      AUC <- coef(ldrw_fit)[1]
      AUC_inverse <- 1/AUC
      volume <- MTT/AUC
      baseline <- coef(ldrw_fit)[5]
      modelling_table <- rbind(modelling_table, c("LDRW",RSE,t0, AUC, MTT,Tp, baseline, AUC_inverse,
                                                  volume))
      
      # fit fpt model
      nls_formula <- as.formula(substitute(y ~ fpt(AUC, t = time, u, lambda, t0, C), list(y = as.name(roi))))
      fpt_fit <- nls_multstart(nls_formula, 
                               data = df,
                               iter = 200,
                               start_lower = c(AUC = 0, u = 0, lambda = 0, t0 = 0, C = 0),
                               start_upper = c(AUC = 10, u = 10, lambda = 10, t0 = 15, C = 80),
                               supp_errors = "Y",
                               lower = c(AUC = 0, u = 1, lambda = 1, t0 = 0, C = 0)
      )
      
      
      df[fpt_plot] <- predict(fpt_fit)
      u <- coef(fpt_fit)[2]
      lambda <- coef(fpt_fit)[3]
      t0 <- coef(fpt_fit)[4]
      MTT <- u
      Tp <- (u/2*lambda)*(sqrt(9+4*lambda^2)-3)
      AUC <- coef(fpt_fit)[1]
      AUC_inverse <- 1/AUC
      volume <- MTT/AUC
      baseline <- coef(fpt_fit)[5]
      modelling_table <- rbind(modelling_table, c("FPT",RSE,t0, AUC, MTT,Tp, baseline, AUC_inverse,
                                                  volume))
      
      
      # save fit stats for each model
      fit_list <- list(lognormal = lognormal_fit,
                              gamma = gammavariate_fit,
                              ldrw = ldrw_fit,
                              fpt = fpt_fit)
      # save the fits from all models to fit list grouped by region
      values$fit_list[[roi]] <- fit_list
      
      modelling_table$MTT <- as.numeric(modelling_table$MTT)
      modelling_table$Tp <- as.numeric(modelling_table$Tp)
      modelling_table$t0 <- as.numeric(modelling_table$t0)
      modelling_table$Baseline <- as.numeric(modelling_table$Baseline)
      modelling_table$Volume <- as.numeric(modelling_table$Volume)
      modelling_table$AUC <- as.numeric(modelling_table$AUC)
      modelling_table$AUC_inverse <- as.numeric(modelling_table$AUC_inverse)
      modelling_table$Peak <- modelling_table$t0+modelling_table$Tp
      modelling_table$Region <- values$region_dict[[roi]]
      modelling_table <- modelling_table %>% select(Region, everything())
      rownames(modelling_table) <- 1:nrow(modelling_table)
      
      if(is.null(values$modelling_table)){
        values$modelling_table <- modelling_table  
      }else{
        values$modelling_table <- rbind(values$modelling_table, modelling_table)
      }
      # set progress bar  
      incProgress(amount = 0.2)
      }
      
      values$modelling_table$Region <- as.factor(values$modelling_table$Region)
      values$modelling_table$Model <- as.factor(values$modelling_table$Model)
      values$modelling_table$RSE <- as.numeric(values$modelling_table$RSE)
      
      dia_df <- df %>% select(time, dia, dia_filt, dia_lognormal, dia_gamma, dia_ldrw, dia_fpt) %>% 
        mutate(region = "Diaphragm",
               raw = dia,
               dia = NULL,
               filtered = dia_filt,
               dia_filt = NULL,
               lognormal = dia_lognormal,
               dia_lognormal = NULL,
               gamma_variate = dia_gamma,
               dia_gamma = NULL,
               ldrw = dia_ldrw,
               dia_ldrw = NULL,
               fpt = dia_fpt,
               dia_fpt = NULL) %>% 
        select(region, everything())
      
      cw_df <- df %>% select(time, cw, cw_filt, cw_lognormal, cw_gamma, cw_ldrw, cw_fpt) %>% 
        mutate(region = "Chest Wall",
               raw = cw,
               cw = NULL,
               filtered = cw_filt,
               cw_filt = NULL,
               lognormal = cw_lognormal,
               cw_lognormal = NULL,
               gamma_variate = cw_gamma,
               cw_gamma = NULL,
               ldrw = cw_ldrw,
               cw_ldrw = NULL,
               fpt = cw_fpt,
               cw_fpt = NULL) %>% 
        select(region, everything())
      
      liver_df <- df %>% select(time, liver, liver_filt, liver_lognormal, liver_gamma, liver_ldrw,
                                liver_fpt) %>% 
        mutate(region = "Liver",
               raw = liver,
               liver = NULL,
               filtered = liver_filt,
               liver_filt = NULL,
               lognormal = liver_lognormal,
               liver_lognormal = NULL,
               gamma_variate = liver_gamma,
               liver_gamma = NULL,
               ldrw = liver_ldrw,
               liver_ldrw = NULL,
               fpt = liver_fpt,
               liver_fpt = NULL) %>% 
        select(region, everything())
      
      modelling_df_long <- rbind(dia_df, cw_df, liver_df)
      modelling_df_long$region <- as.factor(modelling_df_long$region)
      values$modelling_df_long <- modelling_df_long
      values$df <- df
      
      # create residual df
      model_names <- c('lognormal', 'gamma', 'ldrw', 'fpt')
      df_resid <- modelling_df_long %>% select(region,time)
      df_resid_std <- modelling_df_long %>% select(region,time)
      regions <- c('dia','cw','liver')
      for (i in 1:length(values$fit_list)) {
        for (e in 1:length(model_names)) {
          df_resid[which(df_resid$region == values$region_dict[regions[i]]),model_names[e]] <- resid(values$fit_list[[regions[i]]][[model_names[e]]])
          df_resid_std[which(df_resid_std$region == values$region_dict[regions[i]]),model_names[e]] <- resid(values$fit_list[[regions[i]]][[model_names[e]]])/sigma(values$fit_list[[regions[i]]][[model_names[e]]])
        }
      }
      values$df_resid <- df_resid
      values$df_resid_std <- df_resid_std
      
      })
    })
    
    # create standardized residual plot output
    output$standardized_residuals <- renderPlotly({
      req(values$df_resid_std)
      ggplotly(values$df_resid_std %>% ggplot(aes(x = time))+
                 geom_point(aes(y = lognormal, color = "Lognormal"), alpha = 0.1)+
                 geom_point(aes(y = gamma, color = "Gamma Variate"), alpha = 0.1)+
                 geom_point(aes(y = ldrw, color = 'LDRW'), alpha = 0.1)+
                 geom_point(aes(y=fpt, color = "FPT"), alpha = 0.1)+facet_wrap(~region)+
                 geom_hline(yintercept = 0, color = 'black', linetype = 'dashed')+
                 geom_hline(yintercept = 3, color = 'black', linetype = 'dotted', alpha = 0.7)+
                 geom_hline(yintercept = -3, color = 'black', linetype = 'dotted', alpha = 0.7)+
                 ylab("Standardized Residuals")+xlab("Time (s)")+
                 scale_color_manual(values = c('Lognormal' = '#e41a1c', 
                                               'Gamma Variate' = '#377eb8',
                                               'LDRW' = '#4daf4a',
                                               'FPT'=  '#984ea3'))+
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.minor = element_blank(), legend.title = element_blank()))
      
      
    })
    
    # create cumulative error output
    output$cumulative_error <- renderDT({
      req(values$modelling_table)
      table_out <- values$modelling_table %>% select(Model,Region, RSE) %>% 
        group_by(Model) %>% summarize(sum = round(sum(RSE),2), sum1 = round(sum(RSE)/3,2))
      colnames(table_out) <- c("Model", "Cumulative RSE", "Average RSE")
      table_out
    })
    
    # create modelling table output
    output$modelling_table <- renderDT({
      table_out <- values$modelling_table
      table_out[,c('RSE','t0','MTT','AUC','Tp','Baseline','Peak', 'AUC_inverse', 'Volume')] <- round(table_out[,c('RSE','t0','MTT','AUC','Tp','Baseline','Peak', 'AUC_inverse', 'Volume')],2)
      colnames(table_out) <- c('Region', 'Model', 'RSE', 't0', 'AUC','MTT', 'Tp', 'Baseline','1/AUC', 'MTT/AUC', 'Peak time (t0 + Tp)')
      table_out
      })
    
    # create bfi_table output
    observe({
      req(values$start_rise, values$end_rise, values$rise_time, values$max_intensity,
          values$baseline, values$regression_t_zero, values$bfi)
      bfi_table <- data.frame(region = input$roi_selection,
                              bfi = round(values$bfi,2),
                              int_range = round((values$max_intensity - values$baseline),2),
                              baseline = round(values$baseline,2),
                              peak = round(values$max_intensity,2),
                              t0 = round(values$regression_t_zero,2),
                              rise_time = round(values$rise_time,2),
                              start_rise = round(values$start_rise,2),
                              end_rise = round(values$end_rise,2),
                              baseline_range_start = round(input$baseline_slider[1],2),
                              baseline_range_end = round(input$baseline_slider[2],2)
      )
      

      values$bfi_df_temp <- bfi_table
    })
    
    # ---- Observe model selection ----
    observeEvent(input$model_selection, {
      req(values$fit_list, input$model_selection != "All")
      if(input$model_selection == "Lognormal"){
        model <- 'lognormal'
      }else if(input$model_selection == "Gamma Variate"){
       model <- 'gamma'
      }else if(input$model_selection == 'LDRW'){
        model <- 'ldrw'
      }else if(input$model_selection == 'FPT'){
        model <- 'fpt'
      }
      
      
      # create chest wall plot
      plot_var <- paste('cw', model, sep = '_')
      rug_data <- values$modelling_table %>% filter(Model == input$model_selection, Region == "Chest Wall")
      values$cw_plot <- values$df %>% ggplot(aes(x = time))+
        geom_point(aes(y = cw, color = "Raw data"), size = 1)+
        geom_line(aes(y = get(plot_var), color = eval(input$model_selection)), size = 1, alpha = 0.7)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x =  MTT, linetype = "MTT"), color = 'black',
                 alpha = 0.8, length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = Peak, linetype = "t0 + Tp"), color = 'black', alpha = 0.8,
                 length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = t0, linetype = "t0"), color = 'black', alpha = 0.8, size = 1,
                 length = unit(0.3, 'npc'))+
        xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
        ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
        scale_color_manual(values = c('blue','grey', 'black'), aesthetics = c('color'))+
        coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
        scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
      
      # create diaphragm plot
      plot_var_dia <- paste('dia', model, sep = '_')
      rug_data <- values$modelling_table %>% filter(Model == input$model_selection, Region == "Diaphragm")
      values$dia_plot <- values$df %>% ggplot(aes(x = time))+
        geom_point(aes(y = dia, color = "Raw data"), size = 1)+
        geom_line(aes(y = get(plot_var_dia), color = eval(input$model_selection)), size = 1, alpha = 0.7)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x =  MTT, linetype = "MTT"), color = 'black',
                 alpha = 0.8, length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = Peak, linetype = "t0 + Tp"), color = 'black', alpha = 0.8,
                 length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = t0, linetype = "t0"), color = 'black', alpha = 0.8, size = 1,
                 length = unit(0.3, 'npc'))+
        xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
        ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
        scale_color_manual(values = c('blue','grey', 'black'), aesthetics = c('color'))+
        coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
        scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
      
      # create liver plot
      plot_var_liver <- paste('liver', model, sep = '_')
      rug_data <- values$modelling_table %>% filter(Model == input$model_selection, Region == "Liver")
      values$liver_plot <- values$df %>% ggplot(aes(x = time))+
        geom_point(aes(y = liver, color = "Raw data"), size = 1)+
        geom_line(aes(y = get(plot_var_liver), color = eval(input$model_selection)), size = 1, alpha = 0.7)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x =  MTT, linetype = "MTT"), color = 'black',
                 alpha = 0.8, length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = Peak, linetype = "t0 + Tp"), color = 'black', alpha = 0.8,
                 length = unit(1, 'npc'), size = 1)+
        geom_rug(data = rug_data, inherit.aes = F,
                 aes(x = t0, linetype = "t0"), color = 'black', alpha = 0.8, size = 1,
                 length = unit(0.3, 'npc'))+
        xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
        ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
        scale_color_manual(values = c('blue','grey', 'black'), aesthetics = c('color'))+
        coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
        scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
      
      
    })

    # ---- Observe "changes to loaded data or region selection "store bfi" button ----
    observeEvent(input$store_bfi, {
      if(is.null(values$bfi_df)){
        values$bfi_df <- values$bfi_df_temp
      }else{
        if(isTRUE(nrow(values$bfi_df[which(values$bfi_df$region == input$roi_selection),]) == 0)){
          values$bfi_df <- rbind(values$bfi_df, values$bfi_df_temp)
        }else{
          values$bfi_df[which(values$bfi_df$region == input$roi_selection),] <- values$bfi_df_temp[1,]
        }
      }
      
      rownames(values$bfi_df) <- 1:nrow(values$bfi_df)
      }
      )
    
    # write blood flow table
    output$write_bf_table <- downloadHandler(
     
      filename = '_blood_flow_table.csv',
      content = function(file) {
        write.csv(values$bf_table_out, file)
      }
    )
    # write TIC data
    output$write_TIC_data <- downloadHandler(
      
      filename = '_TIC_data.csv',
      content = function(file) {
        write.csv(values$TIC_data, file)
      }
    )
    # write TIC data long 
    output$write_TIC_data_long <- downloadHandler(
      
      filename = '_TIC_data_long.csv',
      content = function(file) {
        write.csv(values$TIC_data_long, file)
      }
    )
    # fit stats
    output$write_model_params <- downloadHandler(
      
      filename = '_model_params.csv',
      content = function(file) {
        write.csv(values$model_fit_params, file)
      }
    )
    # ---- Observe save data button ----
    observeEvent(input$save_data, {
      file_tag <- paste('./output/',input$id,sep = '_',input$condition)
      values$file_tag <- file_tag
      model_selection <- input$model_final
      bfi <- values$bfi_df
      colnames(bfi) <- c('region','bfi','int_range', 'baseline_bfi', 'peak_intensity', 't0_bfi',
                         'rise_time', 'start_rise', 'end_rise', 'bfi_baseline_range_start',
                         'bfi_baseline_range_end')
      model <- values$modelling_table %>% filter(Model == eval(model_selection)) %>% 
        select(Region, MTT, Tp, t0, Baseline, Model)
      colnames(model) <- c('region', 'MTT', 'Tp', 't0_model', 'baseline_model','model_name')
      blood_flow_table <- inner_join(bfi,model, by = 'region')
      
      # write blood flow table
      write.csv(blood_flow_table, file = paste0(file_tag, '_blood_flow_table.csv'))
      # save table for writing
      values$bf_table_out <- blood_flow_table
      
      # write tic data
      regions <- c('cw','dia','liver')
      model_dict <- list('Lognormal' = 'lognormal', 'Gamma Variate' = 'gamma',
                         'LDRW' = 'ldrw', 'FPT' = 'fpt')
      model_column <- model_dict[[model_selection]]
      
      for (i in 1:length(regions)) {
        regions[i] <- paste(regions[i], model_column, sep = '_')
      }
      columns <- c('time', 'dia','cw','liver', 'dia_filt','cw_filt','liver_filt', regions)
      
      df_out <- values$df[,columns]
      write.csv(df_out, file = paste(file_tag, 'TIC_data.csv', sep = '_'))
      # save TIC data for writing
      values$TIC_data <- df_out
      
      # write long tic data
      model_dict <- list('Lognormal' = 'lognormal', 'Gamma Variate' = 'gamma_variate',
                         'LDRW' = 'ldrw', 'FPT' = 'fpt')
      model_column <- model_dict[[model_selection]]
      
      df_out <- values$modelling_df_long %>% select(region, time, raw, filtered, as.name(model_column))
      write.csv(df_out, file = paste(file_tag, 'TIC_data_long.csv', sep = '_'))
      # save long TIC data for writing
      values$TIC_data_long <- df_out
      
      # write fit stats
      regions <- c('cw','dia','liver')
      model_dict <- list('Lognormal' = 'lognormal', 'Gamma Variate' = 'gamma',
                         'LDRW' = 'ldrw', 'FPT' = 'fpt')
      model_selection
      temp <- values$fit_list
      df_out <- NULL
      for (i in regions) {
       #temp1 <- as.data.frame(summary(temp[[i]][[model_dict[[model_selection]]]]))
       temp1 <- temp[[i]][[model_dict[[model_selection]]]]
       temp1 <- summary(temp1)
       temp1 <- temp1$coefficients
       temp1 <- as.data.frame(temp1)
       temp1$region <- i
       if(is.null(df_out)){
         df_out <- temp1
       }else{
         df_out <- rbind(df_out, temp1)
       }
      }
      write.csv(df_out, file = paste(file_tag, 'model_parameter_fit_stats.csv', sep = '_'))
      # save fit stats for writing
      values$model_fit_params <- df_out
    })
    
    # ---- save summary button ----
    output$save_summary <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.rmd")
        file.copy("report.rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(plot1 = values$summary_plot,
                       plot2 = values$summary_plot1,
                       table = values$summary_table,
                       ID = input$id,
                       condition = input$condition)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )

    # build output for bfi table
    output$bfi_table <- renderDT({
      req(values$bfi_df_temp)
      rownames(values$bfi_df_temp) <- "Temporary values"
      if(is.null(values$bfi_df)){
        table_out <- values$bfi_df_temp %>% select(region, bfi, int_range, baseline, peak,
                                                   rise_time, t0, start_rise, end_rise)
      }else{
        table_out <- values$bfi_df
        table_out <- rbind(table_out, values$bfi_df_temp)%>% select(region, bfi, int_range, baseline, peak,
                                                                        rise_time, t0, start_rise, end_rise)
      }
      colnames(table_out) <- c('Region', 'BFI (dB/s)', 'Intensity Range (dB)', 'Baseline (dB)',
                               'Peak Intensity (dB)', 'Rise time (s)','t0 (s)', 'Rise Time Start (s)',
                               'Rise Time End (s)')
      table_out
    })

    # set output for cw plot
    output$cw_plot <- renderPlot({
      req(values$modelling_df_long, values$modelling_table)
      
      if(input$model_selection != 'All'){
        values$cw_plot
      }else{
        # build chest wall plot
        rug_data <- values$modelling_table %>% filter(Region == 'Chest Wall')
        cw_plot <- values$modelling_df_long %>% filter(region == "Chest Wall") %>% 
          ggplot(aes(x = time))+
          geom_point(aes(y = raw, color = "Raw data"), size = 1, alpha = 0.7)+
          geom_line(aes(y = lognormal, color = "Lognormal"), size = 1, alpha = 0.7)+
          geom_line(aes(y = gamma_variate, color = "Gamma Variate"), size = 1, alpha = 0.7)+
          geom_line(aes(y=ldrw, color = "LDRW"), size = 1, alpha = 1)+
          geom_line(aes(y = fpt, color = "FPT"), size = 1, alpha = 0.5)+
          geom_rug(data =rug_data, inherit.aes = F, 
                   aes(x =  MTT, linetype = "MTT", color = Model),
                   alpha = 0.5, length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F, 
                   aes(x = Peak, linetype = "t0 + Tp", color = Model), alpha = 0.5,
                   length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F,
                   aes(x = t0, linetype = "t0", color = Model), alpha = 0.5, size = 1,
                   length = unit(0.3, 'npc'))+
          xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
          ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
          scale_color_manual(values = c('Raw data' = 'grey', 'Lognormal' = '#e41a1c', 
                                        'Gamma Variate' = '#377eb8',
                                        'LDRW' = '#4daf4a',
                                        'FPT'=  '#984ea3'))+
          coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
          scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
        cw_plot
      }
     
    })
    
    # set output for dia plot
    output$dia_plot <- renderPlot({
      req(values$modelling_df_long, values$modelling_table)
      
      if(input$model_selection != 'All'){
        values$dia_plot
      }else{
        # build diaphragm plot
        rug_data <- values$modelling_table %>% filter(Region == 'Diaphragm')
        dia_plot <- values$modelling_df_long %>% filter(region == "Diaphragm") %>% 
          ggplot(aes(x = time))+
          geom_point(aes(y = raw, color = "Raw data"), size = 1, alpha = 0.7)+
          geom_line(aes(y = lognormal, color = "Lognormal"), size = 1, alpha = 0.7)+
          geom_line(aes(y = gamma_variate, color = "Gamma Variate"), size = 1, alpha = 0.7)+
          geom_line(aes(y=ldrw, color = "LDRW"), size = 1, alpha = 1)+
          geom_line(aes(y = fpt, color = "FPT"), size = 1, alpha = 0.5)+
          geom_rug(data = rug_data, inherit.aes = F, 
                   aes(x =  MTT, linetype = "MTT", color = Model),
                   alpha = 0.5, length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F, 
                   aes(x = Peak, linetype = "t0 + Tp", color = Model), alpha = 0.5,
                   length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F,
                   aes(x = t0, linetype = "t0", color = Model), alpha = 0.5, size = 1,
                   length = unit(0.3, 'npc'))+
          xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
          ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
          scale_color_manual(values = c('Raw data' = 'grey', 'Lognormal' = '#e41a1c', 
                                        'Gamma Variate' = '#377eb8',
                                        'LDRW' = '#4daf4a',
                                        'FPT'=  '#984ea3'))+
          coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
          scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
        dia_plot
      }
      
    })
    
    # set output for liver plot
    output$liver_plot <- renderPlot({
      req(values$modelling_df_long, values$modelling_table)
      
      if(input$model_selection != 'All'){
        values$liver_plot
      }else{
        # build liver plot
        rug_data <- values$modelling_table %>% filter(Region == 'Liver')
        liver_plot <- values$modelling_df_long %>% filter(region == "Liver") %>% 
          ggplot(aes(x = time))+
          geom_point(aes(y = raw, color = "Raw data"), size = 1, alpha = 0.7)+
          geom_line(aes(y = lognormal, color = "Lognormal"), size = 1, alpha = 0.7)+
          geom_line(aes(y = gamma_variate, color = "Gamma Variate"), size = 1, alpha = 0.7)+
          geom_line(aes(y=ldrw, color = "LDRW"), size = 1, alpha = 1)+
          geom_line(aes(y = fpt, color = "FPT"), size = 1, alpha = 0.5)+
          geom_rug(data = rug_data, inherit.aes = F, 
                   aes(x =  MTT, linetype = "MTT", color = Model),
                   alpha = 0.5, length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F, 
                   aes(x = Peak, linetype = "t0 + Tp", color = Model), alpha = 0.5,
                   length = unit(1, 'npc'), size = 1)+
          geom_rug(data = rug_data, inherit.aes = F,
                   aes(x = t0, linetype = "t0", color = Model), alpha = 0.5, size = 1,
                   length = unit(0.3, 'npc'))+
          xlim(NA, values$max_time)+ labs(linetype = "Perfusion Parameter", color = "Model")+
          ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
          scale_color_manual(values = c('Raw data' = 'grey', 'Lognormal' = '#e41a1c', 
                                        'Gamma Variate' = '#377eb8',
                                        'LDRW' = '#4daf4a',
                                        'FPT'=  '#984ea3'))+
          coord_cartesian(clip = "off")+theme_classic()+theme(text = element_text(size = 18))+
          scale_linetype_manual(values = c('solid', 'dotdash', 'dotted'))
        liver_plot
      }
      
    })
    
    # build bfi plot
    output$bfi_plot <- renderPlotly({
      
        req(values$df, input$roi_selection, values$raw_data, values$filtered_data, 
            values$df_reg_plot)
        #, values$df_reg_plot[values$regression_data]
        df <- values$df
        
        reg_trace <- values$regression_trace
        fig <- plot_ly(data = df, type = 'scatter', mode = 'markers')  %>%
            add_trace(
                x = ~time,
                y = ~get(values$raw_data),
                marker = list(
                    color = 'rgb(135, 135, 135)'
                ),
                name = "Raw Data",
                opacity = 0.8
            ) %>%
            add_trace(
                x = ~time,
                y = ~get(values$filtered_data),
                mode = 'lines',
                line = list(
                    color = 'rgb(251, 0, 6)'
                ),
                name = "Filtered Data"
            ) %>% 
            layout(annotations = values$bfi_plot_annotations,
                   shapes = values$bfi_plot_shapes,
                   xaxis = list(title = "Time (s)",
                                showline = T,
                                ticklen = 0,
                                zeroline = F,
                                showgrid = F,
                                range = c(0,values$max_time)
                                ),
                   yaxis = list(title = 'Acoustic Intensity (dB)',
                                ticklen = 0,
                                zeroline = F,
                                showline = T,
                                showgrid = F),
                                font = list(size = 15))
        
        if(input$regression_selection == 'Standard'){
            fig <- fig%>% add_trace(
                    data = values$df_reg_plot,
                    x = ~time,
                    y = ~get(values$regression_data),
                    mode = 'lines',
                    line = list(color = toRGB('black'),
                                dash = 'dashdot'),
                    inherit = F,
                    name = "Linear regression"
                ) 
        }else{
            fig <- fig%>% add_trace(
                data = values$df_reg_plot,
                x = ~time,
                y = ~get(values$regression_data),
                mode = 'lines',
                line = list(color = toRGB('black')),
                inherit = F,
                name = "Segmented regression"
            ) %>% add_trace(
                data = values$seg_slope_extrapolation,
                x = ~time,
                y = ~y,
                mode = 'lines',
                line = list(color = toRGB('black'),
                            dash = 'dashdot'),
                inherit = F,
                showlegend = F
            )
        }
        fig
        
    })
    
    
    # create observer to build summary plots
    observeEvent(input$model_final,{
      req(values$modelling_df_long)
      model_dict <- list('Lognormal' = 'lognormal',
                         'Gamma Variate' = 'gamma_variate',
                         'LDRW' = 'ldrw',
                         'FPT' = 'fpt')
      
      model_selection <-  model_dict[[input$model_final]]
      plot_model_string <- paste0("Modelled data (",input$model_final,")")
      plot_df <- values$modelling_df_long %>% select(region, time, raw, filtered, eval(model_selection))
      
      regions <- c("Chest Wall", "Diaphragm", "Liver")
      rug_data <- values$modelling_table %>% filter(Model == input$model_final) %>% 
        select(Region, MTT, Peak)
      colnames(rug_data)[1] <- 'region'
      
      for (i in 1:length(regions)) {
        bfi_temp <- values$bfi_df[which(values$bfi_df$region == regions[i]),]
        
        # subtract t0 calculated from bfi analysis from time, subtract baseline value from other traces
        plot_df[which(plot_df$region == regions[i]),] <- plot_df[which(plot_df$region == regions[i]),] %>%
          mutate(time = time-bfi_temp$t0,
                 raw = raw-bfi_temp$baseline,
                 filtered = filtered-bfi_temp$baseline
          )
        # subtract baseline from modeled trace
        plot_df[which(plot_df$region == regions[i]),model_selection] <- plot_df[which(plot_df$region == regions[i]),model_selection]-bfi_temp$baseline
        
        # subtract t0 from Tp and MTT
        rug_data[which(rug_data$region == regions[i]),] <- rug_data[which(rug_data$region == regions[i]),] %>% 
          mutate(MTT = MTT - bfi_temp$t0, Peak = Peak - bfi_temp$t0)
      }
      
      values$summary_plot1 <- plot_df %>% ggplot(aes(x = time))+geom_point(aes(y = raw, color = region), size = 0.5, alpha = 0.8)+
        geom_line(aes(y = get(model_selection), color = region))+
        labs(color = "Region")+
        ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.minor = element_blank(), panel.background = element_blank())
      
      values$summary_plot <- plot_df %>% ggplot(aes(x = time))+geom_point(aes(y = raw, color = "Raw data"), size = 0.5)+
                             geom_line(aes(y = filtered, color = "Filtered data"))+
                             geom_line(aes(y = get(model_selection), color = plot_model_string))+
                             ylab("Acoustic Intensity (dB)")+xlab("Time (s)")+
                             facet_wrap(~region)+
                             scale_color_manual(values = c('black','red','grey'))+
                             theme(axis.line = element_line(colour = "black"),
                                   panel.grid.minor = element_blank(), panel.background = element_blank(),
                                   legend.title = element_blank())
    })
    
    # build output for summary plot 1
    output$summary_plot <- renderPlotly({
      req(values$summary_plot)
      ggplotly(values$summary_plot)
    })
    
    # build output for summary plot 2
    output$summary_plot1 <- renderPlotly({
      req(values$summary_plot1)
      ggplotly(values$summary_plot1)
    })
    
    # set output for summary table
    output$summary_table <- renderDT({
      req(values$bfi_df, values$modelling_table)
      model_selection <- input$model_final
      bfi <- values$bfi_df %>% select(region, bfi,t0)
      model <- values$modelling_table %>% filter(Model == eval(model_selection)) %>% select(Region, MTT, AUC, Tp, AUC_inverse, Volume)
      colnames(bfi)[1] <- 'Region'
      table_out <- inner_join(bfi,model, by = 'Region')
      table_out <- table_out %>% mutate(MTT_adj = round(MTT - t0,2),
                                        MTT = round(MTT,2),
                                        Tp = round(Tp,2),
                                        AUC = round(AUC,2),
                                        AUC_inverse = signif(AUC_inverse,4),
                                        Volume = signif(Volume,4)) %>%
        select(Region, bfi, MTT, MTT_adj, AUC, AUC_inverse, Volume, Tp)
      colnames(table_out) <- c('Region', 'BFI (dB/s)', 'MTT (s)','Adjusted MTT (s)','AUC', '1/AUC',
                              'MTT/AUC','Tp (s)')
      values$summary_table <- table_out
      table_out
      
    })
    
    output$exclusion_comment <- renderText({
        "The grey shaded region above shows the data that is excluded from the blood flow index (BFI) analysis. 
        All data after the timepoint at which the 'global' peak intensity is reached 
        (i.e., the peak intensity inclusive of both the shaded and unshaded regions) is excluded by default. 
        However, this selection can be manually altered via the slider on the left.
        Once the BFI analysis range is set, the dashed salmon colored vertical line represents the peak intensity
        of the data within that range (i.e., to the left of the grey shaded region)."})
}

# Run the application 
shinyApp(ui = ui, server = server)
