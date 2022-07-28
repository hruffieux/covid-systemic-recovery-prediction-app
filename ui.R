library(shiny)
library(purrr)
library(stringr)
library(shinyjs)
library(shinyWidgets)
library(BiocManager)
options(repos = BiocManager::repositories())

load("predictive_model_fit.RData")


ui <- fluidPage(
  tags$head(tags$style(
    HTML('
         #slider1 {
           // border: 1.4px solid rgba(0,0,0, 0.25);
            background-color: AliceBlue;
         }
         
         #slider1b {
           // border: 1.4px solid rgba(0,0,0, 0.25);
            background-color: AliceBlue;
         }
        
         #slider2 {
           // border: 1.1px solid rgba(0,0,0, 0.25);
            background-color: MintCream;
         }
        
        #slider2b {
           // border: 1.1px solid rgba(0,0,0, 0.25);
            background-color: MintCream;
        }

         #slider3 {
           // border: 1.1px solid rgba(0,0,0, 0.25);
            background-color: WhiteSmoke;
         }
         
         #slider3b {
           // border: 1.1px solid rgba(0,0,0, 0.25);
            background-color: WhiteSmoke;
         }
        
         #slider4 {
           // border: 1.4px solid rgba(0,0,0, 0.25);
            background-color: LightGray;
         }
        
        #slider4b {
           // border: 1.4px solid rgba(0,0,0, 0.25);
            background-color: LightGray;
        }
        
        #slider5 {
           // border: 1.4px solid rgba(0,0,0, 0.25);
            background-color: #f2e4ec;
         }

        body, label, input, button, select { 
          font-family: "Arial";
        }
        
         .sel-green{
           background-color: MediumSeaGreen;
         }
         
         .sel-red{
           background-color: DarkSalmon;
         }
         
          .sel2-orange{
           background-color: NavajoWhite;
         }
         
         .sel2-white{
           background-color: WhiteSmoke;
         }
         
         #textGuide{color: grey;
                                 font-size: 10.6px;
                                 //font-style: italic;
          }
         ')
  )),
  tabPanel("Predictive model", 
           titlePanel('Integrative prediction of systemic recovery from COVID-19'),
           p(paste0("Select patient markers from the two predictive signatures ", 
                    "(cell subsets, polar metabolites, ", 
                    "glyco- & lipo-proteins, diverse metabolic ratios)")),
           fluidRow(column(2, pickerInput(
             inputId = "signature_1",
             label = "Signature 1",
             choices = signature_1_mod,
             selected = signature_1_mod,
             options = list(`actions-box` = T, 
                            `none-selected-text` = "Please make a selection!",
                            `selected-text-format` = "count > 4",
                            selected = signature_1_mod[1]),
             multiple = TRUE
           )),
           column(2, pickerInput(
             inputId = "signature_2",
             label = "Signature 2",
             choices = signature_2_mod,
             selected = signature_2_mod,
             options = list(`actions-box` = T, 
                            # `none-selected-text` = "Please make a selection!",
                            `selected-text-format` = "count > 4", #,
                            selected = NULL #signature_2_mod[1]
             ),
             multiple = TRUE
           ))
           ),
           fluidRow(
             chooseSliderSkin("Flat"),
             shinyjs::useShinyjs(),
             column(2, offset = 0, id = "slider1", lapply(list_signature_1_mod[[1]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value = perc_corresp_hc_iqr_mod["50%", ss], # start at the median of HC levels
                           ticks = FALSE))),
             column(2, offset = 0, id = "slider1b", lapply(list_signature_2_mod[[1]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value = perc_corresp_hc_iqr_mod["50%", ss], 
                           ticks = FALSE)))
           ),
           fluidRow(
             chooseSliderSkin("Flat"),
             shinyjs::useShinyjs(),
             column(2, offset = 0, id = "slider2", lapply(list_signature_1_mod[[2]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE))),
             column(2, offset = 0, id = "slider2b", lapply(list_signature_2_mod[[2]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE)))
           ),
           fluidRow(
             chooseSliderSkin("Flat"),
             shinyjs::useShinyjs(),
             column(2, offset = 0, id = "slider3", lapply(list_signature_1_mod[[3]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE))),
             column(2, offset = 0, id = "slider3b", lapply(list_signature_2_mod[[3]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE)))
           ),
           fluidRow(
             chooseSliderSkin("Flat"),
             shinyjs::useShinyjs(),
             column(2, offset = 0, id = "slider4", lapply(list_signature_1_mod[[4]], function(ss) 
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE))),
             column(2, offset = 0, id = "slider4b", lapply(list_signature_2_mod[[4]], function(ss)
               sliderInput(paste0("slider_", ss),
                           str_c(ss, ':'),
                           min = 0, max = 100, step = 1, value =  perc_corresp_hc_iqr_mod["50%", ss],
                           ticks = FALSE)))
           ),
           fluidRow(
             chooseSliderSkin("Flat"),
             shinyjs::useShinyjs(),
             column(2, offset = 0, id = "slider5", lapply(list_signature_1_mod[[5]], function(ss) {
               if (ss == "Age") {
                 sliderInput(paste0("slider_", ss),
                             str_c(ss, ':'),
                             min = 0, max = 100, step = 1, value =  50, #perc_corresp_hc_iqr_mod["50%", ss],# <------------
                             ticks = FALSE)
               } else {
                 sliderInput(paste0("slider_", ss),
                             paste0(ss, ' F: 0, M: 1'),
                             min = 0, max = 1, step = 1, value =  1, #perc_corresp_hc_iqr_mod["50%", ss],# <------------
                             ticks = FALSE)
               }}))
           ),
           # mainPanel(),
           uiOutput("color_signature_1"),
           uiOutput("color_signature_2"),
           absolutePanel(id = "results", class = "panel panel-default", fixed = TRUE, 
                         style ="background-color: WhiteSmoke;
                         opacity: 0.85;
                         padding: 20px 20px 20px 20px;
                         margin: auto;
                         padding-bottom: 2mm;
                         padding-top: 1mm;",
                         # border-radius: 5pt;
                         # box-shadow: 0pt 0pt 6pt 0px rgba(61,59,61,0.48);
                         # style="padding-left: 10px; padding-right: 10px; padding-top: 0px; padding-bottom: 0px",
                         draggable = TRUE, top = 120, left = 505, right = "auto", bottom = "auto",
                         width = 600, #290, 
                         height = "auto",
                         htmlOutput("sel0"),
                         uiOutput("textBox", width = 10),
                         uiOutput("textBox2", width = 10),
                         htmlOutput("sel3"),
                         plotOutput("plot"),
                         br(),
                         cursor = "move"
           ),
           absolutePanel(id = "guide", class = "panel panel-default", fixed = TRUE, 
                         style ="background-color: WhiteSmoke;
                         opacity: 0.85;
                         padding: 20px 20px 20px 20px;
                         margin: auto;
                         padding-bottom: 2mm;
                         padding-top: 1mm;",
                         draggable = TRUE, top = 120, left = 1125, right = "auto", bottom = "auto",
                         width = 298, #290, 
                         height = "auto",
                         htmlOutput("textGuide"),
                         cursor = "move"
           )
  )
  
)