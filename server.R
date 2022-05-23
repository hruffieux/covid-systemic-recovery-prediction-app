library(shiny)
library(purrr)
library(stringr)
library(shinyjs)
library(shinyWidgets)
library(mixOmics)

server <- function(input, output, session) { 
  
  load("predictive_model_fit.RData")
  source("fun_perf.R")
  
  observeEvent(input$signature_1, {
    
    if(is.null(input$signature_1)){
      
      showNotification("At least one marker from signature 1 should be selected", 
                       duration = 4, type = "error")
      
      updatePickerInput(session = session,
                        inputId = "signature_1",
                        choices = signature_1_mod,
                        selected = signature_1_mod[1]
      )
      input$signature_1
    }
    
  },ignoreNULL = F)
  
  
  observeEvent(input$signature_1, {
    
    # print(input$signature_1)
    
    if (!all(signature_1_mod %in% input$signature_1)) {
      lapply(setdiff(signature_1_mod, input$signature_1), function(ss) disable(paste0("slider_", ss)))
    }
    lapply(input$signature_1, function(ss) enable(paste0("slider_", ss)))
    
  })
  
  
  observeEvent(input$signature_2, {
    
    if (!all(signature_2_mod %in% input$signature_2)) { 
      lapply(setdiff(signature_2_mod, input$signature_2), function(ss) disable(paste0("slider_", ss)))
    }
    
    if (!is.null(input$signature_2)) { # doesn't work if nothing is selected
      lapply(input$signature_2, function(ss) enable(paste0("slider_", ss)))
    }
    
  }, ignoreNULL = F)
  
  
  list_render_color_1 <- lapply(seq_along(signature_1_mod), function(ii) {
    
    marker <- signature_1_mod[ii]
    
    jj <- vec_slider_signature_1[ii]
    
    color <- reactive({
      
      if(input[[paste0("slider_", marker)]][1] <= perc_corresp_hc_iqr_mod["25%", marker]){
        tags$style(HTML(paste0(#".irs {max-width: 50px;}", 
          #" irs-bar {width: 100%; height: 15px; background: black; border-top: 1px solid black; border-bottom: 1px solid black;}",
          ".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Navy}"))) # CornflowerBlue
      }else if(input[[paste0("slider_", marker)]][1] <= perc_corresp_hc_iqr_mod["75%", marker]){
        tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Grey}")))
      }else{
        tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Crimson}"))) # Salmon
      }
      
    })
    
  })
  names(list_render_color_1) <- paste0("color_", signature_1_mod)
  
  
  list_render_color_2 <- lapply(seq_along(signature_2_mod), function(ii) {
    
    marker <- signature_2_mod[ii]
    
    jj <- vec_slider_signature_2[ii]
    
    color <- reactive({
      
      if(input[[paste0("slider_", marker)]][1] <= perc_corresp_hc_iqr_mod["25%", marker]){
        tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Navy}"))) 
      }else if(input[[paste0("slider_", marker)]][1] <= perc_corresp_hc_iqr_mod["75%", marker]){
        tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Grey}")))
      }else{
        tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1, 
          " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Crimson}"))) 
      }
      
    })
    
  })
  names(list_render_color_2) <- paste0("color_", signature_2_mod)
  
  tfc_1 = function(ss) {
    force(ss)
    list_render_color_1[[paste0("color_", ss)]]()
  }
  
  output$color_signature_1 <- renderUI({ 
    lapply(signature_1_mod, function(ss) {
      tfc_1(ss)
    })
  })
  
  tfc_2 = function(ss) {
    force(ss)
    list_render_color_2[[paste0("color_", ss)]]()
  }
  
  output$color_signature_2 <- renderUI({ 
    lapply(signature_2_mod, function(ss) {
      tfc_2(ss)
    })
  })
  
  
  output$plot <- renderPlot({
    
    input_signature_1 <- input$signature_1
    input_signature_1[input$signature_1 %in% "Vg9Vd2 hi gd T"] <- "Vg9Vd2(hi) gd T"
    
    input_signature_2 <- input$signature_2
    input_signature_2[input$signature_2 %in% "IgG Memory"] <- "IgG+ Memory"
    
    input_signatures <- c(input_signature_1, input_signature_2)
    input_signatures_mod <- c(input$signature_1, input$signature_2)
    
    active_markers <- sapply(seq_along(input_signatures), function(ii) {
      
      ss <-  input_signatures[ii]
      ss_mod <- input_signatures_mod[ii]
      
      perc_hc_and_all <- Reduce("cbind", list_perc_hc_and_all)
      perc_hc_and_all[paste0(input[[paste0("slider_", ss_mod)]], "%"), ss]
      
    })
    names(active_markers) <-  input_signatures # input_signature_1
    
    # print(active_markers)
    inactive_markers <- setdiff(c(signature_1, signature_2), 
                                names(active_markers))
    
    # /!\ convert back *_mod to * otherwise top_res won't find "Vg9Vd2 hi gd T"
    list_X_new <- lapply(list_X_test, function(ll) {
      
      ll <- ll[1,, drop = FALSE] # only one patient
      ll[!is.na(ll)] <- NA # initialise all entries to NA by precaution
      
      act_ll <- intersect(names(active_markers), names(ll))
      
      if (!is.null(act_ll)) {
        ll[match(act_ll, names(ll))] <- active_markers[act_ll]
      }
      
      ll*1 # to convert to numeric in case all NA otherwise converted to logical...
      
    })
    
    my_pred_new <- predict(top_res, newdata = list_X_new)
    
    pred_classes <- as.vector(my_pred_new$AveragedPredict.class$max.dist)
    pred_scores <- as.vector(apply(my_pred_new$AveragedPredict, 3, max))
    names(pred_classes) <- names(pred_scores) <- c("comp1", "comp2") # comp2 scores is actually based on both comp1 and comp2

    cols <- c("sel-green", "sel-red")
    names(cols) <- c("i+ii", "iii")
    backgroundColour <<- as.character(cols[pred_classes[2]])

    backgroundColour2 <- ifelse(pred_scores[2] < 0.7, "sel2-orange", "sel2-white")
    
    output$sel <- renderText({
      paste0("<br/> <center> <b> Predicted systemic recovery class: ", pred_classes[2],
             "<br/> ", 
             ifelse(pred_classes[2] == "iii",
                               "<i>Poor prognosis (pilot tool)</i>",
                               "<i>Good prognosis (pilot tool)</i>"), "</b> <br/> </center> ",
            paste0(rep(".", 143), collapse = ""))
    })
    
    output$sel0 <- renderText({
      paste0("<br/>")
    })
    
    output$sel2 <- renderText({
      paste0("<br/> &nbsp; <b>Predicted score: ", format(pred_scores[2], digits = 2), " </b>",
             "(score between 0.5 and 1, the higher the greater the confidence)",
             ifelse(pred_scores[2] < 0.7,
                    "<br/> &nbsp; <i>Warning: low score. Interpret predicted recovery class cautiously. </i> <br/> <br/>  ",
                    "<br/>"))
    })
    
    output$sel3 <- renderText({
      paste0("<br/> &nbsp;  Performance on test set using the selected markers (<b>", 
             length(input_signature_1),
             "</b> & <b>", 
             length(input_signature_2),
             "</b> from signatures 1, resp. 2) <br/><br/> ") 
    })
    
    output$textBox <- renderUI({
      htmlOutput("sel", class=backgroundColour)
    })
    
    output$textBox2 <- renderUI({
      htmlOutput("sel2", class=backgroundColour2)
    })
    
    output$textGuide <- renderText({
      paste0("<br/> <b>Quick Guide</b> <br/> <br/>",
             "<b>**The app has been developed using training-test splits of data ",
             "from a single cohort and therefore DOES NOT constitute an ",
             "externally validated diagnostic tool but rather a pilot to guide ",
             "future clinically-actionable work.**</b> The ambition is to provide an ",
             "estimate of the recovery prognosis for a COVID patient when ", 
             "presented to the clinic at an early disease stage. The statistical ",
             "approach relies on an integrative latent model for `systemic ", 
             "recovery' and on the availability of a selection of blood markers. <br/><br/>", #from the recovery signatures. 
              "The model outputs the estimated systemic recovery class ", #(using the ", 
             # "nomenclature in [REF]) 
             "along with a predicted score, serving as a ",
             "measure of confidence in the predicted class. <br/><br/>",
             "The slider value of a given marker should be set based on where ", 
             "the patient value stands with respect normal levels. It covers the ", 
             "percentiles of healthy controls and COVID-19 subjects enrolled at ",
             "Cambridge hospitals. The colors suggest abnormal/normal ranges: ",
             "the grey zone covers the healthy control interquartile range, ", 
             "the red & blue zones correspond to elevated, resp. reduced values of the ", 
             "marker with respect to the healthy control interquartile range. The ",
             "initial value is the healthy control median. <br/><br/>",
             "The practitioner may only have access to a subset of markers from ", 
             "signatures 1 and 2 for their patient. The drop-down lists allow ", 
             "filtering the markers. The performance on a left-out test set from ",
             "the Cambridge study is assessed by means of ROC curves when the ",
             "prediction is restricted to the subset of selected markers. The end ",
             "should be disregarded when the performance is insufficient (e.g. AUC <0.7).<br/> ", 
             "&nbsp;")
    })
    
    list_X_test_sub <- lapply(seq_along(list_X_test), function(ii) {
      ll <- list_X_test[[ii]]
      
      if (!is.null(inactive_markers)) {
        ll[, colnames(ll) %in% inactive_markers] <- NA
      }
      
      ll*1 # to convert to numeric in case all NA otherwise converted to logical...
    })
    names(list_X_test_sub) <- names(list_X_test)
    
    
    if (length(active_markers)>0) {
      multi_auroc(top_res, newdata = list_X_test_sub, outcome.test = sub_df_info_test[, dep_var],
                  roc.comp = 2, bool_multiple = TRUE, bool_add_avg = TRUE, vec_col = vec_col_diablo)
    }
  })
  
  
}
