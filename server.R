library(shiny)
library(purrr)
library(stringr)
library(shinyjs)
library(shinyWidgets)
library(mixOmics)
library(shinyvalidate)
library(grid)
library(gridExtra) 
library(cowplot)

server <- function(input, output, session) { 
  
  load("predictive_model_fit.RData")
  source("fun_perf.R")
  
  iv <- InputValidator$new()
  iv$add_rule("slider_Age", sv_between(0, 100))
  iv$enable()
  
  observeEvent(input$slider_Age, {
    if(!iv$is_valid()) {
      showNotification("Age must be between 0 and 100.", type = "message")
      updateNumericInput(
        session = session,
        "slider_Age",
        value = 50,
        min = 0,
        max = 100,
        step = 1
      )
    }

  })

  observeEvent(input$signature_0, {
    
    if (!all(signature_0_mod %in% input$signature_0)) { 
      lapply(setdiff(signature_0_mod, input$signature_0), function(ss) disable(paste0("slider_", ss)))
    }
    
    if (!is.null(input$signature_0)) { # doesn't work if nothing is selected
      lapply(input$signature_0, function(ss) enable(paste0("slider_", ss)))
    }
    
  }, ignoreNULL = F)
  
  
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
  
  
  list_render_color_0 <- lapply(seq_along(signature_0_mod), # no colors for age and gender
                                function(ii) {
                                  
                                  marker <- signature_0_mod[ii]
                                  jj <- ii # vec_slider_signature_0[ii]
                                  
                                  color <-  reactive({tags$style(HTML(paste0(".js-irs-", jj-1, " .irs-single, .js-irs-", jj-1,
                                                                             " .irs-bar-edge, .js-irs-", jj-1, " .irs-bar {background: Grey}")))})
                                  
                                })
  names(list_render_color_0) <- paste0("color_", signature_0_mod)
  
  
  
  list_render_color_1 <- lapply(seq_along(signature_1_mod), # no colors for age and gender
                                function(ii) {
    
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
      
    }
    )
    
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

  tfc_0 = function(ss) {
    force(ss)
    list_render_color_0[[paste0("color_", ss)]]()
  }
  
  output$color_signature_0 <- renderUI({ 
    lapply(signature_0_mod, function(ss) {
      tfc_0(ss)
    })
  })
  
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
    
    input_signature_0 <- input$signature_0
    
    input_signature_1 <- input$signature_1
    input_signature_1[input$signature_1 %in% "Vg9Vd2 hi gd T"] <- "Vg9Vd2(hi) gd T"
    
    input_signature_2 <- input$signature_2
    input_signature_2[input$signature_2 %in% "IgG Memory"] <- "IgG+ Memory"
    
    input_signatures <- c(input_signature_0, input_signature_1, input_signature_2)
    input_signatures_mod <- c(input$signature_0, input$signature_1, input$signature_2)
    
    active_markers <- sapply(seq_along(input_signatures), function(ii) {
      
      ss <-  input_signatures[ii]
      ss_mod <- input_signatures_mod[ii]
      
      if (ss %in% c("Age", "Gender")) {
        input[[paste0("slider_", ss_mod)]]
      } else {
        perc_hc_and_all <- Reduce("cbind", list_perc_hc_and_all)
        perc_hc_and_all[paste0(input[[paste0("slider_", ss_mod)]], "%"), ss]
      }
      
    })
    names(active_markers) <-  input_signatures # input_signature_1
    
    inactive_markers <- setdiff(c(signature_0, signature_1, signature_2), 
                                names(active_markers))
    
    # /!\ convert back *_mod to * otherwise top_res_cov won't find "Vg9Vd2 hi gd T"
    list_X_new <- lapply(list_X_test_cov, function(ll) {
      
      ll <- ll[1,, drop = FALSE] # only one patient
      ll[!is.na(ll)] <- NA # initialise all entries to NA by precaution
      
      act_ll <- intersect(names(active_markers), names(ll))
      
      if (!is.null(act_ll)) {
        ll[match(act_ll, names(ll))] <- as.numeric(active_markers[act_ll])
      }
      
      ll*1 # to convert to numeric in case all NA otherwise converted to logical...
      
    })
    
    my_pred_new <- predict(top_res_cov, newdata = list_X_new)
    
    pred_classes <- as.vector(my_pred_new$AveragedPredict.class$max.dist)
    pred_scores <- as.vector(apply(my_pred_new$AveragedPredict, 3, max))
    names(pred_classes) <- names(pred_scores) <- c("comp1", "comp2") # comp2 scores is actually based on both comp1 and comp2

    cols <- c("sel-green", "sel-red")
    names(cols) <- c("i+ii", "iii")
    backgroundColour <<- as.character(cols[pred_classes[2]])

    backgroundColour2 <- ifelse(pred_scores[2] < 0.7, "sel2-orange", "sel2-white")
    
    output$sel <- renderText({
      paste0("<br/> <center> <b> Predicted systemic recovery class: ", 
             ifelse(pred_classes[2] == "iii", "3", "1+2"),
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
             "</b> from signatures 1, resp. 2) <br/><br/> ",
             "<center>  ROC curves per data type on test set <center/>") 
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
             "externally validated diagnostic tool but rather a pilot tool to guide ",
             "future clinically-actionable work.**</b> It is designed to provide an ",
             "estimate of the recovery prognosis for a SARS-CoV-2 infected patient when ", 
             "presented to the clinic at an early disease stage. The statistical ",
             "approach relies on an integrative prediction model for `systemic ", 
             "recovery' and on the availability of a selection of blood markers. <br/><br/>", #from the recovery signatures. 
              "The model outputs the estimated systemic recovery class ", #(using the ", 
             # "nomenclature in [REF]) 
             "along with a predicted score, which serves as a ",
             "measure of confidence in the predicted class. <br/><br/>",
             "The slider value of a given marker should be set based on where ", 
             "the patient levels stand with respect normal levels. It covers the ", 
             "percentiles of healthy controls and SARS-CoV-2 infected individuals enrolled at ",
             "Cambridge hospitals. The colors suggest abnormal/normal ranges: ",
             "the grey zone covers the healthy control interquartile range, ", 
             "the red & blue zones correspond to elevated, resp. reduced values of the ", 
             "marker with respect to the healthy control interquartile range. The ",
             "initial value is the healthy control median. <br/><br/>",
             "The practitioner may only have access to a subset of markers from ", 
             "signatures 1 and 2 for their patient. The drop-down menus allow ", 
             "filtering the markers. The predictive performance is re-assessed on the ",
             "left-out test set when restricting to the selected markers (the ROC ", 
             "curves are updated). In case of poor performance (e.g., AUC < 0.7), ",
             "the predicted class should be disregarded and input for additional ",
             "markers should be supplied where possible.<br/> ", 
             "&nbsp;")
    })
    
    list_X_test_cov_sub <- lapply(seq_along(list_X_test_cov), function(ii) {
      ll <- list_X_test_cov[[ii]]
      
      if (!is.null(inactive_markers)) {
        ll[, colnames(ll) %in% inactive_markers] <- NA
      }
      
      ll*1 # to convert to numeric in case all NA otherwise converted to logical...
    })
    names(list_X_test_cov_sub) <- names(list_X_test_cov)
    
    
    if (length(active_markers)>0) {
      set_null_device("cairo")
      res <- multi_auroc(top_res_cov, newdata = list_X_test_cov_sub, outcome.test = sub_df_info_test[, dep_var],
                  roc.comp = 2, bool_multiple = TRUE, bool_add_avg = TRUE, vec_col = vec_col_diablo)
  
      ff <- factor(res$leg, levels = res$leg)
      
      vec_all_col <- c(vec_col_diablo, "black")
      sub_names <- c("metabolites", "Cell", "Glyco", "Ratios", "gender", "Average")
      
      df <- data.frame(AUC = ff, empty = ff)
      p_empty <- ggplot(df, aes(empty, fill = AUC)) +
        geom_bar() + 
        scale_fill_manual(values = vec_all_col[sapply(sub_names, 
                                                      function(n_cc_i) 
                                                        {tt <- grep(n_cc_i, ff); ifelse(length(tt) == 0, F, T)})])
      legend <- get_legend(p_empty)

      grid.arrange(res$graph, legend, ncol = 2, nrow = 1, heights=4, widths=c(2,1))

    }

  })
  
  
}
