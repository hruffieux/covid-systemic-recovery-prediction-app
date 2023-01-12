multi_auroc <- function(
  object,
  newdata = object$X,
  outcome.test = as.factor(object$Y),
  multilevel = NULL,
  plot = TRUE,
  roc.block = 1,
  roc.comp = 1,
  my_predict = NULL,
  bool_multiple = FALSE,
  bool_add_avg = FALSE,
  bool_add_weight = FALSE,
  vec_col = NULL,
  ...)
{
  
  data=list()
  auc.mean = graph=list()
  data$outcome=factor(outcome.test)
  
  # note here: the dist does not matter as we used the predicted scores only
  
  list.predict  =  predict(object, newdata = newdata, dist = "max.dist", multilevel = multilevel)
  
  
  if (is.null(my_predict)) { # default function
    
    signatures <- c(Reduce("c", sapply(object$loadings, function(ll) { vec <- ll[, 1]; names(vec)[vec!=0]})[-length(object$loadings)]),
                    Reduce("c", sapply(object$loadings, function(ll) { vec <- ll[, 2]; names(vec)[vec!=0]})[-length(object$loadings)]))
    
    
    bool_keep <- sapply(newdata, function(nd) sum(!is.na(nd[, colnames(nd) %in% signatures]))>0) # only keep those data type for which newdata is not all NA 
    
    res.predict <- list.predict$predict[bool_keep]  
    
    if (is.null(vec_col)) {
      vec_col <- c("#999999", "#56B4E9", "#009E73", # "#F0E442", 
                   "#D55E00", "#CC79A7")[seq_along(res.predict)]
    } else {
      vec_col <- vec_col[bool_keep]
    }
    
    if (bool_multiple) {
      
      # controls the ordering of colors
      nn_pred <- names(res.predict)
      nn_pred[nn_pred == "MS"] <- "Polar metabolites"
      nn_pred[nn_pred == "cell_types"] <- "Cell subsets"
      nn_pred[nn_pred == "glyco-lipo-proteins"] <- "Glyco− and Lipo−proteins"
      nn_pred[nn_pred == "ratios"] <- "Ratios"
      nn_pred[nn_pred == "covariates"] <- "Age and gender"
      names(res.predict) <- paste0(" ", nn_pred) 
      names(vec_col) <- names(res.predict) 
      
      if (bool_add_avg) {
        res.predict <- append(res.predict, list(list.predict$AveragedPredict))
        names(res.predict)[length(res.predict)] <- " Average contributions"
        vec_col <- c(vec_col, "black")
        names(vec_col)[length(res.predict)] <- " Average contributions"
      } 
      
      if (bool_add_weight) {
        res.predict <- append(res.predict, list(list.predict$WeightedPredict))
        names(res.predict)[length(res.predict)] <- " Weighted contributions"
        vec_col <- c(vec_col, "black")
        names(vec_col)[length(res.predict)] <- " Weighted contributions"
      }
      
      data$data=lapply(res.predict, function(pp) pp[,,roc.comp])
      title="" 
      out = my_statauc_multiple(data, plot = TRUE, title = title, vec_col = vec_col)
      
    } else {
      
      block.all = names(res.predict)
      block.temp = names(res.predict[roc.block])
      
      for(j in 1:length(res.predict)) # j = data-type index
      {
        vec_col <- vec_col[j]
        for (i in 1:object$ncomp[j]) # ii = component id (can differ between data-types, that's why depends on j)
        {
          data$data=res.predict[[j]][,,i]
          title=paste("ROC Curve\nBlock: ", names(res.predict)[j], ", comp: ",i, sep="")
          
          plot.temp = ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp, plot, FALSE)
          temp = my_statauc(data, plot = plot.temp, title = title, vec_col = vec_col)
          auc.mean[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = temp[[1]]
          graph[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = temp$graph
          
        }
      }
      out = c(auc.mean,graph=graph)
      
    }
    
    
  } else {
    
    vec_col <- "black"
    
    if (my_predict == "Averaged") { # uses the average predicted values over the blocks
      res.predict <- list.predict$AveragedPredict
    } else if (my_predict == "Weighted") { # uses the weighted average of the predicted values over the blocks 
                                           # (using the predict and weights outputs)
      res.predict <- list.predict$WeightedPredict
    } else {
      stop("Predict type not implemented.")
    }
    
    data$data=res.predict[,,roc.comp]
    title=paste("ROC Curve\n", my_predict, " block contributions, comp: ",roc.comp, sep="")
    
    plot.temp = plot
    
    temp = my_statauc(data, plot = plot.temp, title = title, vec_col = vec_col)
    auc.mean[[paste("comp",roc.comp,sep = "")]] = temp[[1]]
    graph[[paste("comp",roc.comp,sep = "")]] = temp$graph
    
    out = c(auc.mean,graph=graph) 
  }
  
  return(invisible(out))
  
}

my_statauc <- function(data = NULL, plot = FALSE, title = NULL, vec_col = NULL){
  res.predict = data$data; outcome = data$outcome
  
  
  ann_text = matrix(ncol=2,nrow=nlevels(outcome))
  colnames(ann_text) = c("AUC", "p-value")
  
  df = NULL; seauc = NULL; zauc = NULL
  for (i in 1 : nlevels(outcome)){
    tempout = outcome
    levels(tempout)[-i] = "Other(s)"
    tempout = factor(tempout, c(levels(outcome)[i], "Other(s)"))
    temp = my_roc.default(response = tempout,
                          predictor = as.matrix(res.predict[, i],ncol=1))
    
    tryCatch({
      temp_plot <- pROC::roc(
        response = tempout, 
        predictor = res.predict[, i],
        smooth = TRUE
      )
    }, error=function(x){
      temp_plot <- pROC::roc(
        response = tempout, 
        predictor = res.predict[, i],
        smooth = FALSE
      )
      
    })
    
    seauc = sqrt((0.25 + (sum(table(tempout)) -2) * 0.083333)/(prod(table(tempout))))
    zauc = (temp$auc-0.5)/seauc
    zauc[zauc < 0] = - zauc[zauc < 0]
    for (k in unique(temp_plot$specificities)){
      temp_plot$sensitivities[which(temp_plot$specificities == k)] = temp_plot$sensitivities[rev(which(temp_plot$specificities == k))]
    }
    if (nlevels(outcome) == 2){
      ann_text = matrix(ncol=2,nrow=1)
      colnames(ann_text) = c("AUC", "p-value")
      df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                           paste(paste(levels(outcome)[1], levels(outcome)[2], sep = " vs "),"\n",
                                 "AUC: ", signif(temp$auc*100, 3), "%")))
      #ann_text[i , 1] =
      ann_text[i , 1] = signif(temp$auc, 3)
      ann_text[i , 2] = signif((1 - pnorm(zauc, 0, 1))*2, 3)
      break
    } else {
      df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                           paste(levels(outcome)[i], "vs Other(s)\n",
                                 "AUC: ", signif(temp$auc*100, 3), "%")))
      #ann_text[i , 1] =
      ann_text[i , 1] = as.numeric(signif(temp$auc, 3))
      ann_text[i , 2] = as.numeric(signif((1 - pnorm(zauc, 0, 1))*2, 3))
    }
  }
  #define rownames for ann_text
  if (nlevels(outcome) == 2){
    rownames(ann_text) = paste(levels(outcome)[1], levels(outcome)[2], sep = " vs ")
  }else{
    rownames(ann_text) = paste(levels(outcome), "vs Other(s)")
  }
  
  df = data.frame(df, stringsAsFactors = FALSE)
  names(df) = c("Specificity", "Sensitivity", "Outcome")
  df$Specificity = 100 - 100*as.numeric(df$Specificity)
  df$Sensitivity = 100*as.numeric(df$Sensitivity)
  
  if(plot)
  {
    
    if(is.null(title))
      title = "ROC Curve"
    else
      title=title
    p = ggplot(df, aes(x=Specificity,
                       y=Sensitivity,
                       group = Outcome,
                       colour = Outcome)) +
      xlab("100 - Specificity (%)") +
      ylab("Sensitivity (%)") +
      geom_line(linewidth = 0.75) +
      scale_x_continuous(breaks=seq(0, 100, by = 10)) +
      scale_y_continuous(breaks=seq(0, 100, by = 10)) +
      theme_bw() 
    
    p = p + geom_abline(intercept = 1) +
      theme(legend.key.size = unit(1, "cm"),
            plot.title = element_blank(), #element_text(lineheight=.8, face="bold"),
            legend.title = element_text(size=14, face="bold")) +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank()) + coord_fixed()
    
    if (is.null(vec_col)) {
      p = p + scale_color_brewer(palette="Spectral")  
    } else {
      p = p + scale_color_manual(values=vec_col) 
    }
    
  } else {
    p=NULL
  }
  return(list(ann_text,graph=p))
  
}




my_statauc_multiple <- function(data = NULL, plot = FALSE, title = NULL, vec_col = NULL){
  
  list_res.predict = data$data; outcome = data$outcome
  
  ann_text = matrix(ncol=2,nrow=nlevels(outcome))
  colnames(ann_text) = c("AUC", "p-value")
  
  list_df <- NULL
  
  for (dn in seq_along(list_res.predict)) {
    
    name_predict <- names(list_res.predict)[dn]
    res.predict <- list_res.predict[[dn]]
    
    df = NULL; seauc = NULL; zauc = NULL
    for (i in 1 : nlevels(outcome)){
      tempout = outcome
      levels(tempout)[-i] = "Other(s)"
      tempout = factor(tempout, c(levels(outcome)[i], "Other(s)"))
      temp = my_roc.default(response = tempout,
                            predictor = as.matrix(res.predict[, i],ncol=1))
      
      temp_plot <- pROC::roc(
        response = tempout, 
        predictor = res.predict[, i],
        smooth = F
      )
      
      
      tryCatch({
        temp_plot <- pROC::roc(
          response = tempout, 
          predictor = res.predict[, i],
          smooth = TRUE
        )
      }, error=function(x){
        temp_plot <- pROC::roc(
          response = tempout, 
          predictor = res.predict[, i],
          smooth = FALSE
        )
      })
      
      seauc = sqrt((0.25 + (sum(table(tempout)) -2) * 0.083333)/(prod(table(tempout))))
      zauc = (temp$auc-0.5)/seauc
      zauc[zauc < 0] = - zauc[zauc < 0]
      for (k in unique(temp_plot$specificities)){
        temp_plot$sensitivities[which(temp_plot$specificities == k)] = temp_plot$sensitivities[rev(which(temp_plot$specificities == k))]
      }
      if (nlevels(outcome) == 2){
        ann_text = matrix(ncol=2,nrow=1)
        colnames(ann_text) = c("AUC", "p-value")
        df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                             paste0(gsub("_", " ", name_predict), ": ", 
                                    signif(temp$auc*100, 3), "%")))
        ann_text[i , 1] = signif(temp$auc, 3)
        ann_text[i , 2] = signif((1 - pnorm(zauc, 0, 1))*2, 3)
        break
      } else {
        df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                             paste0(gsub("_", " ", name_predict), ": ", 
                                   signif(temp$auc*100, 3), "%")))
        ann_text[i , 1] = as.numeric(signif(temp$auc, 3))
        ann_text[i , 2] = as.numeric(signif((1 - pnorm(zauc, 0, 1))*2, 3))
      }
    }
    #define rownames for ann_text
    if (nlevels(outcome) == 2){
      rownames(ann_text) = paste(levels(outcome)[1], levels(outcome)[2], sep = " vs ")
    }else{
      rownames(ann_text) = paste(levels(outcome), "vs Other(s)")
    }
    
    df = data.frame(df, stringsAsFactors = FALSE)
    names(df) = c("Specificity", "Sensitivity", "Outcome")
    df$Specificity = 100 - 100*as.numeric(df$Specificity)
    df$Sensitivity = 100*as.numeric(df$Sensitivity)
    
    list_df <- append(list_df, list(df))

  }
  names(list_df) <- names(list_res.predict)
  
  if(plot)
  {
    if(is.null(title))
      title = "ROC Curve"
    else
      title=title
    p = ggplot(list_df[[1]], aes(x=Specificity,
                                 y=Sensitivity,
                                 group = Outcome,
                                 colour = "white")) +
      xlab("100 - Specificity (%)") +
      ylab("Sensitivity (%)") +
      geom_line(linewidth = 0.75) +
      scale_x_continuous(breaks=seq(0, 100, by = 10)) +
      scale_y_continuous(breaks=seq(0, 100, by = 10)) +
      theme_bw()  
    
    for (ii in 1:length(list_df)) {
 
      if (is.null(vec_col)) {
        p = p + geom_line(list_df[[ii]], mapping = aes(x=Specificity,
                                                       y=Sensitivity,
                                                       group = Outcome,
                                                       colour = Outcome), linewidth = 0.75) 
        
      } else {
        p = p + geom_line(list_df[[ii]], mapping = aes(x=Specificity,
                                                       y=Sensitivity,
                                                       group = Outcome,
                                                       colour = Outcome), linewidth = 0.75, col = vec_col[ii])
      }
    }
    if (is.null(vec_col)) {
      p = p + scale_color_brewer(palette = "Spectral") 
    } else {
      p = p + scale_colour_manual(values=vec_col, labels = unique(as.vector(sapply(list_df, "[[", "Outcome"))))
    }
    
    p = p + geom_abline(intercept = 1) +
      theme(legend.key.size = unit(1, "cm"),
            plot.title = element_text(lineheight=.8),
            legend.title = element_text(size=0),
            legend.text = element_text(size=10)) +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank()) + coord_fixed() 
    p <- p + guides(fill=guide_legend(title="Data type"))
    
    plot(p)
  } else {
    p=NULL
  }
  leg <- unique(as.vector(unlist(lapply(list_df, "[[", "Outcome"))))
  return(list(ann_text,graph=p, leg = leg))
  
}




my_roc.utils.perfs <- function(threshold, controls, cases, direction) {
  if (direction == '>') {
    tp <- sum(cases <= threshold)
    tn <- sum(controls > threshold)
  }
  else if (direction == '<') {
    tp <- sum(cases >= threshold)
    tn <- sum(controls < threshold)
  }
  # return(c(sp, se))
  return(c(sp=tn/length(controls), se=tp/length(cases)))
}

my_roc.default <- function(response, predictor,
                           
                           auc=TRUE,
                           
                           levels=base::levels(as.factor(response))
                           
                           
) {
  
  
  # Response / Predictor
  original.predictor <- predictor # store a copy of the original predictor (before converting ordered to numeric and removing NA)
  original.response <- response # store a copy of the original predictor (before converting ordered to numeric)
  
  # remove NAs if requested
  nas <- is.na(response) | is.na(predictor)
  if (any(nas)) {
    na.action <- grep(TRUE, nas)
    class(na.action) <- "omit"
    response <- response[!nas]
    attr(response, "na.action") <- na.action
    predictor <- predictor[!nas]
    attr(predictor, "na.action") <- na.action
  }
  
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(levels[1])]]
  cases <- splitted[[as.character(levels[2])]]
  
  # Remove patients not in levels
  patients.in.levels <- response %in% levels
  if (!all(patients.in.levels)) {
    response <- response[patients.in.levels]
    predictor <- predictor[patients.in.levels]
  }
  
  # update 13/01/17: first level as control to force directionality:
  # > : if the predictor values for the control group are higher than the values of
  #the case group (controls > t >= cases)
  
  #if (median(controls) <= median(cases))
  #direction <- "<"
  #else if (median(controls) > median(cases))
  direction <- ">"
  
  
  # create the roc object
  roc <- list()
  class(roc) <- "roc"
  
  # compute SE / SP
  thresholds <-((c(-Inf, sort(unique(c(controls, cases)))) + c(sort(unique(c(controls, cases))), +Inf))/2)
  perf.matrix <- sapply(thresholds, my_roc.utils.perfs, controls=controls, cases=cases, direction=direction)
  perfs <- list(se=perf.matrix[2,], sp=perf.matrix[1,])
  
  se <- perfs$se
  sp <- perfs$sp
  
  
  # store the computations in the roc object
  roc$sensitivities <- se
  roc$specificities <- sp
  roc$thresholds <- thresholds
  roc <- my_sort.roc(roc)
  roc$direction <- direction
  roc$cases <- cases
  roc$controls <- controls
  
  # compute AUC
  if (auc)
    roc$auc <- my_auc_roc(roc)
  
  roc$call <- match.call()
  roc$original.predictor <- original.predictor
  roc$original.response <- original.response
  roc$predictor <- predictor
  roc$response <- response
  roc$levels <- levels
  
  
  # return roc
  return(roc)
}

my_sort.roc <- function(roc) {
  if (is.unsorted(roc$specificities)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
    roc$thresholds <- rev(roc$thresholds)
  }
  return(roc)
}

my_auc_roc <- function(roc,
                       # Partial auc definition
                       partial.auc=FALSE, # false (consider total area) or numeric length 2: boundaries of the AUC to consider, between 0 and 1, or 0 and 100 if percent is TRUE
                       partial.auc.focus=c("specificity", "sensitivity"), # if partial.auc is not FALSE: do the boundaries
                       partial.auc.correct=FALSE,
                       allow.invalid.partial.auc.correct = FALSE,
                       ... # unused required to allow roc passing arguments to plot or ci.
) {
  if (!identical(partial.auc, FALSE)) {
    partial.auc.focus <- match.arg(partial.auc.focus)
  }
  
  percent <- FALSE
  
  # Validate partial.auc
  if (! identical(partial.auc, FALSE) & !(is.numeric(partial.auc) && length(partial.auc)==2))
    stop("partial.auc must be either FALSE or a numeric vector of length 2")
  
  # Ensure partial.auc is sorted with partial.auc[1] >= partial.auc[2]
  partial.auc <- sort(partial.auc, decreasing=TRUE)
  # Get and sort the sensitivities and specificities
  roc <- my_sort.roc(roc)
  se <- roc$se
  sp <- roc$sp
  
  # Full area if partial.auc is FALSE
  if (identical(partial.auc, FALSE)) {
    if (is(roc, "smooth.roc") && ! is.null(roc$smoothing.args) && roc$smoothing.args$method == "binormal") {
      coefs <- coefficients(roc$model)
      auc <- unname(pnorm(coefs[1] / sqrt(1+coefs[2]^2)) * ifelse(percent, 100^2, 1))
    }
    else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  }
  # Partial area
  else {
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }
    else {
      x <- sp
      y <- se
    }
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    }
    else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
  }
  
  # In percent, we have 100*100 = 10,000 as maximum area, so we need to divide by a factor 100
  if (percent)
    auc <- auc/100
  
  # Correction according to McClish DC, 1989
  if (all(!identical(partial.auc, FALSE), partial.auc.correct)) { # only for pAUC
    min <- roc.utils.min.partial.auc(partial.auc, percent)
    max <- roc.utils.max.partial.auc(partial.auc, percent)
    # The correction is defined only when auc >= min
    if (!allow.invalid.partial.auc.correct && auc < min) {
      warning("Partial AUC correction not defined for ROC curves below the diagonal.")
      auc <- NA
    }
    else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    }
    else {
      auc <- (1+((auc-min)/(max-min)))/2 # original formula by McClish
    }
  }
  # Prepare the AUC to return with attributes
  auc <- as.vector(auc) # remove potential pre-existing attributes
  attr(auc, "partial.auc") <- partial.auc
  attr(auc, "percent") <- percent
  attr(auc, "roc") <- roc
  if (!identical(partial.auc, FALSE)) {
    attr(auc, "partial.auc.focus") <- partial.auc.focus
    attr(auc, "partial.auc.correct") <- partial.auc.correct
  }
  class(auc) <- c("auc", class(auc))
  return(auc)
}
