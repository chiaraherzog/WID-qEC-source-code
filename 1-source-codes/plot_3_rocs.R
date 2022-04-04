plot_3_rocs <- function(type,
                          index1, index2, index3,
                          col1 = "black", col2 = "blue", col3 = "red",
                          title1 = "Index1", title2 = "Index2", title3 = "Index3",
                          direction = "<",
                        style = "default"){
  
  require(pROC)
  require(ggsci)
  col <- pal_npg(palette = c("nrc"), alpha = 1)(7)
  
  #--------------------------------#  
  # make annotation labels
  auc <- round(as.numeric(roc(type,index1, quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index1, quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index1, quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno1 <- paste('AUC (', title1, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type,index2, quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index2, quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index2, quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno2 <- paste('AUC (', title2, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type,index3, quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index3, quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index3, quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno3 <- paste('AUC (', title3, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "·", anno1)
    anno2 <- gsub("[.]", "·", anno2)
    anno3 <- gsub("[.]", "·", anno3)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type, index1, direction = direction)
  roc2 <- roc(type, index2, direction = direction)
  roc3 <- roc(type, index3, direction = direction)
  title1 <- as.character(title1)
  title2 <- as.character(title2)
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              colour = col2,
              size = 0.7) +
    geom_path(aes(x=1-roc3$specificities,
                  y=(roc3$sensitivities)),
              colour = col3,
              size = 0.7) +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_minimal() +
    theme(plot.title = element_text(size=10),
          panel.grid = element_blank()) +
    annotate(geom='label',
             x=0.55,
             y=0.4,
             label=anno1,
             fill=col1,
             colour='white',
             size=3.1)  +
    annotate(geom='label',
             x=0.55,
             y=0.25,
             label=anno2,
             fill=col2,
             colour='white',
             size=3.1) +
    annotate(geom='label',
             x=0.55,
             y=0.1,
             label=anno3,
             fill=col3,
             colour='white',
             size=3.1) 
  
}
