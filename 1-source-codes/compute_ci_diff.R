# Author: Chiara Herzog
# Purpose: Compute CI of difference for sensitivity, specificity, PPV, NPV
# Note: Vars need to be coded as 0 or 1.
# Disease = disease status (1 diseased, 0 non-diseased.)
# Score = score (1 positive, 0 negative)
# t = t value, from T tables appropriate for sample size.
# pr = population prevalence. In this case always 9% for symptomatic populations but may be modified.
# Date: 18 May 2022

compute_ci_diff <- function(disease, score, variable, t, pr = 0.09){
  
  require(DTComPair)
  require(bdpv)
  
  # subfunctions
  var <- function(se, n){
    var <- (se*sqrt(n))^2
    return(var)
  }
  
  ci_diff <- function(x1, x2, var1, var2, t, n1, n2){
    pooled_var <- ((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2)
    cil <- (x2-x1) - t*sqrt((pooled_var)/n1 + (pooled_var)/n2)
    ciu <- (x2-x1) + t*sqrt((pooled_var)/n1 + (pooled_var)/n2)
    
    return(list(pooled_var = pooled_var,
                diff = x2-x1,
                cil = cil,
                ciu = ciu))
  }
  
  var1 <- as.character(unique(variable[variable != "N/A"]))[1]
  
  cat("Group 1 is ", var1, ", group 2 is ", as.character(unique(variable[variable != "N/A"]))[2],
      sep = "")
  
  # GROUP 1 Tabulation
  tab <- table(score[variable %in% c(var1, "N/A")], disease[variable %in% c(var1, "N/A")])
  tab <- matrix(c(tab[2,2], tab[1,2], tab[2,1], tab[1,1]), ncol = 2)
  colnames(tab) <- c("Case", "Control")
  rownames(tab) <- c("Case", "Control")
  group1 <- BDtest(xmat=as.matrix(tab), pr=pr, conf.level = 0.95)
  
  # GROUP 2 Tabulation
  tab <- table(score[variable != var1 | variable == "N/A"],
               disease[variable != var1 | variable == "N/A"])
  tab <- matrix(c(tab[2,2], tab[1,2], tab[2,1], tab[1,1]), ncol = 2)
  colnames(tab) <- c("Case", "Control")
  rownames(tab) <- c("Case", "Control")
  group2 <- BDtest(xmat=as.matrix(tab), pr=pr, conf.level = 0.95)
  
  # Sensitivity
  sens1 <- group1$SESPDAT$Estimate[1]
  sp1 <- group1$SESPDAT$Estimate[2]
  sens.se1 <- (group1$SESPDAT$`Upper 97.5% limit`[1]-group1$SESPDAT$`Lower 97.5% limit`[1])/(t*2)
  sp.se1 <- (group1$SESPDAT$`Upper 97.5% limit`[2]-group1$SESPDAT$`Lower 97.5% limit`[2])/(t*2)
  n1 <- sum(variable==var1)
  
  sens2 <- group2$SESPDAT$Estimate[1]
  sp2 <- group2$SESPDAT$Estimate[2]
  sens.se2 <- (group2$SESPDAT$`Upper 97.5% limit`[1]-group2$SESPDAT$`Lower 97.5% limit`[1])/(t*2)
  sp.se2 <- (group2$SESPDAT$`Upper 97.5% limit`[2]-group2$SESPDAT$`Lower 97.5% limit`[2])/(t*2)
  n2 <- sum(variable != var1)

  sens.v1 <- var(sens.se1, n1)
  sens.v2 <- var(sens.se2, n2)
  sp.v1 <- var(sp.se1, n1)
  sp.v2 <- var(sp.se2, n2)

  # Compute Differences
  sens.diff <- ci_diff(sens1, sens2, sens.v1, sens.v2, t, n1, n2)
  sens.diff <- as.data.frame(sens.diff)
  rownames(sens.diff) <- "sensitivity"
  
  sp.diff <- ci_diff(sp1, sp2, sp.v1, sp.v2, t, n1, n2)
  sp.diff <- as.data.frame(sp.diff)
  rownames(sp.diff) <- "specificity"
  
  
  #----- PPV/NPV
  
  # grp1
  ppv1 <- group1$PPVNPVDAT$Estimate[2]
  ppv1l <- group1$PPVNPVDAT$`Lower 97.5% limit`[2]
  ppv1h <- group1$PPVNPVDAT$`Upper 97.5% limit`[2]
  ppv.se1 <- (ppv1h-ppv1l)/(t*2)
  ppv.var1 <- var(ppv.se1, n1)
  
  npv1 <- group1$PPVNPVDAT$Estimate[1]
  npv1l <- group1$PPVNPVDAT$`Lower 97.5% limit`[1]
  npv1h <- group1$PPVNPVDAT$`Upper 97.5% limit`[1]
  npv.se1 <- (npv1h-npv1l)/(t*2)
  npv.var1 <- var(npv.se1, n1)

  # grp 2
  ppv2 <- group2$PPVNPVDAT$Estimate[2]
  ppv2l <- group2$PPVNPVDAT$`Lower 97.5% limit`[2]
  ppv2h <- group2$PPVNPVDAT$`Upper 97.5% limit`[2]
  ppv.se2 <- (ppv2h-ppv2l)/(t*2)
  ppv.var2 <- var(ppv.se2, n2)
  
  npv2 <- group2$PPVNPVDAT$Estimate[1]
  npv2l <- group2$PPVNPVDAT$`Lower 97.5% limit`[1]
  npv2h <- group2$PPVNPVDAT$`Upper 97.5% limit`[1]
  npv.se2 <- (npv2h-npv2l)/(t*2)
  npv.var2 <- var(npv.se2, n2)
  
  # compute difference
  ppv.diff <- ci_diff(ppv1, ppv2, ppv.var1, ppv.var2, t, n1, n2)
  ppv.diff <- as.data.frame(ppv.diff)
  rownames(ppv.diff) <- "ppv"
  
  npv.diff <- ci_diff(npv1, npv2, npv.var1, npv.var2, t, n1, n2)
  npv.diff <- as.data.frame(npv.diff)
  rownames(npv.diff) <- "npv"
  
  #----------
  diff <- rbind(sens.diff, sp.diff, ppv.diff, npv.diff)
  return(diff)
}
