library(tableone)
library(kableExtra)
library(tidyverse)
library(data.table)
library(kableExtra)
library(MatchIt)
library(WeightIt)
library(MASS)
library(survey)
library(PSweight)
library(dplyr)

Baseline_Table <-function(DF, Cohort){
  vars_outcome0 <- colnames(DF)[-1]

  baseline_all0 <- CreateTableOne(data = DF,
                                  vars = vars_outcome0,
                                  strata = Cohort,
                                  includeNA = TRUE)

  a0 <- print(baseline_all0, smd = TRUE, showAllLevels = T)

  baseline <- a0%>%as.data.frame()%>%
    mutate(p = ifelse(p=="", "",
                      ifelse(p == "<0.001", "<0.001*",
                             ifelse(as.numeric(p)<0.05,
                                    paste(p, "*", sep = ""), p))))%>%
    mutate(Variables = rownames(a0), .before = 1)%>%
    dplyr::select(-test)

  for(i in 1:nrow(baseline)){
    baseline[i,3] <- gsub(c("\\( ", "\\(  "), "\\(", baseline[i,3])
    baseline[i,4] <- gsub(c("\\( ", "\\(  "), "\\(", baseline[i,4])
  }

  row.names(baseline) <- NULL

  return(baseline)
}


Baseline_Table(DF = WDmatrix, Cohort = "Sex")


