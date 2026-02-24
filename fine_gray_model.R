library(tidyverse)
library(dplyr)
library(tableone)
library(data.table)
library(WeightIt)
library(survey)
library(PSweight)
library(survival)
library(survminer)
library(finalfit)
library(stringi)

FineGray_Model <- function(DF = DF,                      # complete data frame
                           survival_analysis = "Surv(Followup, Mortality)",  # time-to-event
                           variables = c(),
                           variations = "multi",             # multi/uni
                           HR1 = 1,
                           weights = NULL                    # optional weights vector
){

  if (!is.null(weights)) {
    # Fit Cox model with weights if provided
    if (variations == "multi") {
      VARs <- paste(variables, collapse = "+")
      # survival formula and fine-gray model formula
      FML <- paste(survival_analysis, "~", VARs)
      FML_FG <- paste("Surv(fgstart, fgstop, fgstatus) ~ ", VARs)

      #fine-gray hazard model data frame
      FG <- finegray(formula = as.formula(FML),
                     data = DF,
                     na.action= na.pass,
                     weights = weights)

      #Fine-Gray Hazard Model in Cox
      HR_multi <- coxph(as.formula(FML_FG), data=FG)%>%
        fit2df(estimate_suffix = "_0.95CI")

      HR_title <- finalfit(DF,
                           dependent = "Surv(Followup_OS.CSS, Re_CSS)",
                           explanatory = variables)%>%
        as.data.frame()%>%
        dplyr::select(1,2,5)

      HR_title[HR_title$`HR (multivariable)`!="-",3] <- HR_multi[,2]

      DF_HR <- HR_title

    } else if (variations == "uni") {
      HR_uni <- c()

      for(i in 1:length(variables)){
        FML <- paste(survival_analysis, "~", variables[i])
        FML_FG <- paste("Surv(fgstart, fgstop, fgstatus) ~ ", variables[i])

        #fine-gray hazard model data frame
        FG <- finegray(formula = as.formula(FML),
                       data = DF,
                       na.action= na.pass,
                       weights = weights)

        HR_uni <- rbind(HR_uni,
                        coxph(as.formula(FML_FG), data=FG)%>%
                          fit2df(estimate_suffix = "_0.95CI"))
      }
      HR_title <- finalfit(DF,
                           dependent = "Surv(Followup_OS.CSS, Re_CSS)",
                           explanatory = variables)%>%
        as.data.frame()%>%
        dplyr::select(1,2,4)

      HR_title[HR_title$`HR (univariable)`!="-",3] <- HR_uni[,2]

      DF_HR <- HR_title

    } else {
      stop("Error: wrong analysis method, please use 'multi' or 'uni' in variations")
    }
  } else {
    # Use original finalfit call if no weights provided
    DF_HR <- DF %>%
      finalfit(survival_analysis , variables) %>%
      as.data.frame()
  }

  # Add columns for HR and p-value formatted output
  DF_HR <- DF_HR %>%
    mutate(HR = NA, P = NA)

  if(variations == "multi"){
    lev <- nrow(DF_HR)

    for(i in 1:lev){
      if(DF_HR$`HR (multivariable)`[i]=="-"){
        DF_HR[i, "HR"] <- HR1
        DF_HR[i, "P"] <- NA
      }else{
        HR_info <- stri_remove_empty(strsplit(DF_HR$`HR (multivariable)`[i],
                                              split = "[-/(/,/)/p=/p/ ]",
                                              fixed = FALSE)[[1]])
        DF_HR[i, "HR"] <- paste(HR_info[1], " (", HR_info[2], ", ", HR_info[3], ")", sep = "")
        DF_HR[i, "P"] <- ifelse(HR_info[4]=="<0.001", "<0.001*",
                                ifelse(as.numeric(HR_info[4])<0.01, paste(HR_info[4], "*", sep=""),
                                       ifelse(as.numeric(HR_info[4])<0.05, paste(HR_info[4], "*", sep=""), HR_info[4])))
      }
    }
    return(DF_HR[,c(1,2,4,5)])
  } else if(variations == "uni"){
    lev <- nrow(DF_HR)

    for(i in 1:lev){
      if(DF_HR$`HR (univariable)`[i]=="-"){
        DF_HR[i, "HR"] <- HR1
        DF_HR[i, "P"] <- NA
      }else{
        HR_info <- stri_remove_empty(strsplit(DF_HR$`HR (univariable)`[i],
                                              split = "[-/(/,/)/p=/p/ ]",
                                              fixed = FALSE)[[1]])
        DF_HR[i, "HR"] <- paste(HR_info[1], " (", HR_info[2], ", ", HR_info[3], ")", sep = "")
        DF_HR[i, "P"] <- ifelse(HR_info[4]=="<0.001", "<0.001*",
                                ifelse(as.numeric(HR_info[4])<0.01, paste(HR_info[4], "*", sep=""),
                                       ifelse(as.numeric(HR_info[4])<0.05, paste(HR_info[4], "*", sep=""), HR_info[4])))
      }
    }
    return(DF_HR[,c(1,2,4,5)])
  } else {
    stop("Error: wrong analysis method, please use 'multi' or 'uni' in variations")
  }
}
