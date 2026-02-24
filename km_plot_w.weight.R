library(cowplot)
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

#overlap weighted survival curve
SurvCurve_overlap <- function(fit_model,
                              WD,
                              curve_title = "",
                              group_name = "",
                              group_level = names(fit_model$strata)){


  # Extract survival times at which risk is calculated
  time_points <- seq(0, 60, by = 12)

  # Function to get risk table with fractional counts (weighted)
  get_risk_table <- function(fit, times) {
    # Extract number at risk at each time point
    risk_data <- summary(fit, times = times, extend = TRUE)

    # Create a data frame
    df <- data.frame(
      time = rep(risk_data$time, length(fit$strata)),
      strata = rep(names(fit$strata), each = length(times)),
      n.risk = as.numeric(risk_data$n.risk)
    )
    return(df)
  }

  names(fitOS$strata) <- group_level
  risk_df <- get_risk_table(fit_model, time_points)

  # Format n.risk with 1 decimal digit
  risk_df <- risk_df %>%
    mutate(n.risk = round(n.risk, 1))

  surv_plot <- ggsurvplot(fit_model,
                          data = WD,
                          surv.median.line = "hv",
                          legend.title = group_name,
                          legend.labs = group_level,
                          pval = TRUE,
                          conf.int = TRUE,
                          risk.table = FALSE,  # no risk table here
                          palette = c("#E7B800", "#2E9FDF"),
                          ggtheme = theme_bw(),
                          xlim = c(0, 60),
                          break.x.by = 12,
                          xlab = "Time in months",
                          title = curve_title
  )

  risk_table_plot <- ggplot(risk_df, aes(x = time, y = strata, label = n.risk, group = strata, color = strata)) +
    geom_text(size = 4) +
    scale_x_continuous(breaks = time_points, limits = c(0, 60)) +
    labs(x = "Time in months", y = "Weighted of Number at Risk") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 8)) +
    scale_color_manual(values = c("#E7B800", "#2E9FDF"))

  combined_plot <- plot_grid(
    surv_plot$plot,
    risk_table_plot,
    ncol = 1,
    rel_heights = c(3, 1)  # adjust heights as needed
  )

  combined_plot
}


#Example
fitOS <- survfit(Surv(Followup_OS.CSS, Re_OS) ~ Sex, data = SurvCurve, weights = SurvCurve$overlap_weight)
# weights = ps_model1_overlap$overlap_weight
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
KM_OS_overlap <- SurvCurve_overlap(fitOS,
                                   SurvCurve,
                                   curve_title = "Overall Survival Curve (Overlap weighted)",
                                   group_name = "Sex",
                                   group_level = c("1 Male", "2 Female"))

KM_OS_overlap
