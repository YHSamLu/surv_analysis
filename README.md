# Survival Analysis
Survival analysis is widely used to show the risk of patients, and especially Cox proportion hazard (PH) model was used for over 50 years that was proposed by Sir David Cox. Therefore, I used R to construct the model for the hazard ratio (HR) with 95% confidence interval and its p-value as data frame using coxph and other function.

```{r}
library(tidyverse)
library(dplyr)
library(survey)
library(survival)
library(survminer)
library(finalfit)
library(stringi)
```

Import the packages, `tidyverse` and `dplyr` turn the data matrix (data frame) into the form which is avaliable to be analyzed. `survey`, `survival`, and `survminer` are including the all survival model and time-to-event function. 

```{r}
load("surv_functions.RData") # as RData
```


## Fine-Gray Hazard/Cox Proportion Hazard Model 

Cox model takes varaibles to estimate time-to-event outcome that is two value (status and time).  The formula is the same as linear regression, we code that into the variables of the function. Our function export the table with hazard ratio (95% CI) and its p-value. However, Fine-Gray hazard model add a competing risk, the time-to-event outcome used to be a binary event such as "Yes" (1) or "No" (0) for mortality. Like the code shown:

```{r}
#Example for cox model
Var <- colnames(DF)# variables

Cox_Model <- HR_Cox_table(DF = DF,
                         survival_analysis = "Surv(Followup, Re_OS)", #time-to-event
                         variables = Var,
                         variations = "uni", #univariate/multivaraite
                         weights = DF$weight
)

Cox_Model
```

```{r}
#Example for Fine-Gray model
DF <- DF%>%
  mutate(FG_CSS = factor(ifelse(Re_OS==0, 0,
                                ifelse(Re_CSS==1, 1, 2))))#Status as 0 No, 1 Yes, and 2 competing risk

Var_FG <- colnames(DF)

FG_Model <- FineGray_Model(survival_analysis = "Surv(Followup, FG_CSS)", #time-to-event
                                variables = Var_FG,
                                DF = DF,
                                variations = "uni",
                                weights = DF$weight)

FG_Model
```

### Result (Table)

| Variables|Levels     |  HR (95% CI)      | P     |
|----------|:----------|:-----------------:|:-----:|
| Sex      | Male      | 1                 |       |
|          | Female    | 2.67 (1.12, 4.33) | 0.004 |
| BMI      | Normal    | 1                 |       |
|          | Overwieght| 1.02 (0.21, 2.13) | 0.842 |
| Stone    | 0 No      | 1                 |       |
|          | 1 Yes     | 7.67 (4.52, 11.4) |<0.001 |

## Survival Curve

Survival Curve (Kaplan-Meier Plots) could be displayed very easy by using `ggsurvplot` such as:

```{r}
fit_CSS <- survfit(Surv(Followup, Re_CSS) ~ Stone, data = DF)

KM_CSS <- ggsurvplot(fit_CSS,
           data = DF,
           surv.median.line = "hv", # Add medians survival

           # Change legends: title & labels
           legend.title = "Stone",
           legend.labs = c("0 No", "1 Yes"),

           # Add p-value and tervals
           pval = TRUE,
           conf.int = TRUE,

           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           risk.table.title = "Number at Risk",

           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw(), # Change ggplot2 theme

           # 5-years survival
           xlim = c(0,60),
           # breaks as a year (12 monthes)
           break.x.by = 12,
           xlab = "Time from surgery in months",
           ylab = "Cancer-Specific survival probability",
           title = "Cancer-Specific Survival Curve"
)

KM_CSS
```

![](https://github.com/YHSamLu/stone_w.UTUC/blob/main/figures/CSS.wmf)

### Curve with weighted

```{r}
fitCSS<- survfit(Surv(Followup_OS.CSS, Re_CSS) ~ Stone, data = SurvCurve, weights = DF$weight)
# weights = ps_model1_overlap$overlap_weight
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
KM_CSS_overlap <- SurvCurve_overlap(fitCSS,
                  DF,
                  curve_title = "Cancer-Specific Survival Curve (Overlap weighted)",
                  group_name = "Stone",
                  group_level = c("0 No", "1 Yes"))

KM_CSS_overlap
```

![](https://github.com/YHSamLu/stone_w.UTUC/blob/main/figures/CSSoverlap.wmf)