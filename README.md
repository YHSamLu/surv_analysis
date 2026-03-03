# Survival Analysis

Survival analysis is widely used to assess patient risk, with the Cox proportional hazards (PH) model—proposed by Sir David Cox over 50 years ago—being a cornerstone method. Here, I used R to fit the model via `coxph()`, generating a data frame of hazard ratios (HRs) with 95% confidence intervals and p-values.

```{r}
library(tidyverse)
library(dplyr)
library(survey)
library(survival)
library(survminer)
library(finalfit)
library(stringi)library(tableone)
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
load("surv_functions.RData")
```

Import key packages: `tidyverse` and `dplyr` for reshaping data frames into analysis-ready formats; `survey`, `survival`, and `survminer` for survival models and time-to-event analyses. Custom functions are also provided in `surv_functions.RData`.

## Patients Baseline Characteristics

Before analysis, create a baseline table to compare characteristics between cohorts.

```{r}
Baseline_Table(DF, Cohort)
```

- DF	a **data frame** or **matrix** with all variables, ensuring the grouping variable (defining cohorts for analysis) is the first column.
- Cohort	a **character**, name the grouping variable (defining cohorts for analysis) as the first column in your data frame.

### Example

```{r message=FALSE, warning=FALSE}
Baseline <- Baseline_Table(DF = WDmatrix, Cohort = "Sex")

Baseline%>% kbl(caption = "Patient's Baseline") %>%
  kable_paper("hover",full_width=F)
```

| Variables|Levels     |  1 Male  | 2 Female | P    | SMD|
|----------|:----------|:--------:|:--------:|:----:|:--:|
| BMI(%)   | Normal    |316 (55.9)|471 (61.0)| 0.07|0.103|
|          | Overwieght|249 (44.1)|301 (39.0)|      |    |
| Smoke(%) | 0 No      |334 (59.1)|747 (96.8)|**<0.001 \* **|1.02|
|          | 1 Yes     |231 (40.9)|25 (3.2)  |      |    |

p<0.05 *

## Fine-Gray Hazard/Cox Proportion Hazard Model 

The Cox model estimates time-to-event outcomes using covariates, with the outcome defined by two variables: status (event indicator) and time. The formula follows linear regression syntax, specifying variables directly in the function call. Our custom function exports a results table with hazard ratios (95% CI) and p-values.

The Fine-Gray model extends this by incorporating competing risks, where the binary event (e.g., mortality: 1=Yes, 0=No) accounts for alternative outcomes.

```{r}
HR_Cox_table(DF,
             survival_analysis,
             variables,
             variations,
             weights)
             
FineGray_Model(survival_analysis,
               variables,
               DF,
               variations,
               weights)
```

- DF	a **data frame** or **matrix** with all variables, ensuring the grouping variable (defining cohorts for analysis) is the first column.
- survival_analysis	a character string, typically `"Surv(time, event)"`, defines the time-to-event outcome.
- variables	a **vector** containing all predictor variables.
- variations	**"uni"** or **"multi"** representing univariate or multivariate.
- weights	a **vector** for each patient, with length matching the number of rows in the data frame.

### Example: univariate cox model for overall survival (OS)

```{r}
Var <- colnames(DF)

Cox_Model <- HR_Cox_table(DF = DF,
                         survival_analysis = "Surv(Followup, Re_OS)",
                         variables = Var,
                         variations = "uni",
                         weights = DF$weight
)

Cox_Model
```

### Example: multiivariate Fine-Gray hazard model for cancer-specific survival (CSS)

```{r}
DF <- DF%>%
  mutate(FG_CSS = factor(ifelse(Re_OS==0, 0,
                                ifelse(Re_CSS==1, 1, 2))))

Var_FG <- colnames(DF)

FG_Model <- FineGray_Model(survival_analysis = "Surv(Followup, FG_CSS)",
                                variables = Var_FG,
                                DF = DF,
                                variations = "multi",
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
|          | 1 Yes     | 7.67 (4.52, 11.4) |<0.001 \*|

p<0.05 *

## Survival Curve

Survival curves (Kaplan-Meier plots) can be generated easily using `ggsurvplot()` :

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
