---
title: "A Practical Introduction to Propensity Score Analysis using R"
author: "Ehsan Karim [http://ehsank.com/]"
date: "30 Sept 2020: SPPH, UBC"
output:
  beamer_presentation: 
    keep_md: yes
    highlight: tango
  ioslides_presentation:
    widescreen: yes
    smaller: yes
    keep_md: yes
    highlight: tango
    css: slides.css
  slidy_presentation: 
    widescreen: yes
    smaller: yes
    keep_md: yes
    highlight: tango
    css: slides.css
classoption: "aspectratio=169"    
---



\newcommand{\indep}{\perp \!\!\! \perp}

## About this event

- [TI Methods Speaker Series page](https://www.ti.ubc.ca/2020/08/06/sep-30-ti-methods-speaker-series-a-practical-introduction-to-propensity-score-analysis-using-r/): Therapeutics Initiative
    - Dr. Carole Lunny

- [SARGC page](https://ssc.ca/en/students-and-recent-graduates-committee-0): Students and Recent Graduates Committee (SARGC) coordinate activities for the Statistical Society of Canada (SSC)'s student and recent graduate members
    - Md. Erfanul Hoque
    - Janie Coulombe

## Assumptions of the webinar

- *Target audience*: 
  - Familiar with regression
  - Familiar with R
    - will explain some necessary package / functions / arguments
  - have no/minimal idea about propensity score

- *Topics covered*
  - Not covering any new research
  - Not covering statistical theory 
    - implementation being the goal here
  - Not attempting to reach any clinical conclusion

## Format of the webinar

- *Presentation format*
  - Rather informal
  - 1 hr vs. 2 hr
  - Q/A at 
    - 45 min and 
    - at the end

- *Webinar Materials*
  - All reproducible codes provided
    - [ehsanx.github.io/SARGC-TIMethods/](https://ehsanx.github.io/SARGC-TIMethods/)
  - Necessary references cited in respective slides

## Outline

- [1] Data and Regression 
  - (Diagnostics)
- [2] Exact matching 
  - (motivation)
- [3] Propensity score matching 
  - (4 steps)
- [4] Propensity score Reviews in different disease areas 
  - (brief)

## [1] Right Heart Catheterization (RHC) Dataset 

The [dataset](http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.html) that we will use today is from Connors et al. (1996).


\includegraphics[width=0.4\linewidth]{images/citeRHC} \includegraphics[width=0.4\linewidth]{images/rhcvars} 

### Notations

- Outcome `Death` ($Y$)
  - Death at any time up to 180 Days
- Treatment `swang1`  ($A$: Swan-Ganz catheter)
  - Whether or not a patient received a RHC
- [Covariate list](http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/Crhc.html): $L$ (`age`, `sex`, `race` , $\ldots$)
- Analysis strategy: matching RHC patients with non-RHC patients


## [1] Right Heart Catheterization (RHC) Dataset 

- `RHC` is helpful in guiding therapy decision
  - Helps determine the pressures within the heart
  - Popularly beleived that `RHC` is benefitial
  - Conducting RCT is hard (ethical reasons)
  - Benefit of `RHC` was not shown earlier (1996)

- SUPPORT data has 2 phases
  - phase 1: prospective observational study
  - phase 2: cluster RCT
  - Data in this study is combined
  

\includegraphics[width=0.2\linewidth]{images/RHC} 

## [1] Load data



```r
# Load the cleaned up data. 
# Reproducible codes:
# https://ehsanx.github.io/SARGC-TIMethods/ 
analytic.data <- readRDS("data/RHC.Rds")
# Data size and number of variables 
dim(analytic.data)
```

```
## [1] 4767   23
```

```r
# variable names
names(analytic.data)
```

```
##  [1] "age"               "sex"               "race"             
##  [4] "Disease.category"  "DNR.status"        "APACHE.III.score" 
##  [7] "Pr.2mo.survival"   "No.of.comorbidity" "DASI.2wk.prior"   
## [10] "Temperature"       "Heart.rate"        "Blood.pressure"   
## [13] "Respiratory.rate"  "WBC.count"         "PaO2.by.FIO2"     
## [16] "PaCO2"             "pH"                "Creatinine"       
## [19] "Albumin"           "GComa.Score"       "RHC"              
## [22] "Death"             "ID"
```

## [1] Inspecting data: Crude


```r
require(tableone)
# 2 x 2 table
tab0 <- CreateTableOne(vars = "RHC",
               data = analytic.data, 
               strata = "Death")
print(tab0, showAllLevels = TRUE)
```

```
##          Stratified by Death
##           level  No           Yes          p      test
##   n              2013         2754                    
##   RHC (%) No RHC 1315 (65.3)  1268 (46.0)  <0.001     
##           RHC     698 (34.7)  1486 (54.0)
```

## [1] Inspecting data: Some baseline variables


```r
baselinevars <- c("age","sex", "race")
# Table 1
tab1 <- CreateTableOne(vars = baselinevars,
               data = analytic.data, 
               strata = "Death", includeNA = TRUE, 
               test = TRUE, smd = FALSE)
print(tab1, showAllLevels = FALSE, smd = FALSE)
```

```
##                   Stratified by Death
##                    No           Yes          p      test
##   n                2013         2754                    
##   age (%)                                    <0.001     
##      [-Inf,50)      713 (35.4)   400 (14.5)             
##      [50,60)        351 (17.4)   452 (16.4)             
##      [60,70)        426 (21.2)   789 (28.6)             
##      [70,80)        382 (19.0)   750 (27.2)             
##      [80, Inf)      141 ( 7.0)   363 (13.2)             
##   sex = Female (%)  919 (45.7)   865 (31.4)  <0.001     
##   race (%)                                    0.004     
##      white         1554 (77.2)  2225 (80.8)             
##      black          332 (16.5)   403 (14.6)             
##      other          127 ( 6.3)   126 ( 4.6)
```

## [1] Crude regression


```r
# adjust the exposure variable (primary interest)
fit0 <- glm(I(Death=="Yes")~RHC,
            family=binomial, data = analytic.data)
require(Publish)
publish(fit0)
```

```
##  Variable  Units OddsRatio       CI.95 p-value 
##       RHC No RHC       Ref                     
##              RHC      2.21 [1.96;2.49]  <1e-04
```

## [1] Adjusted regression


```r
# adjust the exposure variable + demographics
fit1 <- glm(I(Death=="Yes")~RHC + age + sex + race,
            family=binomial, data = analytic.data)
publish(fit1)
```

```
##  Variable     Units OddsRatio       CI.95     p-value 
##       RHC    No RHC       Ref                         
##                 RHC      2.71 [2.38;3.08]     < 1e-04 
##       age [-Inf,50)       Ref                         
##             [50,60)      3.56 [3.01;4.20]     < 1e-04 
##             [60,70)      0.74 [0.64;0.87]   0.0001274 
##             [70,80)      1.33 [1.15;1.54]     < 1e-04 
##           [80, Inf)      1.06 [0.93;1.21]   0.3633526 
##       sex      Male       Ref                         
##              Female      0.49 [0.43;0.56]     < 1e-04 
##      race     white       Ref                         
##               black      1.10 [0.93;1.31]   0.2666157 
##               other      0.97 [0.73;1.28]   0.8068800
```

## [1] Why adjust?

### Continuous outcome ($Y$)

- treated group $A=1$ (`RHC`)
- control group $A=0$ (`no RHC`)

Treatment effect = $E[Y|A=1]$ vs. $E[Y|A=0]$

- Would only work if 2 groups are comparable / exchangeable / ignorable treatment assignment
- Randomization with enough sample size is one 

### Binary outcome ($Y$)

Treatment effect = $prob[Y = 1|A=1]$ vs. $prob[Y=1|A=0]$

## [1] Why adjust?

In absence of randomization, 

$E[Y|A=1] - E[Y|A=0]$

includes

- Treatment effect
- Systematic differences in 2 groups (‘confounding’)
  - Doctors may prescribe tx more to frail and older age patients. 
  - In here, $L$ = age is a confounder.

## [1] Why adjust?

In absence of randomization, if age is a known issue

### Causal effect for young

- $E[Y|A=1, L =$ `younger age`$]$ - $E[Y|A=0, L =$ `younger age`$]$

### Causal effect for old

- $E[Y|A=1, L =$ `older age`$]$ - $E[Y|A=0, L =$ `older age`$]$

Conditional exchangeability; only works if $L$ is measured

## [1] Why adjust?




\includegraphics[width=0.5\linewidth]{images/dag1} 

This was not a completely randomized data; some observational data was combined.

## [1] Adjusted regression (v2)


```r
# adjust the exposure variable + adjustment variables
baselinevars <- c("age","sex", "race","Disease.category", 
                  "DNR.status", "APACHE.III.score",
                  "Pr.2mo.survival","No.of.comorbidity",
                  "DASI.2wk.prior","Temperature",
                  "Heart.rate", "Blood.pressure",
                  "Respiratory.rate", "WBC.count",
                  "PaO2.by.FIO2","PaCO2","pH",
                  "Creatinine","Albumin","GComa.Score")
out.formula <- as.formula(paste("I(Death=='Yes')", "~", 
                               paste(baselinevars, 
                                     collapse = "+")))
out.formula
```

```
## I(Death == "Yes") ~ age + sex + race + Disease.category + DNR.status + 
##     APACHE.III.score + Pr.2mo.survival + No.of.comorbidity + 
##     DASI.2wk.prior + Temperature + Heart.rate + Blood.pressure + 
##     Respiratory.rate + WBC.count + PaO2.by.FIO2 + PaCO2 + pH + 
##     Creatinine + Albumin + GComa.Score
```

## [1] Adjusted regression (v2)


```r
fit2 <- glm(out.formula,
            family=binomial, data = analytic.data)
publish(fit2) 
```

```
##           Variable     Units OddsRatio       CI.95     p-value 
##                age [-Inf,50)       Ref                         
##                      [50,60)      1.54 [1.27;1.87]     < 1e-04 
##                      [60,70)      0.73 [0.62;0.87]   0.0002517 
##                      [70,80)      1.13 [0.97;1.32]   0.1142952 
##                    [80, Inf)      1.12 [0.97;1.29]   0.1117607 
##                sex      Male       Ref                         
##                       Female      0.48 [0.42;0.55]     < 1e-04 
##               race     white       Ref                         
##                        black      1.12 [0.93;1.36]   0.2314302 
##                        other      1.00 [0.74;1.35]   0.9881736 
##   Disease.category       ARF       Ref                         
##                          CHF      1.64 [1.24;2.16]   0.0004557 
##                         MOSF      1.06 [0.89;1.26]   0.4986705 
##                        Other      1.14 [0.93;1.41]   0.2070459 
##         DNR.status        No       Ref                         
##                          Yes      2.20 [1.67;2.89]     < 1e-04 
##   APACHE.III.score                1.01 [1.00;1.01]   0.0115427 
##    Pr.2mo.survival                0.01 [0.01;0.02]     < 1e-04 
##  No.of.comorbidity                1.20 [1.13;1.29]     < 1e-04 
##     DASI.2wk.prior                0.96 [0.94;0.97]     < 1e-04 
##        Temperature                0.94 [0.90;0.98]   0.0038857 
##         Heart.rate                1.00 [1.00;1.00]   0.0191355 
##     Blood.pressure                1.00 [1.00;1.00]   0.5000114 
##   Respiratory.rate                1.00 [0.99;1.00]   0.7569105 
##          WBC.count                1.01 [1.00;1.01]   0.0360557 
##       PaO2.by.FIO2                1.00 [1.00;1.00]   0.2173154 
##              PaCO2                0.99 [0.99;1.00]   0.0902647 
##                 pH                0.88 [0.39;1.94]   0.7433254 
##         Creatinine                1.03 [1.00;1.07]   0.0874115 
##            Albumin                0.98 [0.90;1.07]   0.6650240 
##        GComa.Score                1.00 [0.99;1.00]   0.0260942
```
## [1] Adjusted regression (v2)


```r
plot(fit2, which =1)
```


\includegraphics[width=0.5\linewidth]{slidePS_files/figure-beamer/reg3bc1-1} 

- curvilinear trends?
  - logistic regression IS curvilinear by nature
  
## [1] Adjusted regression (v2)


```r
plot(fit2, which =3)
```


\includegraphics[width=0.5\linewidth]{slidePS_files/figure-beamer/reg3bc3-1} 

- red line is approximately horizontal?
- points have approximately equal spread around the red line? 
  - more about detecting heteroscedasticity?
  
## [1] Adjusted regression (v2)


```r
plot(fit2, which =4)
```


\includegraphics[width=0.5\linewidth]{slidePS_files/figure-beamer/reg3bc3jr-1} 

-  Cook's D estimates the influence of data points

## [2] Alternate to Regression


How sure are you about the model-specification? 

- Interaction?
- Polynomial?

- Potential solution?
  - Exact Matching

## [2] Exact Matching: 2 variables


```r
var.comb <- do.call('paste0', 
                    analytic.data[, c('race', 'sex')])
length(table(var.comb))
```

```
## [1] 6
```

```r
table(var.comb)
```

```
## var.comb
## blackFemale   blackMale otherFemale   otherMale whiteFemale   whiteMale 
##         331         404         113         140        1340        2439
```

```r
table(analytic.data$RHC,var.comb)
```

```
##         var.comb
##          blackFemale blackMale otherFemale otherMale whiteFemale whiteMale
##   No RHC         161       239          50        61         667      1405
##   RHC            170       165          63        79         673      1034
```

## [2] Exact Matching: 2 variables


```r
require(MatchIt)
# exact match by sex and race
m.out = matchit (RHC=="RHC" ~ sex + race,
                 data = analytic.data, 
                 method = "exact") 
m.out$nn
```

```
##           Control Treated
## All          2583    2184
## Matched      2583    2184
## Unmatched       0       0
## Discarded       0       0
```

## [2] Exact Matching: 3 variables


```r
var.comb <- do.call('paste0', 
                    analytic.data[, c('race', 'sex', 'age')])
length(table(var.comb))
```

```
## [1] 30
```

```r
table(analytic.data$RHC,var.comb=="otherMale[80, Inf)")
```

```
##         
##          FALSE TRUE
##   No RHC  2580    3
##   RHC     2183    1
```

```r
table(analytic.data$RHC,var.comb=="otherFemale[80, Inf)")
```

```
##         
##          FALSE TRUE
##   No RHC  2581    2
##   RHC     2184    0
```

## [2] Exact Matching: 3 variables


```r
# exact match by age, sex and race
m.out = matchit (RHC=="RHC" ~ age + sex + race,
                 data = analytic.data, 
                 method = "exact") 
m.out$nn
```

```
##           Control Treated
## All          2583    2184
## Matched      2581    2184
## Unmatched       0       0
## Discarded       0       0
```

## [2] Exact Matching: 3 variables


```r
matched.data <- match.data(m.out)
dim(matched.data)
```

```
## [1] 4765   25
```

```r
nrow(analytic.data)-nrow(matched.data) # subjects deleted
```

```
## [1] 2
```

```r
# Not taking into account of matched sets
fit1m <- glm(I(Death=="Yes")~RHC, 
            family=binomial, data = matched.data)
publish(fit1m)
```

```
##  Variable  Units OddsRatio       CI.95 p-value 
##       RHC No RHC       Ref                     
##              RHC      2.21 [1.96;2.49]  <1e-04
```

## [2] Exact Matching: many categorical variables


```r
m.out = matchit (RHC=="RHC" ~ age + sex + race + 
                   Disease.category + DNR.status,
                 data = analytic.data, 
                 method = "exact") 
m.out$nn
```

```
##           Control Treated
## All          2583    2184
## Matched      2524    2150
## Unmatched       0       0
## Discarded       0       0
```

## [2] Exact Matching: many categorical variables


```r
matched.data <- match.data(m.out)
dim(matched.data)
```

```
## [1] 4674   25
```

```r
fit2m <- glm(I(Death=="Yes")~RHC,
            family=binomial, data = matched.data)
publish(fit2m)
```

```
##  Variable  Units OddsRatio       CI.95 p-value 
##       RHC No RHC       Ref                     
##              RHC      2.23 [1.98;2.51]  <1e-04
```

## [2] Exact Matching: including a continuous variable


```r
m.out = matchit (RHC=="RHC" ~ age + sex + race + 
                   Disease.category + DNR.status+ 
                   Heart.rate, # continuous
                 data = analytic.data, 
                 method = "exact") 
m.out$nn
```

```
##           Control Treated
## All          2583    2184
## Matched       929     947
## Unmatched       0       0
## Discarded       0       0
```

## [2] Exact Matching: including more continuous variables


```r
m.out = matchit (RHC=="RHC" ~ age + sex + race + 
                   Disease.category + DNR.status+ 
                   Heart.rate + Blood.pressure + 
                   Temperature, 
                 data = analytic.data, 
                 method = "exact") 
m.out$nn
```

```
##           Control Treated
## All          2583    2184
## Matched         3       3
## Unmatched       0       0
## Discarded       0       0
```

## [2] Exact Matching: including more continuous variables


```r
matched.data <- match.data(m.out)
dim(matched.data)
```

```
## [1]  6 25
```

```r
nrow(analytic.data)-nrow(matched.data) # subjects deleted
```

```
## [1] 4761
```

```r
fit3m <- glm(I(Death=="Yes")~RHC,
            family=binomial, data = matched.data)
publish(fit3m)
```

```
##  Variable  Units OddsRatio        CI.95 p-value 
##       RHC No RHC       Ref                      
##              RHC      1.00 [0.03;29.81]       1
```

## [3] Propensity Score 

### Defining Propensity score (PS)

  - Conditional Probability of getting treatment, given the observed covariates
  - Prob(treatment: `A` = 1 | baseline or pre-treatment covariates: `L`)
      - Prob(`RHC` = treated/RHC group | `age`, `sex`, `race`, etc.)
      - f(L) = Prob(A=1|L)


```r
baselinevars
```

```
##  [1] "age"               "sex"               "race"             
##  [4] "Disease.category"  "DNR.status"        "APACHE.III.score" 
##  [7] "Pr.2mo.survival"   "No.of.comorbidity" "DASI.2wk.prior"   
## [10] "Temperature"       "Heart.rate"        "Blood.pressure"   
## [13] "Respiratory.rate"  "WBC.count"         "PaO2.by.FIO2"     
## [16] "PaCO2"             "pH"                "Creatinine"       
## [19] "Albumin"           "GComa.Score"
```

## [3] Propensity Score 


\includegraphics[width=0.5\linewidth]{images/citePS} 

### Theoretical result

**Rosenbaum, Rubin (1983)** showed: 

- For potential outcomes $(Y^0, Y^1)$, if you have sufficient observed covariate list $L$ to reduce confounding (`strong ignoribility'):  $A$ being treatment assignment here: 
  - i.e., if $(Y^0, Y^1) \indep A | L$ (Note that is this NOT $Y \indep A | L$)

- then
  - $(Y^0, Y^1) \indep A | PS$ and 
  - $A \indep L | PS$

## [3] Propensity Score 

### Assumptions

- no unmeasured confounding
- positivity ($ 0 < PS < 1 $)
- well-defined treatment
- sufficient overlap
- model-specification

## [3] Propensity Score 


\includegraphics[width=0.5\linewidth]{images/psvar} \includegraphics[width=0.5\linewidth]{images/psdesign} 

- Observed covariates are used to fix design
- Which covariates should be selected:
  - known to be a confounder (causes of `Death` and `RHC`)
  - known to be a cause of the outcome (risk factors of `Death`)
  - avoid known instruments or noise variables: **SE suffers**
  - mediating factors should be avoided (total effect = goal) 

- Stepwise (p-value or criterion based) not recommended
  - depending on sample size, different values can get selected
  - may select variables highly associated with $A$
- Don't look at the outcome (`Death`) in your data to select covariates

## [3] Propensity Score

Many ways to use propensity scores (PS) in the analysis

- **PS matching** [our focus today: intuitive!]
- PS weighting
- PS stratification
- PS used as a covariate

## [3] Propensity Score Matching


\includegraphics[width=0.5\linewidth]{images/citeaustin} 

### Propensity score matching has 4 steps

- Stage 1: exposure modelling: $PS = Prob(A=1|L)$
- Stage 2: Match by $PS$
- Stage 2: Assess balance and overlap ($PS$ and $L$)
- Stage 4: outcome modelling: $Prob(Y=1|A=1)$


## [3] Propensity Score Matching


\includegraphics[width=0.5\linewidth]{images/citeaustin0} 

- Assessment of Balance in the whole data 
  - balance = similarity of the covariate distributions
  - $d$ or $SMD > 0.1$ can be considered as imbalance


\includegraphics[width=0.2\linewidth]{images/d1} 

\includegraphics[width=0.4\linewidth]{images/d2} 


```r
tab1e <- CreateTableOne(vars = baselinevars,
               data = analytic.data, strata = "RHC", 
               includeNA = TRUE, 
               test = FALSE, smd = TRUE)
```


## [3] Propensity Score Matching

```r
print(tab1e, smd = TRUE)
```

```
##                                Stratified by RHC
##                                 No RHC          RHC             SMD   
##   n                               2583            2184                
##   age (%)                                                        0.181
##      [-Inf,50)                     573 (22.2)      540 (24.7)         
##      [50,60)                       432 (16.7)      371 (17.0)         
##      [60,70)                       638 (24.7)      577 (26.4)         
##      [70,80)                       603 (23.3)      529 (24.2)         
##      [80, Inf)                     337 (13.0)      167 ( 7.6)         
##   sex = Female (%)                 878 (34.0)      906 (41.5)    0.155
##   race (%)                                                       0.098
##      white                        2072 (80.2)     1707 (78.2)         
##      black                         400 (15.5)      335 (15.3)         
##      other                         111 ( 4.3)      142 ( 6.5)         
##   Disease.category (%)                                           0.557
##      ARF                          1206 (46.7)      909 (41.6)         
##      CHF                           194 ( 7.5)      209 ( 9.6)         
##      MOSF                          515 (19.9)      858 (39.3)         
##      Other                         668 (25.9)      208 ( 9.5)         
##   DNR.status = Yes (%)             317 (12.3)      155 ( 7.1)    0.176
##   APACHE.III.score (mean (SD))   50.06 (18.28)   60.74 (20.27)   0.553
##   Pr.2mo.survival (mean (SD))     0.63 (0.18)     0.57 (0.20)    0.301
##   No.of.comorbidity (mean (SD))   1.51 (1.17)     1.48 (1.13)    0.022
##   DASI.2wk.prior (mean (SD))     20.72 (5.68)    20.70 (5.03)    0.003
##   Temperature (mean (SD))        37.66 (1.73)    37.59 (1.83)    0.035
##   Heart.rate (mean (SD))        112.29 (40.13)  118.93 (41.47)   0.163
##   Blood.pressure (mean (SD))     85.13 (38.54)   68.20 (34.24)   0.464
##   Respiratory.rate (mean (SD))   28.93 (13.69)   26.65 (14.17)   0.163
##   WBC.count (mean (SD))          15.29 (11.22)   16.27 (12.55)   0.082
##   PaO2.by.FIO2 (mean (SD))      236.55 (113.17) 192.43 (105.54)  0.403
##   PaCO2 (mean (SD))              40.13 (14.23)   36.79 (10.97)   0.263
##   pH (mean (SD))                  7.39 (0.11)     7.38 (0.11)    0.133
##   Creatinine (mean (SD))          1.94 (2.06)     2.47 (2.05)    0.260
##   Albumin (mean (SD))             3.18 (0.65)     2.98 (0.93)    0.251
##   GComa.Score (mean (SD))        20.53 (30.22)   18.97 (28.26)   0.053
```

## [3] Propensity Score Matching

### Step 1: PS estimation

Specify the propensity score model to estimate propensity scores, and fit the model


```r
ps.formula <- as.formula(paste("I(RHC == 'RHC')", "~", 
                paste(baselinevars, collapse = "+")))
ps.formula
```

```
## I(RHC == "RHC") ~ age + sex + race + Disease.category + DNR.status + 
##     APACHE.III.score + Pr.2mo.survival + No.of.comorbidity + 
##     DASI.2wk.prior + Temperature + Heart.rate + Blood.pressure + 
##     Respiratory.rate + WBC.count + PaO2.by.FIO2 + PaCO2 + pH + 
##     Creatinine + Albumin + GComa.Score
```
- Coef of PS model fit is not of concern 
- Model can be rich: to the extent that prediction is better
- But look for multi-collinearity issues
  - SE too high?


## [3] Propensity score Matching

While PS has balancing property, PS is unknown and needs to be estimated:


```r
# fit logistic regression to estimate propensity scores
PS.fit <- glm(ps.formula,family="binomial", 
              data=analytic.data)
# extract estimated propensity scores from the fit
analytic.data$PS <- predict(PS.fit, 
                            newdata = analytic.data, type="response")
```

- Other machine learning alternatives are possible to use instead of logistic regression.
  - tree based methods have better ability to detect non-linearity / non-additivity (model-specification aspect)
  - shrinkage methods - lasso / elastic net may better deal with multi-collinearity
  - ensemble learners / super learners were successfully used
  - shallow/deep learning!

## [3] Propensity score Matching

- Don't loose sight that better *balance* is the ultimate goal for propensity score
- Prediction of $A$ is just a means to that end (as true PS is unknown).
- May attract variables highly associated with $A$


\includegraphics[width=0.5\linewidth]{images/citesuper0} \includegraphics[width=0.5\linewidth]{images/citesuper} \includegraphics[width=0.5\linewidth]{images/psalt} \includegraphics[width=0.5\linewidth]{images/psml} 

## [3] Propensity score Matching

### Step 1

```r
# summarize propensity scores
summary(analytic.data$PS)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.009182 0.268112 0.454924 0.458150 0.640362 0.975476
```

```r
# summarize propensity scores by exposure group
tapply(analytic.data$PS, analytic.data$RHC, summary)
```

```
## $`No RHC`
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.009182 0.184909 0.330687 0.357838 0.504012 0.974095 
## 
## $RHC
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.05156 0.42874 0.59400 0.57679 0.74044 0.97548
```

## [3] Propensity Score Matching

### Step 2: PS matching


```r
logitPS <-  -log(1/analytic.data$PS - 1) 
# logit of the propensity score
.2*sd(logitPS) # suggested in the literature
```

```
## [1] 0.2382708
```

```r
0.1*sd(logitPS) # we are using this
```

```
## [1] 0.1191354
```

```r
# choosing too strict PS has unintended consequences 
```

\includegraphics[width=0.5\linewidth]{images/citecapiler} \includegraphics[width=0.5\linewidth]{images/pscal} 


## [3] Propensity Score Matching

### Step 2

Match using estimates propensity scores

- nearest-neighbor (NN) matching
- without replacement
- with caliper = .1*SD of logit of propensity score
- with 1:1 ratio (pair-matching)


\includegraphics[width=0.3\linewidth]{images/nn} 

## [3] Propensity Score Matching

### Step 2

Match using estimates propensity scores


```r
set.seed(123)
match.obj <- matchit(ps.formula, data = analytic.data,
                     distance = analytic.data$PS, 
                     method = "nearest", replace=FALSE,
                     caliper = .1*sd(logitPS), ratio = 1)
# see matchit function options here
# https://www.rdocumentation.org/packages/MatchIt/versions/1.0-1/topics/matchit
analytic.data$PS <- match.obj$distance
summary(match.obj$distance)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.009182 0.268112 0.454924 0.458150 0.640362 0.975476
```

## [3] Propensity Score Matching

### Step 2


```r
match.obj$nn
```

```
##           Control Treated
## All          2583    2184
## Matched      1519    1519
## Unmatched    1064     665
## Discarded       0       0
```

## [3] Propensity Score Matching

### Step 2

Step 1 and 2 can be done together by specifying `distance`


```r
match.obj <- matchit(ps.formula, data = analytic.data,
                     distance = 'logit', 
                     method = "nearest", 
                     replace=FALSE,
                     caliper = .1*sd(logitPS), 
                     ratio = 1)
analytic.data$PS <- match.obj$distance
summary(match.obj$distance)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.009182 0.268112 0.454924 0.458150 0.640362 0.975476
```
## [3] Propensity Score Matching

### Step 2: Taking a closer look at the matches


```r
# Ref: https://lists.gking.harvard.edu/pipermail/matchit/2013-October/000559.html
matches <- as.data.frame(match.obj$match.matrix)
colnames(matches)<-c("matched_unit")
matches$matched_unit<-as.numeric(
  as.character(matches$matched_unit))
matches$treated_unit<-as.numeric(rownames(matches))
matches.only<-matches[!is.na(matches$matched_unit),]
head(matches.only)
```

```
##    matched_unit treated_unit
## 5           438            5
## 10         2385           10
## 12         4177           12
## 13         4429           13
## 17         5228           17
## 22         1009           22
```

## [3] Propensity Score Matching

### Step 2: Taking a closer look at the matches (1st pair)


```r
analytic.data[analytic.data$ID %in% 
                as.numeric(matches.only[1,]),]
```

```
##           age  sex  race Disease.category DNR.status APACHE.III.score
## 5     [60,70) Male white             MOSF        Yes               72
## 438 [80, Inf) Male white              ARF        Yes               93
##     Pr.2mo.survival No.of.comorbidity DASI.2wk.prior Temperature Heart.rate
## 5        0.43699980                 0       21.05078    34.79688        125
## 438      0.01399994                 2       15.95312    34.89844          0
##     Blood.pressure Respiratory.rate WBC.count PaO2.by.FIO2 PaCO2       pH
## 5               65               27 29.699219     478.0000    17 7.229492
## 438             46                0  8.699219     138.0938    82 7.019531
##     Creatinine  Albumin GComa.Score    RHC Death  ID        PS
## 5     3.599609 3.500000          41    RHC   Yes   5 0.4213309
## 438   3.500000 2.799805         100 No RHC   Yes 438 0.4125774
```

## [3] Propensity Score Matching

### Step 2: Taking a closer look at the matches (2nd pair)


```r
analytic.data[analytic.data$ID %in% 
                as.numeric(matches.only[2,]),]
```

```
##            age    sex  race Disease.category DNR.status APACHE.III.score
## 10   [-Inf,50) Female white              ARF         No               48
## 2385   [70,80) Female white              ARF         No               36
##      Pr.2mo.survival No.of.comorbidity DASI.2wk.prior Temperature Heart.rate
## 10         0.6689453                 1       23.25781        38.5        141
## 2385       0.6219997                 1       18.35156        36.5        120
##      Blood.pressure Respiratory.rate WBC.count PaO2.by.FIO2 PaCO2       pH
## 10               73               40 20.597656      68.0000    30 7.349609
## 2385             67                9  3.199707     168.5625    34 7.429688
##      Creatinine Albumin GComa.Score    RHC Death   ID        PS
## 10     0.500000     2.5           0    RHC    No   10 0.6565211
## 2385   1.599854     3.5           0 No RHC    No 2385 0.6292985
```


## [3] Propensity score Matching

### Step 3: Assessing balance and overlap 

**Balance is more important than prediction**! 

- Criteria to assess success of step 2: PS estimation
  - better balance
  - better overlap [no extrapolation!]
  - PS = 0 or PS = 1 needs close inspection


```r
boxplot(PS ~ RHC=='RHC', data = analytic.data, 
        lwd = 2, ylab = 'PS')
stripchart(PS ~ RHC=='RHC', vertical = TRUE, 
           data = analytic.data, method = "jitter", 
           add = TRUE, pch = 20, col = 'blue')
```


\includegraphics[width=0.3\linewidth]{slidePS_files/figure-beamer/ps3-1} 

## [3] Propensity Score Matching

### Step 3

Vizualization


```r
plot(match.obj, type = "jitter")
```


\includegraphics[width=0.5\linewidth]{slidePS_files/figure-beamer/ps8-1} 

```
## [1] "To identify the units, use first mouse button; to stop, use second."
```

```
## integer(0)
```


## [3] Propensity Score Matching

### Step 3

Vizualization for assessing overlap issues


```r
plot(match.obj, type = "hist")
```


\includegraphics[width=0.5\linewidth]{slidePS_files/figure-beamer/ps9-1} 

## [3] Propensity Score Matching

### Step 3

Assessment of Balance: Better than regression diagnostics!


```r
matched.data <- match.data(match.obj)
tab1m <- CreateTableOne(vars = baselinevars,
               data = matched.data, strata = "RHC", 
               includeNA = TRUE, 
               test = TRUE, smd = TRUE)
```

Compare the similarity of baseline characteristics between treated and untreated subjects in a the propensity score-matched sample. 

- In this case, we will compare SMD < 0.1 or not. 
- In some literature, other generous values (0.25) are proposed.


\includegraphics[width=0.5\linewidth]{images/citeaustin0} \includegraphics[width=0.5\linewidth]{images/smdcut} 

## [3] Propensity Score Matching

### Step 3


```r
print(tab1m, showAllLevels = FALSE, smd = TRUE, test = FALSE) 
```

```
##                                Stratified by RHC
##                                 No RHC          RHC             SMD   
##   n                               1519            1519                
##   age (%)                                                        0.042
##      [-Inf,50)                     367 (24.2)      365 (24.0)         
##      [50,60)                       255 (16.8)      272 (17.9)         
##      [60,70)                       385 (25.3)      395 (26.0)         
##      [70,80)                       368 (24.2)      351 (23.1)         
##      [80, Inf)                     144 ( 9.5)      136 ( 9.0)         
##   sex = Female (%)                 560 (36.9)      572 (37.7)    0.016
##   race (%)                                                       0.017
##      white                        1212 (79.8)     1212 (79.8)         
##      black                         229 (15.1)      234 (15.4)         
##      other                          78 ( 5.1)       73 ( 4.8)         
##   Disease.category (%)                                           0.043
##      ARF                           709 (46.7)      689 (45.4)         
##      CHF                           172 (11.3)      162 (10.7)         
##      MOSF                          452 (29.8)      479 (31.5)         
##      Other                         186 (12.2)      189 (12.4)         
##   DNR.status = Yes (%)             123 ( 8.1)      126 ( 8.3)    0.007
##   APACHE.III.score (mean (SD))   55.09 (18.55)   56.41 (19.33)   0.070
##   Pr.2mo.survival (mean (SD))     0.61 (0.19)     0.59 (0.19)    0.059
##   No.of.comorbidity (mean (SD))   1.53 (1.18)     1.52 (1.16)    0.005
##   DASI.2wk.prior (mean (SD))     20.75 (5.62)    20.75 (5.12)    0.001
##   Temperature (mean (SD))        37.67 (1.85)    37.65 (1.70)    0.010
##   Heart.rate (mean (SD))        115.75 (41.06)  115.73 (39.83)   0.001
##   Blood.pressure (mean (SD))     75.41 (36.03)   73.74 (35.75)   0.047
##   Respiratory.rate (mean (SD))   28.25 (13.69)   27.74 (14.13)   0.036
##   WBC.count (mean (SD))          15.93 (11.99)   15.76 (12.15)   0.014
##   PaO2.by.FIO2 (mean (SD))      214.17 (110.19) 210.06 (108.71)  0.038
##   PaCO2 (mean (SD))              37.57 (10.98)   37.47 (11.47)   0.009
##   pH (mean (SD))                  7.39 (0.11)     7.39 (0.11)    0.029
##   Creatinine (mean (SD))          2.24 (2.41)     2.28 (1.88)    0.019
##   Albumin (mean (SD))             3.07 (0.68)     3.05 (0.97)    0.029
##   GComa.Score (mean (SD))        18.50 (28.53)   18.68 (28.31)   0.006
```


## [3] Propensity Score Matching

### Step 3 

Possible to get p-values to check balance: but strongly discouraged

- P-value based balance assessment can be influenced by sample size


\includegraphics[width=0.5\linewidth]{images/citeaustin} 


```r
print(tab1m, showAllLevels = FALSE, smd = FALSE, test = TRUE) 
```

```
##                                Stratified by RHC
##                                 No RHC          RHC             p      test
##   n                               1519            1519                     
##   age (%)                                                        0.859     
##      [-Inf,50)                     367 (24.2)      365 (24.0)              
##      [50,60)                       255 (16.8)      272 (17.9)              
##      [60,70)                       385 (25.3)      395 (26.0)              
##      [70,80)                       368 (24.2)      351 (23.1)              
##      [80, Inf)                     144 ( 9.5)      136 ( 9.0)              
##   sex = Female (%)                 560 (36.9)      572 (37.7)    0.680     
##   race (%)                                                       0.896     
##      white                        1212 (79.8)     1212 (79.8)              
##      black                         229 (15.1)      234 (15.4)              
##      other                          78 ( 5.1)       73 ( 4.8)              
##   Disease.category (%)                                           0.707     
##      ARF                           709 (46.7)      689 (45.4)              
##      CHF                           172 (11.3)      162 (10.7)              
##      MOSF                          452 (29.8)      479 (31.5)              
##      Other                         186 (12.2)      189 (12.4)              
##   DNR.status = Yes (%)             123 ( 8.1)      126 ( 8.3)    0.895     
##   APACHE.III.score (mean (SD))   55.09 (18.55)   56.41 (19.33)   0.055     
##   Pr.2mo.survival (mean (SD))     0.61 (0.19)     0.59 (0.19)    0.102     
##   No.of.comorbidity (mean (SD))   1.53 (1.18)     1.52 (1.16)    0.889     
##   DASI.2wk.prior (mean (SD))     20.75 (5.62)    20.75 (5.12)    0.985     
##   Temperature (mean (SD))        37.67 (1.85)    37.65 (1.70)    0.786     
##   Heart.rate (mean (SD))        115.75 (41.06)  115.73 (39.83)   0.989     
##   Blood.pressure (mean (SD))     75.41 (36.03)   73.74 (35.75)   0.200     
##   Respiratory.rate (mean (SD))   28.25 (13.69)   27.74 (14.13)   0.315     
##   WBC.count (mean (SD))          15.93 (11.99)   15.76 (12.15)   0.703     
##   PaO2.by.FIO2 (mean (SD))      214.17 (110.19) 210.06 (108.71)  0.301     
##   PaCO2 (mean (SD))              37.57 (10.98)   37.47 (11.47)   0.801     
##   pH (mean (SD))                  7.39 (0.11)     7.39 (0.11)    0.425     
##   Creatinine (mean (SD))          2.24 (2.41)     2.28 (1.88)    0.593     
##   Albumin (mean (SD))             3.07 (0.68)     3.05 (0.97)    0.425     
##   GComa.Score (mean (SD))        18.50 (28.53)   18.68 (28.31)   0.868
```

## [3] Propensity Score Matching

### Step 3

Assessment of balance in the matched data


```r
smd.res <- ExtractSmd(tab1m)
t(round(smd.res,2))
```

```
##         age  sex race Disease.category DNR.status APACHE.III.score
## 1 vs 2 0.04 0.02 0.02             0.04       0.01             0.07
##        Pr.2mo.survival No.of.comorbidity DASI.2wk.prior Temperature Heart.rate
## 1 vs 2            0.06              0.01              0        0.01          0
##        Blood.pressure Respiratory.rate WBC.count PaO2.by.FIO2 PaCO2   pH
## 1 vs 2           0.05             0.04      0.01         0.04  0.01 0.03
##        Creatinine Albumin GComa.Score
## 1 vs 2       0.02    0.03        0.01
```

## [3] Propensity Score Matching

### Step 3: Variance ratio

- Variance ratios $\sim$ 1 means: 
 - equal variances in groups 
 - group balance
 - could vary from 1/2 to 2
 - other cut-points are suggested as well (0.8 to 1.2)


\includegraphics[width=0.5\linewidth]{images/psbal} \includegraphics[width=0.5\linewidth]{images/vr} 

## [3] Propensity Score Matching

### Step 3: Variance ratio


```r
require(cobalt)
baltab.res <- bal.tab(x = match.obj, data = analytic.data, 
                      treat = analytic.data$RHC, 
                      disp.v.ratio = TRUE)
```

```
## Note: 's.d.denom' not specified; assuming pooled.
```

```r
baltab.res$Balance$V.Ratio.Adj
```

```
##  [1] 1.0990553        NA        NA        NA        NA        NA        NA
##  [8]        NA        NA        NA        NA        NA        NA        NA
## [15]        NA 1.0867497 0.9714495 0.9605864 0.8305596 0.8395535 0.9408913
## [22] 0.9841995 1.0655834 1.0262382 0.9733399 1.0919443 1.0916685 0.6100881
## [29] 2.0325397 0.9847091
```

## [3] Propensity Score Matching

### Step 4: Outcome modelling

- Some flexibility in choosing outcome model 
  - considered independent of exposure modelling
  - some propose double robust approach
  - adjusting imbalanced covariates only?

Estimate the effect of treatment on outcomes using propensity score-matched sample

```r
fit3 <- glm(I(Death=="Yes")~RHC,
            family=binomial, data = matched.data)
publish(fit3)
```

```
##  Variable  Units OddsRatio       CI.95 p-value 
##       RHC No RHC       Ref                     
##              RHC      2.08 [1.80;2.41]  <1e-04
```

## [3] Propensity Score Matching

### Step 4: Outcome modelling


```r
out.formula
```

```
## I(Death == "Yes") ~ age + sex + race + Disease.category + DNR.status + 
##     APACHE.III.score + Pr.2mo.survival + No.of.comorbidity + 
##     DASI.2wk.prior + Temperature + Heart.rate + Blood.pressure + 
##     Respiratory.rate + WBC.count + PaO2.by.FIO2 + PaCO2 + pH + 
##     Creatinine + Albumin + GComa.Score
```

```r
fit3b <- glm(out.formula,
            family=binomial, data = matched.data)
publish(fit3b)
```

```
##           Variable     Units OddsRatio       CI.95     p-value 
##                age [-Inf,50)       Ref                         
##                      [50,60)      1.62 [1.26;2.09]   0.0001639 
##                      [60,70)      0.76 [0.61;0.94]   0.0131863 
##                      [70,80)      1.14 [0.94;1.38]   0.1935493 
##                    [80, Inf)      1.26 [1.06;1.51]   0.0084480 
##                sex      Male       Ref                         
##                       Female      0.45 [0.38;0.54]     < 1e-04 
##               race     white       Ref                         
##                        black      1.13 [0.89;1.44]   0.3224165 
##                        other      0.87 [0.59;1.28]   0.4779060 
##   Disease.category       ARF       Ref                         
##                          CHF      1.47 [1.05;2.06]   0.0251017 
##                         MOSF      1.02 [0.83;1.26]   0.8354884 
##                        Other      1.28 [0.94;1.73]   0.1134749 
##         DNR.status        No       Ref                         
##                          Yes      2.51 [1.69;3.75]     < 1e-04 
##   APACHE.III.score                1.00 [1.00;1.01]   0.2170832 
##    Pr.2mo.survival                0.01 [0.00;0.01]     < 1e-04 
##  No.of.comorbidity                1.18 [1.09;1.28]     < 1e-04 
##     DASI.2wk.prior                0.95 [0.94;0.97]     < 1e-04 
##        Temperature                0.95 [0.90;1.00]   0.0594218 
##         Heart.rate                1.00 [1.00;1.00]   0.0802537 
##     Blood.pressure                1.00 [1.00;1.00]   0.3146312 
##   Respiratory.rate                1.00 [1.00;1.01]   0.3492897 
##          WBC.count                1.00 [0.99;1.01]   0.5707415 
##       PaO2.by.FIO2                1.00 [1.00;1.00]   0.0128897 
##              PaCO2                1.00 [0.99;1.01]   0.8801280 
##                 pH                0.92 [0.33;2.54]   0.8731477 
##         Creatinine                1.03 [0.98;1.08]   0.1897378 
##            Albumin                1.01 [0.91;1.13]   0.8175665 
##        GComa.Score                1.00 [0.99;1.00]   0.0123163
```


## [3] Propensity Score Matching

### Step 4: Other cosiderations for outcome model

The above analysis do not take matched pair into consideration while regressing. Literature proposes different strategies: 

- do not control for pairs / clusters 
  - use `glm` as is
- control for pairs / clusters
  - use `cluster` option or GEE or conditional logistic
- Bootstrap for matched pairfor WOR
  - may not be appropriate for WR


\includegraphics[width=0.5\linewidth]{images/boot} 

## [3] Propensity Score Matching

### Step 4

- The example compared `RHC` (a treated group; target) vs `No RHC` (untreated).
- Thc corresponding treatment effect estimate is known as
  - Average Treatment Effects on the Treated (ATT) 
- Other estimates from PS analysis are possible that compared the whole population
  - what if everyone treated vs. what if nobody was treated (ATE)

## [3] Propensity Score Matching

### Other matching algorithms

- Optimal
- genetic matching
- CEM
- variable ratio NN

## [3] Propensity Score Matching

- MatchIt
- Matching

Other useful packages

- cobalt
- twang

Outdated package

- nonrandom

## [4] Discipline-specific PS Systematic Reviews 

- Propensity score matching most popular
  - Cardiovascular / Infective endocarditis / Intensive care 
  - Critical care / anesthesiology / Sepsis / Psychology
  - Cancer / Multiple sclerosis 
- Not meta-analysis; but reviews of usage of PS methods in different disciplines

\includegraphics[width=0.3\linewidth]{images/r1} \includegraphics[width=0.3\linewidth]{images/r2} \includegraphics[width=0.3\linewidth]{images/r3} \includegraphics[width=0.3\linewidth]{images/r4} \includegraphics[width=0.3\linewidth]{images/r5} \includegraphics[width=0.3\linewidth]{images/r6} \includegraphics[width=0.3\linewidth]{images/r7} \includegraphics[width=0.3\linewidth]{images/r8} \includegraphics[width=0.3\linewidth]{images/r9} 

## [4] Discipline-specific PS Systematic Reviews 

### Reporting Guideline

- Be specific about population of interest
  - ATT vs. ATE
  - exclusion criteria
- Be specific about exposure
  - no multiple version of treatment
  - no interference
  - comparator
- Report clearly about missing data
  - how handled
- Why PS matching (or other approach) was selected?  
- Software

## [4] Discipline-specific PS Systematic Reviews 

### Reporting Guideline

- How variables selected
- Any important variables not measured
  - proxy
- Model selection
  - interaction or polynomials
  - logistic vs. machine learning
- Overlap vs. balance 
  - numeric and visual



\includegraphics[width=0.3\linewidth]{images/books} 

## [4] Discipline-specific PS Systematic Reviews 

### Reporting Guideline

- Reduction % of the matched data: main objection against this method!
- Residual imbalance
  - refit PS model
- Subgroup analysis
  - Refit within each group for matching
- Sensitivity analysis
  - unmeasured confounder / hdPS
  - any positivity issue? Deleting extremes has consequences!
    - ad-hoc methods: truncation / trimming: bias-variance trade-off


\includegraphics[width=0.5\linewidth]{images/sub} \includegraphics[width=0.5\linewidth]{images/hdps} 

## Further Reading


\includegraphics[width=0.2\linewidth]{images/book1} 

\includegraphics[width=0.5\linewidth]{images/book} 

Companion site: [study.sagepub.com/leite](https://study.sagepub.com/leite)

## Thank you!

<center>
<font size="50">[ehsank.com/workshops/](https://ehsank.com/workshops/)</font>
</center>
