### BONUS_example.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 21 2022 (11:20) 
## Version: 
## Last-Updated: jun 21 2022 (11:21) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(data.table)
library(survival)
library(prodlim)

n.T <- 2000
n.C <- 2000

max.time <- 60
HR.high <- 0.05 
HR.mid <- 0.05*0.5
HR.low <- 0.05*0.05
HR.high2low <- c(HR.high,0.75*HR.high+0.25*HR.low,0.5*HR.high+0.5*HR.low,0.25*HR.high+0.75*HR.low,HR.low)
time <- c(0, 12, 18, 24, 36)
time.early <- c(0, 12, 16, 20, 24)

set.seed(10)
dt.scenario1 <- simBuyseTest(n.C = n.C, n.T = n.T, 
                             argsTTE = list(scale.C = 1/HR.high, shape.C = 1, dist.C = "weibull",
                                            scale.T = 1/HR.mid, shape.T = 1, dist.T = "weibull",
                                            scale.censoring.C = 0, shape.censoring.C = max.time, dist.censoring.C = "uniform",
                                            scale.censoring.T = 0, shape.censoring.T = max.time, dist.censoring.T = "uniform"))

dt.scenario2 <- simBuyseTest(n.C = n.C, n.T = n.T, 
                             argsTTE = list(scale.C = 1/HR.high, shape.C = 1, dist.C = "weibull",
                                            scale.T = list(1/HR.high2low), shape.T = list(time), dist.T = "piecewiseExp",
                                            scale.censoring.C = 0, shape.censoring.C = max.time, dist.censoring.C = "uniform",
                                            scale.censoring.T = 0, shape.censoring.T = max.time, dist.censoring.T = "uniform"))

dt.scenario3 <- simBuyseTest(n.C = n.C, n.T = n.T, 
                             argsTTE = list(scale.C = list(1/HR.high2low), shape.C = list(time), dist.C = "piecewiseExp",
                                            scale.T = list(1/HR.high2low), shape.T = list(time.early), dist.T = "piecewiseExp",
                                            scale.censoring.C = 0, shape.censoring.C = max.time, dist.censoring.C = "uniform",
                                            scale.censoring.T = 0, shape.censoring.T = max.time, dist.censoring.T = "uniform"))

## pdf("Figures/visual-scenario-surv.pdf", width = 15)
par(mfrow = c(1,3))
plot(prodlim(Hist(eventtime,status)~treatment, data = dt.scenario1)); title("scenario 1")
plot(prodlim(Hist(eventtime,status)~treatment, data = dt.scenario2)); title("scenario 2")
plot(prodlim(Hist(eventtime,status)~treatment, data = dt.scenario3)); title("scenario 3")
## dev.off()

e.BT1 <- BuyseTest(treatment ~ tte(eventtime, status), data = dt.scenario1, method.inference = "none")
summary(e.BT1)
 ##  endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta
 ## eventtime      100        65.44          33.15          0     1.41 0.3229

e.BT2 <- BuyseTest(treatment ~ tte(eventtime, status), data = dt.scenario2, method.inference = "none")
summary(e.BT2)
##  endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta
## eventtime      100        50.89          47.37          0     1.74 0.0351

e.BT3 <- BuyseTest(treatment ~ tte(eventtime, status), data = dt.scenario3, method.inference = "none")
summary(e.BT3)
##  endpoint total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta
## eventtime      100        46.24          42.73          0    11.04 0.0351


##----------------------------------------------------------------------
### BONUS_example.R ends here
