### FIGURE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  9 2022 (09:53) 
## Version: 
## Last-Updated: nov 18 2022 (11:27) 
##           By: Brice Ozenne
##     Update #: 81
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(ggplot2)

if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    setwd("x:/GPC/article-restricted/")
}else{
    setwd("C:/Users/max/Desktop/simulation_peron/simulation-article-restricted")
}

## * Load results
dtS.sc1 <- readRDS(file = "Results/simSummary-ChemoVSChemo.rds")
dtS.sc2 <- readRDS(file = "Results/simSummary-ChemoVSImmuno.rds")
dtS.sc3 <- readRDS(file = "Results/simSummary-ImmunoVSImmuno.rds")
dtS.sc4 <- readRDS(file = "Results/simSummary-type1.rds")

relabel.estimator <- c("nbPeron" = "Net benefit",
                       "nbPeronTox1" = "Net benefit \n (equal toxicity)",
                       "nbPeronTox2" = "Net benefit \n (unequal toxicity)",
                       "rnbPeron24" = "Restricted net benefit \n (24 months)",
                       "rnbPeron36" = "Restricted net benefit \n (36 months)",
                       "rnbPeron48" = "Restricted net benefit \n (48 months)",
                       "logrank" = "Log-rank test",
                       "wlogrank" = "Weighted log-rank test",
                       "rmstDiff" = "Difference in RMST",
                       "rmstRatio" = "Ratio of RMST")
shape.estimator <- c("nbPeron" = 1,
                     "nbPeronTox1" = 2,
                     "nbPeronTox2" = 2,
                     "rnbPeron24" = 3,
                     "rnbPeron36" = 3,
                     "rnbPeron48" = 3,
                     "logrank" = 5,
                     "wlogrank" = 5,
                     "rmstDiff" = 6,
                     "rmstRatio" = 6)
color.estimator <- c("nbPeron" = "black",
                     "nbPeronTox1" = "red",
                     "nbPeronTox2" = "sienna",
                     "rnbPeron24" = "cyan",
                     "rnbPeron36" = "deepskyblue",
                     "rnbPeron48" = "blue",
                     "logrank" = "olivedrab3",
                     "wlogrank" = "forestgreen",
                     "rmstDiff" = "gold1",
                     "rmstRatio" = "darkorange")
relabel.scenario <- c("0" = "0",
                      "3" = "3 (th= 0.2 FU)",
                      "4" = "4 (oo FU)",
                      "5" = "5 (th = 0.3 FU)",
                      "6" = "6 (th = 0.4 FU)")

type.power <- "power5."
legend.power <- "Power (5% significance threshold)"
legend.type1 <- "Type 1 error (5% significance threshold)"

## * Generate plots
## ** ChemoVSChemo
dtEstimate.sc1 <- melt(dtS.sc1, id.vars = c("rep","censure","scenario","threshold","followUp"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc1[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc1[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtEstimate.sc1[,shape := factor(estimator,names(shape.estimator), shape.estimator)]
dtEstimate.sc1[,color := factor(estimator,names(color.estimator), color.estimator)]
dtEstimate.sc1[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc1 <- melt(dtS.sc1, id.vars = c("rep","censure","scenario","threshold","followUp"),
                    measure = patterns(type.power),
                    value.name = c("power"), variable.name = "estimator")
dtPower.sc1[,estimator := gsub(type.power,"",estimator)]
dtPower.sc1[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtPower.sc1[,shape := as.numeric(factor(estimator,names(shape.estimator), shape.estimator))]
dtPower.sc1[,color := factor(estimator,names(color.estimator), color.estimator)]
dtPower.sc1[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## dtPower.sc1[scenario==0 & threshold == 0 & followUp == 60]
##       rep   censure scenario threshold followUp                             estimator     power
##  1: 10000 0.0258685        0         0       60                         Log-rank test 0.9390000
##  2: 10000 0.0258685        0         0       60                Weighted log-rank test 0.8582157
##  3: 10000 0.0258685        0         0       60                    Difference in RMST 0.9387000
##  4: 10000 0.0258685        0         0       60                         Ratio of RMST 0.9393000
##  5: 10000 0.0258685        0         0       60                           Net benefit 0.8654000
##  6: 10000 0.0258685        0         0       60       Net benefit \n (equal toxicity) 0.8653000
##  7: 10000 0.0258685        0         0       60     Net benefit \n (unequal toxicity) 0.8659000
##  8: 10000 0.0258685        0         0       60 Restricted net benefit \n (24 months) 0.8458000
##  9: 10000 0.0258685        0         0       60 Restricted net benefit \n (36 months) 0.8622000
## 10: 10000 0.0258685        0         0       60 Restricted net benefit \n (48 months) 0.8653000

## petit changement de nom du graphe pour ne pas melanger
ggBenefit1 <- ggplot(dtEstimate.sc1[scenario == 0], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape = estimator))
ggBenefit1 <- ggBenefit1 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit1 <- ggBenefit1 + xlab("Follow-up time (months)") + ylab("Estimate")
ggBenefit1 <- ggBenefit1 + scale_shape_manual(name="",
                                              breaks = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggBenefit1 <- ggBenefit1 + scale_color_manual(name="",
                                              breaks = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])

ggBenefit1 <- ggBenefit1 + geom_hline(data = dtEstimate.sc1[scenario == 4][,scenario := 0], aes(yintercept = estimate, color = estimator))
ggBenefit1 <- ggBenefit1 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))

ggBenefit1.bis <- ggplot(dtEstimate.sc1[scenario %in% c(1:3,5:6)], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape=estimator))
ggBenefit1.bis <- ggBenefit1.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both, nrow = 1)
ggBenefit1.bis <- ggBenefit1.bis + xlab("Follow-up time (months)") + ylab("Estimate") + labs(color = "")
ggBenefit1.bis <- ggBenefit1.bis + scale_shape_manual(name="",
                                                      breaks = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggBenefit1.bis <- ggBenefit1.bis + scale_color_manual(name="",
                                                      breaks = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggBenefit1.bis <- ggBenefit1.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))


ggPower1 <- ggplot(dtPower.sc1[scenario==0], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower1 <- ggPower1 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggPower1 <- ggPower1 + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower1 <- ggPower1 + scale_shape_manual(name="",
                                          breaks = dtPower.sc1[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc1[scenario == 0][!duplicated(estimator),as.numeric(shape)])
ggPower1 <- ggPower1 + scale_color_manual(name="",
                                          breaks = dtPower.sc1[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc1[scenario == 0][!duplicated(estimator),as.character(color)])
ggPower1 <- ggPower1 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower1.bis <- ggplot(dtPower.sc1[scenario %in% c(1:3,5:6)], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower1.bis <- ggPower1.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both, nrow = 1)
ggPower1.bis <- ggPower1.bis + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower1.bis <- ggPower1.bis + scale_shape_manual(name="",
                                                  breaks = dtPower.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggPower1.bis <- ggPower1.bis + scale_color_manual(name="",
                                                  breaks = dtPower.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc1[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggPower1.bis <- ggPower1.bis + theme(text = element_text(size=15),
                                     axis.line = element_line(size = 1),
                                     axis.ticks = element_line(size = 1),
                                     axis.ticks.length=unit(.25, "cm"),
                                     legend.key.width = unit(3,"line"),
                                     legend.key.height = unit(2.5,"line"),
                                     legend.position = "bottom",
                                     panel.spacing = unit(1, "lines"))

## ** ChemoVSImmuno
dtEstimate.sc2 <- melt(dtS.sc2, id.vars = c("rep","censure","scenario","threshold","followUp"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc2[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc2[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtEstimate.sc2[,shape := factor(estimator,names(shape.estimator), shape.estimator)]
dtEstimate.sc2[,color := factor(estimator,names(color.estimator), color.estimator)]
dtEstimate.sc2[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc2 <- melt(dtS.sc2, id.vars = c("rep","censure","scenario","threshold","followUp"),
                 measure = patterns(type.power),
                 value.name = c("power"), variable.name = "estimator")
dtPower.sc2[,estimator := gsub(type.power,"",estimator)]
dtPower.sc2[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtPower.sc2[,shape := as.numeric(factor(estimator,names(shape.estimator), shape.estimator))]
dtPower.sc2[,color := factor(estimator,names(color.estimator), color.estimator)]
dtPower.sc2[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## petit changement de nom du graphe pour ne pas melanger
ggBenefit2 <- ggplot(dtEstimate.sc2[scenario==0], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape = estimator))
ggBenefit2 <- ggBenefit2 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit2 <- ggBenefit2 + xlab("Follow-up time (months)") + ylab("Estimate")
ggBenefit2 <- ggBenefit2 + scale_shape_manual(name="",
                                              breaks = dtEstimate.sc2[scenario == 0][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc2[scenario == 0][!duplicated(estimator),as.numeric(shape)])
ggBenefit2 <- ggBenefit2 + scale_color_manual(name="",
                                              breaks = dtEstimate.sc2[scenario == 0][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc2[scenario == 0][!duplicated(estimator),as.character(color)])
ggBenefit2 <- ggBenefit2 + geom_hline(data = dtEstimate.sc2[scenario == 4][,scenario := 0], aes(yintercept = estimate, color = estimator))
ggBenefit2 <- ggBenefit2 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))


ggBenefit2.bis <- ggplot(dtEstimate.sc2[scenario %in% c(1:3,5:6)], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape = estimator))
ggBenefit2.bis <- ggBenefit2.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both)
ggBenefit2.bis <- ggBenefit2.bis + xlab("Follow-up time (months)") + ylab("Estimate")
ggBenefit2.bis <- ggBenefit2.bis + scale_shape_manual(name="",
                                                      breaks = dtEstimate.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggBenefit2.bis <- ggBenefit2.bis + scale_color_manual(name="",
                                                      breaks = dtEstimate.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggBenefit2.bis <- ggBenefit2.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))

ggPower2 <- ggplot(dtPower.sc2[scenario==0], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower2 <- ggPower2 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggPower2 <- ggPower2 + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower2 <- ggPower2 + scale_shape_manual(name="",
                                          breaks = dtPower.sc2[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc2[scenario == 0][!duplicated(estimator),as.numeric(shape)])
ggPower2 <- ggPower2 + scale_color_manual(name="",
                                          breaks = dtPower.sc2[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc2[scenario == 0][!duplicated(estimator),as.character(color)])
ggPower2 <- ggPower2 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower2.bis <- ggplot(dtPower.sc2[scenario %in% c(1:3,5:6)], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower2.bis <- ggPower2.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both)
ggPower2.bis <- ggPower2.bis + scale_shape_manual(name="",
                                                  breaks = dtPower.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggPower2.bis <- ggPower2.bis + scale_color_manual(name="",
                                                  breaks = dtPower.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc2[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggPower2.bis <- ggPower2.bis + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower2.bis <- ggPower2.bis + theme(text = element_text(size=15),
                                     axis.line = element_line(size = 1),
                                     axis.ticks = element_line(size = 1),
                                     axis.ticks.length=unit(.25, "cm"),
                                     legend.key.width = unit(3,"line"),
                                     legend.key.height = unit(2.5,"line"),
                                     legend.position = "bottom",
                                     panel.spacing = unit(1, "lines"))

## ** ImmunoVSImmuno
dtEstimate.sc3 <- melt(dtS.sc3, id.vars = c("rep","censure","scenario","threshold","followUp"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc3[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc3[,shape := factor(estimator,names(shape.estimator), shape.estimator)]
dtEstimate.sc3[,color := factor(estimator,names(color.estimator), color.estimator)]
dtEstimate.sc3[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtEstimate.sc3[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc3 <- melt(dtS.sc3, id.vars = c("rep","censure","scenario","threshold","followUp"),
                 measure = patterns(type.power),
                 value.name = c("power"), variable.name = "estimator")
dtPower.sc3[,estimator := gsub(type.power,"",estimator)]
dtPower.sc3[,shape := as.numeric(factor(estimator,names(shape.estimator), shape.estimator))]
dtPower.sc3[,color := factor(estimator,names(color.estimator), color.estimator)]
dtPower.sc3[,scenario.f := factor(scenario,names(relabel.scenario), relabel.scenario)]
dtPower.sc3[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## petit changement de nom du graphe pour ne pas melanger
ggBenefit3 <- ggplot(dtEstimate.sc3[scenario==0], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape = estimator))
ggBenefit3 <- ggBenefit3 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit3 <- ggBenefit3 + xlab("Follow-up time (months)") + ylab("Estimate")
ggBenefit3 <- ggBenefit3 + scale_shape_manual(name="",
                                              breaks = dtEstimate.sc3[scenario == 0][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc3[scenario == 0][!duplicated(estimator),as.numeric(shape)])
ggBenefit3 <- ggBenefit3 + scale_color_manual(name="",
                                              breaks = dtEstimate.sc3[scenario == 0][!duplicated(estimator),estimator],
                                              values = dtEstimate.sc3[scenario == 0][!duplicated(estimator),as.character(color)])
ggBenefit3 <- ggBenefit3 + geom_hline(data = dtEstimate.sc3[scenario == 4][,scenario := 0], aes(yintercept = estimate, color = estimator))
ggBenefit3 <- ggBenefit3 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))

ggBenefit3.bis <- ggplot(dtEstimate.sc3[scenario %in% c(1:3,5:6)], aes(x = followUp, y = estimate, group = estimator, color = estimator, shape = estimator))
ggBenefit3.bis <- ggBenefit3.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both, nrow = 1)
ggBenefit3.bis <- ggBenefit3.bis + xlab("Follow-up time (months)") + ylab("Estimate")
ggBenefit3.bis <- ggBenefit3.bis + scale_shape_manual(name="",
                                                      breaks = dtEstimate.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggBenefit3.bis <- ggBenefit3.bis + scale_color_manual(name="",
                                                      breaks = dtEstimate.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                      values = dtEstimate.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggBenefit3.bis <- ggBenefit3.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))

ggPower3 <- ggplot(dtPower.sc3[scenario==0], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower3 <- ggPower3 + geom_point(size = 2) + geom_line(size = 1) + facet_grid(scenario~threshold, labeller = label_both)
ggPower3 <- ggPower3 + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower3 <- ggPower3 + scale_shape_manual(name="",
                                          breaks = dtPower.sc3[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc3[scenario == 0][!duplicated(estimator),as.numeric(shape)])
ggPower3 <- ggPower3 + scale_color_manual(name="",
                                          breaks = dtPower.sc3[scenario == 0][!duplicated(estimator),estimator],
                                          values = dtPower.sc3[scenario == 0][!duplicated(estimator),as.character(color)])
ggPower3 <- ggPower3 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower3.bis <- ggplot(dtPower.sc3[scenario %in% c(1:3,5:6)], aes(x = followUp, y = power, group = estimator, color = estimator, shape = estimator))
ggPower3.bis <- ggPower3.bis + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~scenario.f, labeller = label_both, nrow = 1)
ggPower3.bis <- ggPower3.bis + xlab("Follow-up time (months)") + ylab(legend.power)
ggPower3.bis <- ggPower3.bis + scale_shape_manual(name="",
                                                  breaks = dtPower.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.numeric(shape)])
ggPower3.bis <- ggPower3.bis + scale_color_manual(name="",
                                                  breaks = dtPower.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),estimator],
                                                  values = dtPower.sc3[scenario %in% c(1:3,5:6)][!duplicated(estimator),as.character(color)])
ggPower3.bis <- ggPower3.bis + theme(text = element_text(size=15),
                                     axis.line = element_line(size = 1),
                                     axis.ticks = element_line(size = 1),
                                     axis.ticks.length=unit(.25, "cm"),
                                     legend.key.width = unit(3,"line"),
                                     legend.key.height = unit(2.5,"line"),
                                     legend.position = "bottom",
                                     panel.spacing = unit(1, "lines"))

## ** scenario 4
range(dtS.sc4$power5.logrank)
range(dtS.sc4$power5.wlogrank)
range(dtS.sc4$power5.rmstDiff)
range(dtS.sc4$power5.nbPeron)
range(dtS.sc4$power5.rnbPeron)

ggType1 <- ggplot(dtS.sc4[scenario==0], aes_string(x = "rtime", y = paste0(type.power,"rnbPeron")))
ggType1 <- ggType1 + geom_abline(slope = 0, intercept = 0.05, color = "red")
ggType1 <- ggType1 + geom_point(size = 2) + geom_line(size = 1) + facet_wrap(~threshold, labeller = label_both)
ggType1 <- ggType1 + xlab("Follow-up time (months)") + ylab(legend.type1) + ggtitle("Restricted net benefit (Peron scoring rule)")
ggType1 <- ggType1 + theme(text = element_text(size=15),
                           axis.line = element_line(size = 1),
                           axis.ticks = element_line(size = 1),
                           axis.ticks.length=unit(.25, "cm"),
                           panel.spacing = unit(1, "lines"))

## * export
for(iExtension in c("pdf","png")){ ## iExtension <- "pdf"
    ggsave(ggBenefit1, filename = paste0("Figures/ggBenefit-ChemoVSChemo-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggBenefit1.bis, filename = paste0("Figures/ggBenefit-ChemoVSChemo-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggPower1, filename = paste0("Figures/ggPower-ChemoVSChemo-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggPower1.bis, filename = paste0("Figures/ggPower-ChemoVSChemo-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggBenefit2, filename = paste0("Figures/ggBenefit-ChemoVSImmuno-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggBenefit2.bis, filename = paste0("Figures/ggBenefit-ChemoVSImmuno-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggPower2, filename = paste0("Figures/ggPower-ChemoVSImmuno-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggPower2.bis, filename = paste0("Figures/ggPower-ChemoVSImmuno-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggBenefit3, filename = paste0("Figures/ggBenefit-ImmunoVSImmuno-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggBenefit3.bis, filename = paste0("Figures/ggBenefit-ImmunoVSImmuno-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggPower3, filename = paste0("Figures/ggPower-ImmunoVSImmuno-scenario0.",iExtension), height = 6, width = 12)
    ggsave(ggPower3.bis, filename = paste0("Figures/ggPower-ImmunoVSImmuno-scenario123.",iExtension), height = 6, width = 12)

    ggsave(ggType1, filename = paste0("Figures/ggType1.",iExtension), height = 6, width = 8)
}
##----------------------------------------------------------------------
### FIGURE.R ends here
