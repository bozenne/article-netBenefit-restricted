### FIGURE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  9 2022 (09:53) 
## Version: 
## Last-Updated: sep 20 2022 (09:54) 
##           By: Brice Ozenne
##     Update #: 43
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

relabel.estimator <- c("logrank" = "Log-rank test",
                       "wlogrank" = "Weighted log-rank test",
                       "rmstDiff" = "Difference in RMST",
                       "rmstRatio" = "Ratio of RMST",
                       "nbPeron" = "Net benefit \n (Peron scoring rule)",
                       "nbPeronTox" = "Net benefit H2 \n (Peron scoring rule)",
                       "rnbPeron" = "Restricted net benefit \n (Peron scoring rule)")

type.power <- "power5."
legend.power <- "Power (5% significance threshold)"
legend.type1 <- "Type 1 error (5% significance threshold)"

## * Generate plots
## ** ChemoVSChemo
dtEstimate.sc1 <- melt(dtS.sc1, id.vars = c("rep","censure","scenario","threshold","rtime"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc1[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc1[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc1 <- melt(dtS.sc1, id.vars = c("rep","censure","scenario","threshold","rtime"),
                    measure = patterns(type.power),
                    value.name = c("power"), variable.name = "estimator")
dtPower.sc1[,estimator := gsub(type.power,"",estimator)]
dtPower.sc1[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## dtPower.sc1[scenario==0 & threshold == 0 & rtime == 60]
##      rep  censure scenario threshold rtime                                      estimator     power
## 1: 10000 0.025755        0         0    60                                  Log-rank test 0.9392000
## 2: 10000 0.025755        0         0    60                         Weighted log-rank test 0.8616617
## 3: 10000 0.025755        0         0    60                             Difference in RMST 0.9369000
## 4: 10000 0.025755        0         0    60                                  Ratio of RMST 0.9386000
## 5: 10000 0.025755        0         0    60            Net benefit \n (Gehan scoring rule) 0.8693000
## 6: 10000 0.025755        0         0    60 Restricted net benefit \n (Peron scoring rule) 0.8695000

## petit changement de nom du graphe pour ne pas melanger
ggBenefit1 <- ggplot(dtEstimate.sc1[scenario == 0], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit1 <- ggBenefit1 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit1 <- ggBenefit1 + xlab("Follow-up time (months)") + ylab("Estimate") + labs(color = "")
ggBenefit1 <- ggBenefit1 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))

ggBenefit1.bis <- ggplot( dtEstimate.sc1[scenario %in% c(1:3,5:6)], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit1.bis <- ggBenefit1.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both, nrow = 1)
ggBenefit1.bis <- ggBenefit1.bis + xlab("Follow-up time (months)") + ylab("Estimate") + labs(color = "")
ggBenefit1.bis <- ggBenefit1.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))


ggPower1 <- ggplot(dtPower.sc1[scenario==0], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower1 <- ggPower1 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggPower1 <- ggPower1 + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
ggPower1 <- ggPower1 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower1.bis <- ggplot(dtPower.sc1[scenario %in% c(1:3,5:6)], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower1.bis <- ggPower1.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both, nrow = 1)
ggPower1.bis <- ggPower1.bis + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
ggPower1.bis <- ggPower1.bis + theme(text = element_text(size=15),
                                     axis.line = element_line(size = 1),
                                     axis.ticks = element_line(size = 1),
                                     axis.ticks.length=unit(.25, "cm"),
                                     legend.key.width = unit(3,"line"),
                                     legend.key.height = unit(2.5,"line"),
                                     legend.position = "bottom",
                                     panel.spacing = unit(1, "lines"))

## ** ChemoVSImmuno
dtEstimate.sc2 <- melt(dtS.sc2, id.vars = c("rep","censure","scenario","threshold","rtime"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc2[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc2[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc2 <- melt(dtS.sc2, id.vars = c("rep","censure","scenario","threshold","rtime"),
                 measure = patterns(type.power),
                 value.name = c("power"), variable.name = "estimator")
dtPower.sc2[,estimator := gsub(type.power,"",estimator)]
dtPower.sc2[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## petit changement de nom du graphe pour ne pas melanger
ggBenefit2 <- ggplot(dtEstimate.sc2[scenario==0], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit2 <- ggBenefit2 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit2 <- ggBenefit2 + xlab("Follow-up time (months)") + ylab("Estimate") + ylab("") + labs(color = "")
ggBenefit2 <- ggBenefit2 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))


ggBenefit2.bis <- ggplot(dtEstimate.sc2[scenario %in% c(1:3,5:6)], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit2.bis <- ggBenefit2.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both)
ggBenefit2.bis <- ggBenefit2.bis + xlab("Follow-up time (months)") + ylab("Estimate") + ylab("") + labs(color = "")
ggBenefit2.bis <- ggBenefit2.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))

ggPower2 <- ggplot(dtPower.sc2[scenario==0], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower2 <- ggPower2 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggPower2 <- ggPower2 + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
ggPower2 <- ggPower2 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower2.bis <- ggplot(dtPower.sc2[scenario %in% c(1:3,5:6)], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower2.bis <- ggPower2.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both)
ggPower2.bis <- ggPower2.bis + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
ggPower2.bis <- ggPower2.bis + theme(text = element_text(size=15),
                                     axis.line = element_line(size = 1),
                                     axis.ticks = element_line(size = 1),
                                     axis.ticks.length=unit(.25, "cm"),
                                     legend.key.width = unit(3,"line"),
                                     legend.key.height = unit(2.5,"line"),
                                     legend.position = "bottom",
                                     panel.spacing = unit(1, "lines"))

## ** ImmunoVSImmuno
dtEstimate.sc3 <- melt(dtS.sc3, id.vars = c("rep","censure","scenario","threshold","rtime"),
                       measure = patterns("estimate."),
                       value.name = c("estimate"), variable.name = "estimator")
dtEstimate.sc3[,estimator := gsub("estimate.","",estimator)]
dtEstimate.sc3[,estimator := factor(estimator, levels = unique(estimator), labels = relabel.estimator[unique(estimator)])]

dtPower.sc3 <- melt(dtS.sc3, id.vars = c("rep","censure","scenario","threshold","rtime"),
                 measure = patterns(type.power),
                 value.name = c("power"), variable.name = "estimator")
dtPower.sc3[,estimator := gsub(type.power,"",estimator)]
dtPower.sc3[,estimator := factor(estimator, levels = names(relabel.estimator), labels = relabel.estimator)]

## petit changement de nom du graphe pour ne pas melanger
ggBenefit3 <- ggplot(dtEstimate.sc3[scenario==0], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit3 <- ggBenefit3 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggBenefit3 <- ggBenefit3 + xlab("Follow-up time (months)") + ylab("") + labs(color = "")
ggBenefit3 <- ggBenefit3 + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.key.width = unit(3,"line"),
                                 legend.key.height = unit(2.5,"line"),
                                 legend.position = "bottom",
                                 panel.spacing = unit(1, "lines"))

ggBenefit3.bis <- ggplot(dtEstimate.sc3[scenario %in% c(1:3,5:6)], aes(x = rtime, y = estimate, group = estimator, color = estimator))
ggBenefit3.bis <- ggBenefit3.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both, nrow = 1)
ggBenefit3.bis <- ggBenefit3.bis + xlab("Follow-up time (months)") + ylab("") + labs(color = "")
ggBenefit3.bis <- ggBenefit3.bis + theme(text = element_text(size=15),
                                         axis.line = element_line(size = 1),
                                         axis.ticks = element_line(size = 1),
                                         axis.ticks.length=unit(.25, "cm"),
                                         legend.key.width = unit(3,"line"),
                                         legend.key.height = unit(2.5,"line"),
                                         legend.position = "bottom",
                                         panel.spacing = unit(1, "lines"))

ggPower3 <- ggplot(dtPower.sc3[scenario==0], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower3 <- ggPower3 + geom_point() + geom_line() + facet_grid(scenario~threshold, labeller = label_both)
ggPower3 <- ggPower3 + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
ggPower3 <- ggPower3 + theme(text = element_text(size=15),
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.width = unit(3,"line"),
                             legend.key.height = unit(2.5,"line"),
                             legend.position = "bottom",
                             panel.spacing = unit(1, "lines"))


ggPower3.bis <- ggplot(dtPower.sc3[scenario %in% c(1:3,5:6)], aes(x = rtime, y = power, group = estimator, color = estimator))
ggPower3.bis <- ggPower3.bis + geom_point() + geom_line() + facet_wrap(~scenario, labeller = label_both, nrow = 1)
ggPower3.bis <- ggPower3.bis + xlab("Follow-up time (months)") + ylab(legend.power) + labs(color = "")
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
ggType1 <- ggType1 + geom_point() + geom_line() + facet_wrap(~threshold, labeller = label_both)
ggType1 <- ggType1 + xlab("Follow-up time (months)") + ylab(legend.type1) + ggtitle("Restricted net benefit (Peron scoring rule)")
ggType1 <- ggType1 + theme(text = element_text(size=15),
                           axis.line = element_line(size = 1),
                           axis.ticks = element_line(size = 1),
                           axis.ticks.length=unit(.25, "cm"),
                           panel.spacing = unit(1, "lines"))

## * export
for(iExtension in c("pdf","png")){ ## iExtension <- "pdf"
    ggsave(ggBenefit1, filename = paste0("Figures/ggBenefit-ChemoVSChemo-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggBenefit1.bis, filename = paste0("Figures/ggBenefit-ChemoVSChemo-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggPower1, filename = paste0("Figures/ggPower-ChemoVSChemo-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggPower1.bis, filename = paste0("Figures/ggPower-ChemoVSChemo-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggBenefit2, filename = paste0("Figures/ggBenefit-ChemoVSImmuno-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggBenefit2.bis, filename = paste0("Figures/ggBenefit-ChemoVSImmuno-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggPower2, filename = paste0("Figures/ggPower-ChemoVSImmuno-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggPower2.bis, filename = paste0("Figures/ggPower-ChemoVSImmuno-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggBenefit3, filename = paste0("Figures/ggBenefit-ImmunoVSImmuno-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggBenefit3.bis, filename = paste0("Figures/ggBenefit-ImmunoVSImmuno-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggPower3, filename = paste0("Figures/ggPower-ImmunoVSImmuno-scenario0.",iExtension), height = 6, width = 10)
    ggsave(ggPower3.bis, filename = paste0("Figures/ggPower-ImmunoVSImmuno-scenario123.",iExtension), height = 6, width = 10)

    ggsave(ggType1, filename = paste0("Figures/ggType1.",iExtension), height = 6, width = 8)
}
##----------------------------------------------------------------------
### FIGURE.R ends here
