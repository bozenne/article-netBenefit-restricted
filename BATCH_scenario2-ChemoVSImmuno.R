## * Header
## ** interactive
## path <- "h:/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/"
## setwd(path)
## source("BATCH_scenario2-ChemoVSImmuno.R")
## ** slurm
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## sbatch -a 1-1 -J 'scenario2-ChemoVSImmuno' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_scenario2-ChemoVSImmuno.R /dev/null 
## ** BATCH
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## R CMD BATCH --vanilla '--args iter_sim=1 n.iter_sim=10' BATCH_scenario2-ChemoVSImmuno.R output/scenario2-ChemoVSImmuno/R-ChemoVSImmuno-1.Rout &
## ** BATCH loop
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## for ITER in `seq 1 10`;
## do
## eval 'R CMD BATCH --vanilla "--args iter_sim='$ITER' n.iter_sim=10" BATCH_scenario2-ChemoVSImmuno.R output/scenario2-ChemoVSImmuno/R-ChemoVSImmuno-'$ITER'.Rout &'
## done

## [32] 3353501
## [33] 3353502
## [34] 3353503
## [35] 3353504
## [36] 3353505
## [37] 3353506
## [38] 3353507
## [39] 3353508
## [40] 3353509
## [41] 3353510

rm(list = ls())
gc()

## * Seed
args <- commandArgs(TRUE) ## BATCH MODE
if(length(args)>0){
    for (arg in args){
        eval(parse(text=arg))
    }
}else{ ## SLUMR
    iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
} ## interactive
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * Prepare export
path <- "."
path.res <- file.path(path,"Results","scenario2-ChemoVSImmuno")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","scenario2-ChemoVSImmuno")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * Libraries
require(survival)
require(BuyseTest) ## install.packages("BuyseTest")
require(survRM2) ## install.packages("survRM2")
suppressMessages(require(FHtest))

## * Settings
n.sim <- 500

Tps.inclusion <- 12 
Restriction.time_list <- c(12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60) ## every 3 months
Threshold_list <- c(0,6,12,18,24) ## 0,6,12
TpsFin <- 60

grid <- expand.grid(restrictionTime = Restriction.time_list,
                    threshold = Threshold_list,
                    scenario = 0)

grid <- rbind(grid,
              cbind(restrictionTime = Restriction.time_list,
                    threshold = 0.1*Restriction.time_list,
                    scenario = 1))

## pour faire varier les ratio threshold/restriction time
grid <- rbind(grid,
              cbind(restrictionTime = Restriction.time_list,
                    threshold = 0.05*Restriction.time_list,
                    scenario = 2))

grid <- rbind(grid,
              cbind(restrictionTime = Restriction.time_list,
                    threshold = 0.2*Restriction.time_list,
                    scenario = 3))

## pour obtenir la valeur exacte du net benefit sans censure
grid <- rbind(grid,
              cbind(restrictionTime = 1000,
                    threshold = Threshold_list,
                    scenario = 4))


## * Loop
res <- NULL
for(iSim in 1:n.sim){ ## iSim <- 1
    cat(iSim,": ")
    for (iGrid in 1:NROW(grid)){ ## iGrid <- 1
        cat("*")

        iThreshold <- grid$threshold[iGrid]
        iTime <- grid$restrictionTime[iGrid]
        iScenario <- grid$scenario[iGrid]

        ## ** Generate data
        HazC <- 0.1
	TpsFin <- iTime
        HazT2 <- HazC*0.95
        HazT3 <- HazC*0.90
        HazT4 <- HazC*0.65
        HazT5 <- HazC*0.1
        t1 <- 1
        t2 <- 3
        t3 <- 9
        t4 <- 24
        n.Treatment <- 200
        n.Control <- 200
        n <- n.Treatment+n.Control
        group <- c(rep(1, n.Treatment),rep(0, n.Control))
      
        TimeEvent.Ctr <- rexp(n.Control,HazC)
        TimeEvent.Tr1 <- rexp(n.Control,HazC)
        TimeEvent.Tr2 <- rexp(n.Control,HazT2)
        TimeEvent.Tr3 <- rexp(n.Control,HazT3)
        TimeEvent.Tr4 <- rexp(n.Control,HazT4)
        TimeEvent.Tr5 <- rexp(n.Control,HazT5)
      
        TimeEvent.Tr<-ifelse(TimeEvent.Tr1<t1,TimeEvent.Tr1,
                      ifelse(t1+TimeEvent.Tr2<t2,t1+TimeEvent.Tr2,
                      ifelse(t2+TimeEvent.Tr3<t3,t2+TimeEvent.Tr3,
                      ifelse(t3+TimeEvent.Tr4<t4,t3+TimeEvent.Tr4,
                             t4+TimeEvent.Tr5))))
      
        TimeEvent <- c(TimeEvent.Tr,TimeEvent.Ctr)
        Time.Cens <- runif(n,TpsFin-Tps.inclusion,TpsFin) ## varier temps de censure
        Time <- pmin(Time.Cens,TimeEvent)
        Event <- Time==TimeEvent
        Event <- as.numeric(Event)
        tab <- data.frame(group,Time,Event)
        ## Taux.cens.reel<-1-mean(Event)
  
        ## ** Analysis using LR
        LR <- (survdiff(Surv(time=Time, event=Event) ~ group, data=tab, rho=0))
        pval.LR <- 1 - pchisq(LR$chisq, 1) 
        Taux.cens.reel <- 1-mean(Event)
  
        ## ** Analysis using NBPeron
        NBPeron <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        NBPeron.confint <- confint(NBPeron)

        ## ** Analysis using RMST
        RMST <- rmst2(time=Time, status=Event, arm=group, tau = NULL, covariates = NULL, alpha = 0.05)
        pval.RMSTdif <- RMST[["unadjusted.result"]][1,4]
        pval.RMSTratio <- RMST[["unadjusted.result"]][2,4]
  
        ## ** Analysis using WLR
        WLR <- try(FHtestrcc(Surv(time=Time, event=Event) ~ group, data = tab, rho = 0,lambda = 1))
        if(inherits(WLR,"try-error")){
            pval.WLR <- NA
        }else{
            pval.WLR <- WLR$pvalue
        }

        ## ** Analysis using RNBPeron
        RNBPeron <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold, restriction=iTime),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        RNBPeron.confint <- confint(RNBPeron)
  
        ## ** Gather results
        res <- rbind(res,c(iteration = iSim,
                           Restriction_time = iTime,
                           scenario = iScenario,
                           Threshold = iThreshold,
                           "Taux censure reel" = Taux.cens.reel,
                           estimate.NBPeron = NBPeron.confint[,"estimate"],
                           se.NBPeron = NBPeron.confint[,"se"],
                           lower.NBPeron = NBPeron.confint[,"lower.ci"],
                           upper.NBPeron = NBPeron.confint[,"upper.ci"],
                           pval.NBPeron = NBPeron.confint[,"p.value"],
                           estimate.RNBPeron = RNBPeron.confint[,"estimate"],
                           se.RNBPeron = RNBPeron.confint[,"se"],
                           lower.RNBPeron = RNBPeron.confint[,"lower.ci"],
                           upper.RNBPeron = RNBPeron.confint[,"upper.ci"],
                           pval.RNBPeron = RNBPeron.confint[,"p.value"],
                           pval.LOGRANK = pval.LR,
                           pval.WeightedLOGRANK = pval.WLR,
                           pval.RMSTdif = pval.RMSTdif,
                           pval.RMSTratio = pval.RMSTratio))

    }
    saveRDS(res, file = file.path(path.res,paste0("simul_ChemoVSImmuno_",iter_sim,"(tempo).rds")))
    cat("\n")
}

## * Export
saveRDS(res, file = file.path(path.res,paste0("simul_ChemoVSImmuno_",iter_sim,".rds")))

## * R version
print(sessionInfo())


	
