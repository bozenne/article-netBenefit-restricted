## * Header
## ** interactive
## path <- "h:/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/"
## setwd(path)
## source("BATCH_scenario1.R")
## ** slurm
## cd ucph/hdir/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/
## sbatch -a 1-1 -J 'scenario1-ChemoVSChemo' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_scenario1-ChemoVSChemo.R /dev/null 
## ** BATCH
## cd ucph/hdir/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/
## R CMD BATCH --vanilla '--args iter_sim=1 n.iter_sim=10' BATCH_scenario1-ChemoVSChemo.R output/scenario1-ChemoVSChemo/R-ChemoVSChemo-1.Rout &
## ** BATCH loop
## cd ucph/hdir/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/
## for ITER in `seq 1 10`;
## do
## eval 'R CMD BATCH --vanilla "--args iter_sim='$ITER' n.iter_sim=10" BATCH_scenario1-ChemoVSChemo.R output/scenario1-ChemoVSChemo/R-ChemoVSChemo-'$ITER'.Rout &'
## done

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
path.res <- file.path(path,"Results","scenario1-ChemoVSChemo")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","scenario1-ChemoVSChemo")
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
require(FHtest)

## * Settings
n.sim <- 500

Tps.inclusion <- 12 
Restriction.time_list <- c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60) ## every 3 months
Threshold_list <- c(0,6,12,18) ## 0,6,12
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
for(iSim in 1:n.sim){
    cat(iSim,": ")
    for (iGrid in 1:NROW(grid)){ ## Restriction.time <- 12
        cat("*")
        iThreshold <- grid$threshold[iGrid]
        iTime <- grid$restrictionTime[iGrid]
        iScenario <- grid$scenario[iGrid]
        
        ## ** Generate data 
        TpsFin <- iTime
        HazC <- 0.085
        HazT <- 0.0595
        n.Treatment <- 200
        n.Control <- 200
        n <- n.Treatment+n.Control
        group <- c(rep(1, n.Treatment),rep(0, n.Control))
        TimeEvent.Ctr <- rexp(n.Control,HazC)
        TimeEvent.Tr <- rexp(n.Control,HazT)
  
        TimeEvent <- c(TimeEvent.Tr,TimeEvent.Ctr)
        Time.Cens <- runif(n,TpsFin-Tps.inclusion,TpsFin) #varier temps de censure
        Time <- pmin(Time.Cens,TimeEvent)
        Event <- Time==TimeEvent
        Event <- as.numeric(Event)
        tab <- data.frame(group,Time,Event)
        Taux.cens.reel <- 1-mean(Event)
  
        ## ** Analysis using LR
        LR <- (survdiff(Surv(time=Time, event=Event) ~ group, data=tab,rho=0))
        pval.LR <- 1 - pchisq(LR$chisq, 1) 
        Taux.cens.reel <- 1-mean(Event)
  
        ## ** Analysis using RNBGehan
        NBGehan <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold),
                              method.inference = "u-statistic", scoring.rule = "Gehan", trace = 0)
        NBGehan.confint <- confint(NBGehan)
        
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
                           estimate.NBGehan = NBGehan.confint[,"estimate"],
                           se.NBGehan = NBGehan.confint[,"se"],
                           lower.NBGehan = NBGehan.confint[,"lower.ci"],
                           upper.NBGehan = NBGehan.confint[,"upper.ci"],
                           pval.NBGehan = NBGehan.confint[,"p.value"],
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
    saveRDS(res, file = file.path(path.res,paste0("simul_ChemoVSChemo_",iter_sim,"(tempo).rds")))
    cat("\n")
}

## * Export
saveRDS(res, file = file.path(path.res,paste0("simul_ChemoVSChemo_",iter_sim,".rds")))

## * R version
print(sessionInfo())

