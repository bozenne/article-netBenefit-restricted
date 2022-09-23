## * Header
## ** interactive
## path <- "h:/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/"
## setwd(path)
## source("BATCH_scenario1-ChemoVSChemo.R")
## ** slurm
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## sbatch -a 1-1 -J 'scenario1-ChemoVSChemo' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_scenario1-ChemoVSChemo.R /dev/null 
## ** BATCH
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## R CMD BATCH --vanilla '--args iter_sim=1 n.iter_sim=10' BATCH_scenario1-ChemoVSChemo.R output/scenario1-ChemoVSChemo/R-ChemoVSChemo-1.Rout &
## ** BATCH loop
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
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
suppressMessages(require(FHtest))

## * Settings
n.sim <- 100

Tps.inclusion <- 12 
FollowUp.time_list <- c(12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60) ## every 3 months
Threshold_list <- c(0,6,12,18,24) ## every 6 months

grid <- expand.grid(FollowUpTime = FollowUp.time_list,
                    threshold = Threshold_list,
                    scenario = 0)

## pour faire varier les ratio threshold/restriction time
grid <- rbind(grid,
              cbind(FollowUpTime = FollowUp.time_list,
                  threshold = 0.2*FollowUp.time_list,
                    scenario = 3))
grid <- rbind(grid,
              cbind(FollowUpTime = FollowUp.time_list,
                    threshold = 0.3*FollowUp.time_list,
                    scenario = 5))

grid <- rbind(grid,
              cbind(FollowUpTime = FollowUp.time_list,
                    threshold = 0.4*FollowUp.time_list,
                    scenario = 6))

## pour obtenir la valeur exacte du net benefit sans censure
grid <- rbind(grid,
              cbind(FollowUpTime = 1000,
                    threshold = Threshold_list,
                    scenario = 4))


## * Loop
res <- NULL
for(iSim in 1:n.sim){ ## iSim <- 1 
    cat(iSim,": ")
    for (iGrid in 1:NROW(grid)){ ## iGrid <- 1
        cat("*")
        iThreshold <- grid$threshold[iGrid]
        iFollowUpTime <- grid$FollowUpTime[iGrid]
        iScenario <- grid$scenario[iGrid]
        
        ## ** Generate data 
        HazC <- 0.085
        HazT <- 0.0595
        
        ## Tox
        ptoxC <- c(0.3,0.3)
        ptoxT <- c(0.3,0.1)
        
        n.Treatment <- 200
        n.Control <- 200
        n <- n.Treatment+n.Control
        group <- c(rep(1, n.Treatment),rep(0, n.Control))
        TimeEvent.Ctr <- rexp(n.Control,HazC)
        TimeEvent.Tr <- rexp(n.Control,HazT)
        
        Toxevent1.Ctr <- rbinom(n.Control,1,ptoxC[1])
        Toxevent1.Tr <- rbinom(n.Treatment,1,ptoxT[1])
        Toxevent2.Ctr <- rbinom(n.Control,1,ptoxC[2])
        Toxevent2.Tr <- rbinom(n.Treatment,1,ptoxT[2])
  
        TimeEvent <- c(TimeEvent.Tr,TimeEvent.Ctr)
        Toxevent1 <- c(Toxevent1.Tr,Toxevent1.Ctr)
        Toxevent2 <- c(Toxevent2.Tr,Toxevent2.Ctr)
        
        Time.Cens <- runif(n,iFollowUpTime-Tps.inclusion,iFollowUpTime) #varier temps de censure
        Time <- pmin(Time.Cens,TimeEvent)
        Event <- Time==TimeEvent
        Event <- as.numeric(Event)
        tab <- data.frame(group = group,
                          Time = Time,
                          Event = Event,
                          Toxevent1 = Toxevent1,
                          Toxevent2 = Toxevent2)
        Taux.cens.reel <- 1-mean(Event)
  
        ## ** Analysis using LR
        LR <- (survdiff(Surv(time=Time, event=Event) ~ group, data=tab,rho=0))
        pval.LR <- 1 - pchisq(LR$chisq, 1) 
        Taux.cens.reel <- 1-mean(Event)
  
        ## ** Analysis using NBPeron
        NBPeron <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        NBPeron.confint <- confint(NBPeron)
        
        ## ** Analysis using NBPeron + toxicity
        NBPeronTox1 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold)+ B(Toxevent1, operator = "<0"),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        NBPeronTox1.confint <- confint(NBPeronTox1)
        
        NBPeronTox2 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold)+ B(Toxevent2, operator = "<0"),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        NBPeronTox2.confint <- confint(NBPeronTox2)
        
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
        if(iFollowUpTime<=24){
            RNBPeron24.confint <- NBPeron.confint
        }else{
            RNBPeron24 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold, restriction=24),
                                    method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
            RNBPeron24.confint <- confint(RNBPeron24)
        }

        if(iFollowUpTime<=36){
            RNBPeron36.confint <- NBPeron.confint
        }else{
            RNBPeron36 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold, restriction=36),
                                    method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
            RNBPeron36.confint <- confint(RNBPeron36)
        }
  
        if(iFollowUpTime<=48){
            RNBPeron48.confint <- NBPeron.confint
        }else{
            RNBPeron48 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold, restriction=48),
                                    method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
            RNBPeron48.confint <- confint(RNBPeron48)
        }
        
        ## ** Gather results
        res <- rbind(res,c(iteration = iSim,
                           FollowUp_time = iFollowUpTime,
                           scenario = iScenario,
                           Threshold = iThreshold,
                           "Taux censure reel" = Taux.cens.reel,
                           ## NBPeron
                           estimate.NBPeron = NBPeron.confint[,"estimate"],
                           se.NBPeron = NBPeron.confint[,"se"],
                           lower.NBPeron = NBPeron.confint[,"lower.ci"],
                           upper.NBPeron = NBPeron.confint[,"upper.ci"],
                           pval.NBPeron = NBPeron.confint[,"p.value"],
                           ## NBPeronTox1
                           estimate.NBPeronTox1 = NBPeronTox1.confint[2,"estimate"],
                           se.NBPeronTox1 = NBPeronTox1.confint[2,"se"],
                           lower.NBPeronTox1 = NBPeronTox1.confint[2,"lower.ci"],
                           upper.NBPeronTox1 = NBPeronTox1.confint[2,"upper.ci"],
                           pval.NBPeronTox1 = NBPeronTox1.confint[2,"p.value"],
                           ## NBPeronTox2
                           estimate.NBPeronTox2 = NBPeronTox2.confint[2,"estimate"],
                           se.NBPeronTox2 = NBPeronTox2.confint[2,"se"],
                           lower.NBPeronTox2 = NBPeronTox2.confint[2,"lower.ci"],
                           upper.NBPeronTox2 = NBPeronTox2.confint[2,"upper.ci"],
                           pval.NBPeronTox2 = NBPeronTox2.confint[2,"p.value"],
                           ## RNBPeron24
                           estimate.RNBPeron24 = RNBPeron24.confint[,"estimate"],
                           se.RNBPeron24 = RNBPeron24.confint[,"se"],
                           lower.RNBPeron24 = RNBPeron24.confint[,"lower.ci"],
                           upper.RNBPeron24 = RNBPeron24.confint[,"upper.ci"],
                           pval.RNBPeron24 = RNBPeron24.confint[,"p.value"],
                           ## RNBPeron36
                           estimate.RNBPeron36 = RNBPeron36.confint[,"estimate"],
                           se.RNBPeron36 = RNBPeron36.confint[,"se"],
                           lower.RNBPeron36 = RNBPeron36.confint[,"lower.ci"],
                           upper.RNBPeron36 = RNBPeron36.confint[,"upper.ci"],
                           pval.RNBPeron36 = RNBPeron36.confint[,"p.value"],
                           ## RNBPeron48
                           estimate.RNBPeron48 = RNBPeron48.confint[,"estimate"],
                           se.RNBPeron48 = RNBPeron48.confint[,"se"],
                           lower.RNBPeron48 = RNBPeron48.confint[,"lower.ci"],
                           upper.RNBPeron48 = RNBPeron48.confint[,"upper.ci"],
                           pval.RNBPeron48 = RNBPeron48.confint[,"p.value"],
                           ## other
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

