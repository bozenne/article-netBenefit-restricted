## * Header
## ** interactive
## path <- "h:/SundKonsolidering_BioStatHome/Cluster/GPC/article-restricted/"
## setwd(path)
## source("BATCH_scenario2bis-ChemoVSImmuno.R")
## ** slurm
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## sbatch -a 1-1 -J 'scenario2bis-ChemoVSImmuno' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_scenario2bis-ChemoVSImmuno.R /dev/null 
## ** BATCH
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## R CMD BATCH --vanilla '--args iter_sim=1 n.iter_sim=10' BATCH_scenario2bis-ChemoVSImmuno.R output/scenario2bis-ChemoVSImmuno/R-ChemoVSImmuno-1.Rout &
## ** BATCH loop
## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## for ITER in `seq 1 10`;
## do
## eval 'R CMD BATCH --vanilla "--args iter_sim='$ITER' n.iter_sim=10" BATCH_scenario2bis-ChemoVSImmuno.R output/scenario2bis-ChemoVSImmuno/R-ChemoVSImmuno-'$ITER'.Rout &'
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

## * Arguments 
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

## * Prepare export
path <- "."
path.res <- file.path(path,"Results","scenario2bis-ChemoVSImmuno")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","scenario2bis-ChemoVSImmuno")
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
Threshold_list <- c(0,6,12,18,24) ## 0,6,12

grid <- expand.grid(restrictionTime = FollowUp.time_list,
                    threshold = Threshold_list,
                    scenario = 0)

## pour faire varier les ratio threshold/restriction time
grid <- rbind(grid,
              cbind(restrictionTime = FollowUp.time_list,
                    threshold = 0.2*FollowUp.time_list,
                    scenario = 3))

grid <- rbind(grid,
              cbind(restrictionTime = FollowUp.time_list,
                    threshold = 0.3*FollowUp.time_list,
                    scenario = 5))

grid <- rbind(grid,
              cbind(restrictionTime = FollowUp.time_list,
                    threshold = 0.4*FollowUp.time_list,
                    scenario = 6))

## pour obtenir la valeur exacte du net benefit sans censure
grid <- rbind(grid,
              cbind(restrictionTime = 1000,
                    threshold = Threshold_list,
                    scenario = 4))

## * Seed
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:1e5,size=n.iter_sim*n.sim,replace=FALSE)
## any(duplicated(seqSeed))
iSeed <- seqSeed[(n.sim*(iter_sim-1)+1):(n.sim*iter_sim)]

## * Loop
res <- NULL
for(iSim in 1:n.sim){ ## iSim <- 1
    cat(iSim," (seed=",iSeed[iSim],"): ",sep="")
    for (iGrid in 1:NROW(grid)){ ## iGrid <- 1
        cat("*")
        set.seed(iSeed[iSim])

        iThreshold <- grid$threshold[iGrid]
        iFollowUpTime <- grid$restrictionTime[iGrid]
        iScenario <- grid$scenario[iGrid]

        ## ** Generate data
        HazC <- 0.1
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
      
        TimeEvent.Ctr <- rexp(n.Control,HazC)
        TimeEvent.Tr1 <- rexp(n.Control,HazC)
        TimeEvent.Tr2 <- rexp(n.Control,HazT2)
        TimeEvent.Tr3 <- rexp(n.Control,HazT3)
        TimeEvent.Tr4 <- rexp(n.Control,HazT4)
        TimeEvent.Tr5 <- rexp(n.Control,HazT5)
      
        TimeEvent.Tr <- ifelse(TimeEvent.Tr1<t1,TimeEvent.Tr1,
                        ifelse(t1+TimeEvent.Tr2<t2,t1+TimeEvent.Tr2,
                        ifelse(t2+TimeEvent.Tr3<t3,t2+TimeEvent.Tr3,
                        ifelse(t3+TimeEvent.Tr4<t4,t3+TimeEvent.Tr4,
                               t4+TimeEvent.Tr5))))
	    
        ## Tox
        ptoxC <- c(0.3,0.3)
        ptoxT <- c(0.3,0.2)

        Toxevent1.Ctr <- rbinom(n.Control,1,ptoxC[1])
        Toxevent1.Tr <- rbinom(n.Treatment,1,ptoxT[1])
        Toxevent2.Ctr <- rbinom(n.Control,1,ptoxC[2])
        Toxevent2.Tr <- rbinom(n.Treatment,1,ptoxT[2])
      
        TimeEvent <- c(TimeEvent.Tr,TimeEvent.Ctr)
        Toxevent1 <- c(Toxevent1.Tr,Toxevent1.Ctr)
        Toxevent2 <- c(Toxevent2.Tr,Toxevent2.Ctr)
	 
        inclusion.date.Tr <- seq(0,n.Treatment-1,by=1)*2/3*iFollowUpTime/n.Treatment ## calendar date at which each treated patient starts
        inclusion.date.Ctr <- seq(0,n.Control-1,by=1)*2/3*iFollowUpTime/n.Control ## calendar date at which each control patient starts
        Time.cens.Tr <- iFollowUpTime-pmax(0,runif(n.Treatment,-2,2)+inclusion.date.Tr)
        Time.cens.Ctr <- iFollowUpTime-pmax(0,runif(n.Control,-2,2)+inclusion.date.Ctr)
        Time.cens <- c(Time.cens.Tr,Time.cens.Ctr)

        tab <- data.frame(group = c(rep(1, n.Treatment),rep(0, n.Control)),
                          Time0 = TimeEvent,
                          Event0 = 1,
                          Time = pmin(Time.Cens,TimeEvent),
                          Event = NA,
                          Toxevent1 = Toxevent1,
                          Toxevent2 = Toxevent2)
        tab$Event <- as.numeric(tab$Time==TimeEvent)
        Taux.cens.reel <- 1 - mean(tab$Event)
  
        ## ** Analysis using LR
        LR <- (survdiff(Surv(time=Time, event=Event) ~ group, data=tab, rho=0))
        pval.LR <- 1 - pchisq(LR$chisq, 1) 
  
        ## ** Analysis using NBPeron
        NBPeron <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold),
                              method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
        NBPeron.confint <- confint(NBPeron)

        ## without censoring
        NB.noCensoring <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold),
                                    method.inference = "none", trace = 0)
	    
       ## ** Analysis using NBPeron + toxicity
    	NBPeronTox1 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold) + B(Toxevent1, operator = "<0"),
                         method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
   	NBPeronTox1.confint <- confint(NBPeronTox1)

    	NBPeronTox2 <- BuyseTest(data=tab,group ~ TTE(Time, status=Event, iThreshold) + B(Toxevent2, operator = "<0"),
                         method.inference = "u-statistic", scoring.rule = "Peron", trace = 0)
   	NBPeronTox2.confint <- confint(NBPeronTox2)

        ## without censoring
        NBTox1.noCensoring <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold) + B(Toxevent1, operator = "<0"),
                                        method.inference = "none", trace = 0)
        NBTox2.noCensoring <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold) + B(Toxevent2, operator = "<0"),
                                        method.inference = "none", trace = 0)

        ## ** Analysis using RMST
        RMST <- rmst2(time=tab$Time, status=tab$Event, arm= tab$group, tau = NULL, covariates = NULL, alpha = 0.05)
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

        ## without censoring
        rNB.noCensoring24 <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold, restriction = 24),
                                       method.inference = "none", trace = 0)
        rNB.noCensoring36 <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold, restriction = 36),
                                       method.inference = "none", trace = 0)
        rNB.noCensoring48 <- BuyseTest(data=tab,group ~ TTE(Time0, status=Event0, iThreshold, restriction = 48),
                                       method.inference = "none", trace = 0)

        ## ** Gather results
        res <- rbind(res,c(iteration = iSim,
                           seed = iSeed[iSim],
                           FollowUp_time = iFollowUpTime,
                           scenario = iScenario,
                           Threshold = iThreshold,
                           "Taux censure reel" = Taux.cens.reel,
                           ## NBPeron
                           GS.NBPeron = as.numeric(coef(NB.noCensoring)),
                           estimate.NBPeron = NBPeron.confint[,"estimate"],
                           se.NBPeron = NBPeron.confint[,"se"],
                           lower.NBPeron = NBPeron.confint[,"lower.ci"],
                           upper.NBPeron = NBPeron.confint[,"upper.ci"],
                           pval.NBPeron = NBPeron.confint[,"p.value"],
                           ## NBPeronTox1
                           GS.NBPeronTox1 = as.numeric(utils::tail(coef(NBTox1.noCensoring),1)),
			   estimate.NBPeronTox1 = NBPeronTox1.confint[2,"estimate"],
                           se.NBPeronTox1 = NBPeronTox1.confint[2,"se"],
                           lower.NBPeronTox1 = NBPeronTox1.confint[2,"lower.ci"],
                           upper.NBPeronTox1 = NBPeronTox1.confint[2,"upper.ci"],
                           pval.NBPeronTox1 = NBPeronTox1.confint[2,"p.value"],
                           ## NBPeronTox2
                           GS.NBPeronTox2 = as.numeric(utils::tail(coef(NBTox2.noCensoring),1)),
                           estimate.NBPeronTox2 = NBPeronTox2.confint[2,"estimate"],
                           se.NBPeronTox2 = NBPeronTox2.confint[2,"se"],
                           lower.NBPeronTox2 = NBPeronTox2.confint[2,"lower.ci"],
                           upper.NBPeronTox2 = NBPeronTox2.confint[2,"upper.ci"],
                           pval.NBPeronTox2 = NBPeronTox2.confint[2,"p.value"],
                           ## RNBPeron24
                           GS.rNBPeron24 = as.numeric(coef(rNB.noCensoring24)),
                           estimate.RNBPeron24 = RNBPeron24.confint[,"estimate"],
                           se.RNBPeron24 = RNBPeron24.confint[,"se"],
                           lower.RNBPeron24 = RNBPeron24.confint[,"lower.ci"],
                           upper.RNBPeron24 = RNBPeron24.confint[,"upper.ci"],
                           pval.RNBPeron24 = RNBPeron24.confint[,"p.value"],
                           ## RNBPero36
                           GS.rNBPeron36 = as.numeric(coef(rNB.noCensoring36)),
                           estimate.RNBPeron36 = RNBPeron36.confint[,"estimate"],
                           se.RNBPeron36 = RNBPeron36.confint[,"se"],
                           lower.RNBPeron36 = RNBPeron36.confint[,"lower.ci"],
                           upper.RNBPeron36 = RNBPeron36.confint[,"upper.ci"],
                           pval.RNBPeron36 = RNBPeron36.confint[,"p.value"],
                           ## RNBPero48
                           GS.rNBPeron48 = as.numeric(coef(rNB.noCensoring48)),
                           estimate.RNBPeron48 = RNBPeron48.confint[,"estimate"],
                           se.RNBPeron48 = RNBPeron48.confint[,"se"],
                           lower.RNBPeron48 = RNBPeron48.confint[,"lower.ci"],
                           upper.RNBPeron48 = RNBPeron48.confint[,"upper.ci"],
                           pval.RNBPeron48 = RNBPeron48.confint[,"p.value"],
                           ## others
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


	
