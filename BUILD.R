### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  7 2022 (11:10) 
## Version: 
## Last-Updated: nov 21 2022 (14:17) 
##           By: Brice Ozenne
##     Update #: 70
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)

if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    setwd("x:/GPC/article-restricted/")    
}else if(system("whoami",intern=TRUE)=="hpl802"){
    ## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
}else{
    setwd("C:/Users/max/Desktop/simulation_peron/simulation-article-restricted")
}

## * function used to collect the results from different files
loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = TRUE){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)

    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
        return(iOut)
        }
    }))
    out <- do.call(rbind, ls.out)
    return(out)
}

## * Scenario 1
dt.sc1 <- loadRes("Results/scenario1-ChemoVSChemo")
setnames(dt.sc1, new = "censure", old = "Taux censure reel")
setnames(dt.sc1, new = "threshold", old = "Threshold")
setnames(dt.sc1, new = "followUp", old = "FollowUp_time")

dt.sc1[,.(estimate.NBPeron,estimate.RNBPeron24)]
## any(duplicated(dt.sc1[,unique(seed),by=c("iteration","file")]$V1))

dtS.sc1 <- dt.sc1[, .(rep = .N, censure = mean(censure),
                      ## log-rank
                      power5.logrank = mean(pval.LOGRANK<=0.05),
                      power1.logrank = mean(pval.LOGRANK<=0.01),
                      ## weighted log-rank
                      power5.wlogrank = mean(pval.WeightedLOGRANK<=0.05, na.rm = TRUE),
                      power1.wlogrank = mean(pval.WeightedLOGRANK<=0.01, na.rm = TRUE),
                      ## restricted mean survival time
                      power5.rmstDiff = mean(pval.RMSTdif<=0.05),
                      power1.rmstDiff = mean(pval.RMSTdif<=0.01),
                      power5.rmstRatio = mean(pval.RMSTratio<=0.05),
                      power1.rmstRatio = mean(pval.RMSTratio<=0.01),
                      ## net benefit
                      power5.nbPeron = mean(pval.NBPeron<=0.05),
                      power1.nbPeron = mean(pval.NBPeron<=0.01),
                      GS.nbPeron = mean(GS.NBPeron),
                      estimate.nbPeron = mean(estimate.NBPeron),
                      coverage.nbPeron = mean((lower.NBPeron <= mean(estimate.NBPeron))*(mean(upper.NBPeron) <= estimate.NBPeron)),
                      ## net benefit (hierarchical)
                      power5.nbPeronTox1 = mean(pval.NBPeronTox1<=0.05), ## under the null
                      power1.nbPeronTox1 = mean(pval.NBPeronTox1<=0.01),
                      GS.nbPeronTox1 = mean(GS.NBPeronTox1),
                      estimate.nbPeronTox1 = mean(estimate.NBPeronTox1),
                      coverage.nbPeronTox1 = mean((lower.NBPeronTox1 <= mean(estimate.NBPeronTox1))*(mean(upper.NBPeronTox1) <= estimate.NBPeronTox1)),
                      power5.nbPeronTox2 = mean(pval.NBPeronTox2<=0.05), ## under the alternative
                      power1.nbPeronTox2 = mean(pval.NBPeronTox2<=0.01),
                      GS.nbPeronTox2 = mean(GS.NBPeronTox2),
                      estimate.nbPeronTox2 = mean(estimate.NBPeronTox2),
                      coverage.nbPeronTox2 = mean((lower.NBPeronTox2 <= mean(estimate.NBPeronTox2))*(mean(upper.NBPeronTox2) <= estimate.NBPeronTox2)),
                      ## restricted net benefit
                      power5.rnbPeron24 = mean(pval.RNBPeron24<=0.05), ## at 24 months
                      power1.rnbPeron24 = mean(pval.RNBPeron24<=0.01),
                      GS.rnbPeron24 = mean(GS.rNBPeron24),
                      estimate.rnbPeron24 = mean(estimate.RNBPeron24),
                      coverage.nbPeron24 = mean((lower.RNBPeron24 <= mean(estimate.RNBPeron24))*(mean(estimate.RNBPeron24) <= upper.RNBPeron24)),
                      power5.rnbPeron36 = mean(pval.RNBPeron36<=0.05), ## at 36 months
                      power1.rnbPeron36 = mean(pval.RNBPeron36<=0.01),
                      GS.rnbPeron36 = mean(GS.rNBPeron36),
                      estimate.rnbPeron36 = mean(estimate.RNBPeron36),
                      coverage.nbPeron36 = mean((lower.RNBPeron36 <= mean(estimate.RNBPeron36))*(mean(estimate.RNBPeron36) <= upper.RNBPeron36)),
                      power5.rnbPeron48 = mean(pval.RNBPeron48<=0.05), ## at 48 months
                      power1.rnbPeron48 = mean(pval.RNBPeron48<=0.01),
                      GS.rnbPeron48 = mean(GS.rNBPeron48),
                      estimate.rnbPeron48 = mean(estimate.RNBPeron48),
                      coverage.nbPeron48 = mean((lower.RNBPeron48 <= mean(estimate.RNBPeron48))*(mean(estimate.RNBPeron48) <= upper.RNBPeron48))
                      ), by = c("scenario","threshold","followUp") ]

dtS.sc1[dtS.sc1$threshold==0 & dtS.sc1$followUp==1000, ]
##    scenario threshold rtime   rep censure power.logrank power.wlogrank
## 1:        4         0  1000 10000       0        0.9479      0.8703222
##    power.rmstDiff power.rmstRatio power.nbGehan estimate.nbGehan
## 1:         0.9495          0.9512        0.8721        0.1771988
##    coverage.nbGehan power.rnbPeron estimate.rnbPeron coverage.nbPeron
## 1:            0.953         0.8721         0.1771988            0.953

dtS.sc1[scenario==4, .(threshold,followUp,rep,censure,estimate.nbPeron,estimate.rnbPeron24,estimate.rnbPeron36,estimate.rnbPeron48)]

## dtS.sc1[scenario == 0 & threshold == 0, .(followUp, estimate.nbPeron,estimate.rnbPeron24,estimate.rnbPeron36,estimate.rnbPeron48)]
## dtS.sc1[scenario == 4 & threshold == 0, .(followUp, estimate.nbPeron,estimate.rnbPeron24,estimate.rnbPeron36,estimate.rnbPeron48)]


## * Scenario 2
dt.sc2 <- loadRes("Results/scenario2-ChemoVSImmuno", tempo.file = TRUE)
setnames(dt.sc2, new = "censure", old = "Taux censure reel")
setnames(dt.sc2, new = "threshold", old = "Threshold")
setnames(dt.sc2, new = "followUp", old = "FollowUp_time")


dtS.sc2 <- dt.sc2[, .(rep = .N, censure = mean(censure),
                      ## log-rank                      
                      power5.logrank = mean(pval.LOGRANK<=0.05),
                      power1.logrank = mean(pval.LOGRANK<=0.01),
                      ## weighted log-rank
                      power5.wlogrank = mean(pval.WeightedLOGRANK<=0.05, na.rm = TRUE),
                      power1.wlogrank = mean(pval.WeightedLOGRANK<=0.01, na.rm = TRUE),
                      ## restricted mean survival time
                      power5.rmstDiff = mean(pval.RMSTdif<=0.05),
                      power1.rmstDiff = mean(pval.RMSTdif<=0.01),
                      power5.rmstRatio = mean(pval.RMSTratio<=0.05),
                      power1.rmstRatio = mean(pval.RMSTratio<=0.01),
                      ## net benefit
                      power5.nbPeron = mean(pval.NBPeron<=0.05),
                      power1.nbPeron = mean(pval.NBPeron<=0.01),
                      GS.nbPeron = mean(GS.NBPeron),
                      estimate.nbPeron = mean(estimate.NBPeron),
                      coverage.nbPeron = mean((lower.NBPeron <= mean(estimate.NBPeron))*(mean(upper.NBPeron) <= estimate.NBPeron)),
                      ## net benefit (hierarchical)
                      power5.nbPeronTox1 = mean(pval.NBPeronTox1<=0.05), ## under the null
                      power1.nbPeronTox1 = mean(pval.NBPeronTox1<=0.01),
                      GS.nbPeronTox1 = mean(GS.NBPeronTox1),
                      GS.nbPeronTox1 = mean(GS.NBPeronTox1),
                      estimate.nbPeronTox1 = mean(estimate.NBPeronTox1),
                      coverage.nbPeronTox1 = mean((lower.NBPeronTox1 <= mean(estimate.NBPeronTox1))*(mean(upper.NBPeronTox1) <= estimate.NBPeronTox1)),
                      power5.nbPeronTox2 = mean(pval.NBPeronTox2<=0.05), ## under the alternative
                      power1.nbPeronTox2 = mean(pval.NBPeronTox2<=0.01),
                      GS.nbPeronTox2 = mean(GS.NBPeronTox2),
                      estimate.nbPeronTox2 = mean(estimate.NBPeronTox2),
                      coverage.nbPeronTox2 = mean((lower.NBPeronTox2 <= mean(estimate.NBPeronTox2))*(mean(upper.NBPeronTox2) <= estimate.NBPeronTox2)),
                      ## restricted net benefit
                      power5.rnbPeron24 = mean(pval.RNBPeron24<=0.05), ## at 24 months
                      power1.rnbPeron24 = mean(pval.RNBPeron24<=0.01),
                      GS.rnbPeron24 = mean(GS.rNBPeron24),
                      estimate.rnbPeron24 = mean(estimate.RNBPeron24),
                      coverage.nbPeron24 = mean((lower.RNBPeron24 <= mean(estimate.RNBPeron24))*(mean(estimate.RNBPeron24) <= upper.RNBPeron24)),
                      power5.rnbPeron36 = mean(pval.RNBPeron36<=0.05), ## at 36 months
                      power1.rnbPeron36 = mean(pval.RNBPeron36<=0.01),
                      GS.rnbPeron36 = mean(GS.rNBPeron36),
                      estimate.rnbPeron36 = mean(estimate.RNBPeron36),
                      coverage.nbPeron36 = mean((lower.RNBPeron36 <= mean(estimate.RNBPeron36))*(mean(estimate.RNBPeron36) <= upper.RNBPeron36)),
                      power5.rnbPeron48 = mean(pval.RNBPeron48<=0.05), ## at 48 months
                      power1.rnbPeron48 = mean(pval.RNBPeron48<=0.01),
                      GS.rnbPeron48 = mean(GS.rNBPeron48),
                      estimate.rnbPeron48 = mean(estimate.RNBPeron48),
                      coverage.nbPeron48 = mean((lower.RNBPeron48 <= mean(estimate.RNBPeron48))*(mean(estimate.RNBPeron48) <= upper.RNBPeron48))
                      ), by = c("scenario","threshold","followUp") ]

dtS.sc2[dtS.sc2$threshold==0 & dtS.sc2$followUp==1000, ]
##    scenario threshold rtime   rep  censure power.logrank power.wlogrank
## 1:        4         0  1000 10000 6.75e-06        0.9502      0.9996993
##    power.rmstDiff power.rmstRatio power.nbGehan estimate.nbGehan power.rnbPeron
## 1:         0.9998          0.9999        0.3079       0.08495964         0.3079
##    estimate.rnbPeron
## 1:        0.08495964

## * Scenario 3
dt.sc3 <- loadRes("Results/scenario3-ImmunoVSImmuno")
setnames(dt.sc3, new = "censure", old = "Taux censure reel")
setnames(dt.sc3, new = "threshold", old = "Threshold")
setnames(dt.sc3, new = "followUp", old = "FollowUp_time")

dtS.sc3 <- dt.sc3[, .(rep = .N, censure = mean(censure),
                      ## log-rank                      
                      power5.logrank = mean(pval.LOGRANK<=0.05),
                      power1.logrank = mean(pval.LOGRANK<=0.01),
                      ## weighted log-rank
                      power5.wlogrank = mean(pval.WeightedLOGRANK<=0.05, na.rm = TRUE),
                      power1.wlogrank = mean(pval.WeightedLOGRANK<=0.01, na.rm = TRUE),
                      ## restricted mean survival time
                      power5.rmstDiff = mean(pval.RMSTdif<=0.05),
                      power1.rmstDiff = mean(pval.RMSTdif<=0.01),
                      power5.rmstRatio = mean(pval.RMSTratio<=0.05),
                      power1.rmstRatio = mean(pval.RMSTratio<=0.01),
                      ## net benefit
                      power5.nbPeron = mean(pval.NBPeron<=0.05),
                      power1.nbPeron = mean(pval.NBPeron<=0.01),
                      GS.nbPeron = mean(GS.NBPeron),
                      estimate.nbPeron = mean(estimate.NBPeron),
                      coverage.nbPeron = mean((lower.NBPeron <= mean(estimate.NBPeron))*(mean(upper.NBPeron) <= estimate.NBPeron)),
                      ## net benefit (hierarchical)
                      power5.nbPeronTox1 = mean(pval.NBPeronTox1<=0.05), ## under the null
                      power1.nbPeronTox1 = mean(pval.NBPeronTox1<=0.01),
                      GS.nbPeronTox1 = mean(GS.NBPeronTox1),
                      estimate.nbPeronTox1 = mean(estimate.NBPeronTox1),
                      coverage.nbPeronTox1 = mean((lower.NBPeronTox1 <= mean(estimate.NBPeronTox1))*(mean(upper.NBPeronTox1) <= estimate.NBPeronTox1)),
                      power5.nbPeronTox2 = mean(pval.NBPeronTox2<=0.05), ## under the alternative
                      power1.nbPeronTox2 = mean(pval.NBPeronTox2<=0.01),
                      GS.nbPeronTox2 = mean(GS.NBPeronTox2),
                      estimate.nbPeronTox2 = mean(estimate.NBPeronTox2),
                      coverage.nbPeronTox2 = mean((lower.NBPeronTox2 <= mean(estimate.NBPeronTox2))*(mean(upper.NBPeronTox2) <= estimate.NBPeronTox2)),
                      ## restricted net benefit
                      power5.rnbPeron24 = mean(pval.RNBPeron24<=0.05), ## at 24 months
                      power1.rnbPeron24 = mean(pval.RNBPeron24<=0.01),
                      GS.rnbPeron24 = mean(GS.rNBPeron24),
                      estimate.rnbPeron24 = mean(estimate.RNBPeron24),
                      coverage.nbPeron24 = mean((lower.RNBPeron24 <= mean(estimate.RNBPeron24))*(mean(estimate.RNBPeron24) <= upper.RNBPeron24)),
                      power5.rnbPeron36 = mean(pval.RNBPeron36<=0.05), ## at 36 months
                      power1.rnbPeron36 = mean(pval.RNBPeron36<=0.01),
                      GS.rnbPeron36 = mean(GS.rNBPeron36),
                      estimate.rnbPeron36 = mean(estimate.RNBPeron36),
                      coverage.nbPeron36 = mean((lower.RNBPeron36 <= mean(estimate.RNBPeron36))*(mean(estimate.RNBPeron36) <= upper.RNBPeron36)),
                      power5.rnbPeron48 = mean(pval.RNBPeron48<=0.05), ## at 48 months
                      power1.rnbPeron48 = mean(pval.RNBPeron48<=0.01),
                      GS.rnbPeron48 = mean(GS.rNBPeron48),
                      estimate.rnbPeron48 = mean(estimate.RNBPeron48),
                      coverage.nbPeron48 = mean((lower.RNBPeron48 <= mean(estimate.RNBPeron48))*(mean(estimate.RNBPeron48) <= upper.RNBPeron48))
                      ),
                  by = c("scenario","threshold","followUp") ]

dtS.sc3[dtS.sc3$threshold==0 & dtS.sc3$followUp==1000, ]
##    scenario threshold rtime   rep    censure power.logrank power.wlogrank
## 1:        4         0  1000 10000 0.00020325          0.94      0.8674699
##    power.rmstDiff power.rmstRatio power.nbGehan estimate.nbGehan power.rnbPeron
## 1:         0.9173          0.9202        0.8675        0.1765435         0.8675
##    estimate.rnbPeron sd.rnbPeron
## 1:         0.1765435   0.0569079

## * Scenario 4
dt.sc4 <- loadRes("Results/type1_error")
setnames(dt.sc4, new = "censure", old = "Taux censure reel")
setnames(dt.sc4, new = "threshold", old = "Threshold")
setnames(dt.sc4, new = "rtime", old = "Restriction_time")

dtS.sc4 <- dt.sc4[, .(rep = .N, censure = mean(censure),
                      power5.logrank = mean(pval.LOGRANK<=0.05),
                      power1.logrank = mean(pval.LOGRANK<=0.01),
                      power5.wlogrank = mean(pval.WeightedLOGRANK<=0.05, na.rm = TRUE),
                      power1.wlogrank = mean(pval.WeightedLOGRANK<=0.01, na.rm = TRUE),
                      power5.rmstDiff = mean(pval.RMSTdif<=0.05),
                      power1.rmstDiff = mean(pval.RMSTdif<=0.01),
                      power5.rmstRatio = mean(pval.RMSTratio<=0.05),
                      power1.rmstRatio = mean(pval.RMSTratio<=0.01),
                      power5.rnbPeron = mean(pval.RNBPeron<=0.05),
                      power1.rnbPeron = mean(pval.RNBPeron<=0.01),
                      estimate.rnbPeron = mean(estimate.RNBPeron),
                      coverage.nbPeron = mean((lower.RNBPeron <= 0)*(0 <= upper.RNBPeron))
                      ), by = c("scenario","threshold","rtime") ]

## * export
saveRDS(dtS.sc1, file = "Results/simSummary-ChemoVSChemo.rds")
saveRDS(dtS.sc2, file = "Results/simSummary-ChemoVSImmuno.rds")
saveRDS(dtS.sc3, file = "Results/simSummary-ImmunoVSImmuno.rds")
saveRDS(dtS.sc4, file = "Results/simSummary-type1.rds")

saveRDS(dt.sc1, file = "Results/sim-ChemoVSChemo.rds")
saveRDS(dt.sc2, file = "Results/sim-ChemoVSImmuno.rds")
saveRDS(dt.sc3, file = "Results/sim-ImmunoVSImmuno.rds")
saveRDS(dt.sc4, file = "Results/sim-type1.rds")

##----------------------------------------------------------------------
### BUILD.R ends here
