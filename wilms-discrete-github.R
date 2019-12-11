
nwts <- read.table('./nwts-expanded.txt',header=TRUE)


stage34 <- ifelse(nwts$stage%in%c(1,2),0,1)
nwts <- data.frame(nwts,stage34=stage34)

head(nwts)

X <- data.frame(nwts$instit,nwts$histol,nwts$stage34,nwts$age,nwts$tumdiam,nwts$histol*nwts$stage34)
colnames(X) <- c('instit','histol','stage34','age','tumdiam','histol:stage34')


dim(nwts)
with(nwts, table(instit,histol))
1-sum(nwts$instit==nwts$histol)/nrow(nwts)


indMC <- 1

M <- 1000

# optVar <- 1 # histol
# optVar <- 2 # stage34
# optVar <- 3 # age
# optVar <- 4 # tumdiam
optVar <- 5 # histol:stage34

link <- 'cloglog'
n <- 400
  

### source functions
source('./source-functions-github.R')



obsY0 <- nwts$trel
obsD0 <- nwts$relaps
obsX0 <- cbind(nwts$histol,nwts$stage34,nwts$age,nwts$tumdiam,nwts$histol*nwts$stage34)
obsZ0 <- matrix(nwts$instit)

mean(obsY0[which(obsD0==1)]<=3)
mean(obsY0[which(obsD0==0)]<=3)

mean(obsY0<=3)
mean(obsY0<=3 & obsD0==1)
mean(obsY0<=3 & obsD0==0)


indOut <- which(obsY0<=3 & obsD0==0)
indOut <- nrow(nwts)+1

obsY <- obsY0[-indOut]
obsD <- obsD0[-indOut]
obsX <- obsX0[-indOut,]
obsZ <- as.matrix(obsZ0[-indOut,])


dfWilms <- nwts[-indOut,]


indTauCut <- 6
tauCut <- 1:indTauCut

obsY <- ceiling(2*obsY)

head(obsY)
table(obsY)
head(obsD)
head(obsX)
head(obsZ)

table(obsX[,1],obsZ)
1-sum(obsX[,1]==obsZ[,1])/length(obsY)

N <- length(obsY)

# n <- 400
r <- 1-n/N

n0 <- n/2
r0 <- 1-n/N/2

trueAlpha <- rep(0,indTauCut)
trueBeta <- rep(0,ncol(obsX))


### m-th MC loop

alphaMsRand <- alphaCpRand <- alphaIpwRand <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsRand <- betaCpRand <- betaIpwRand <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

alphaMsBal <- alphaCpBal <- alphaIpwBal <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsBal <- betaCpBal <- betaIpwBal <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

alphaMsPilot <- alphaCpPilot <- alphaIpwPilot <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsPilot <- betaCpPilot <- betaIpwPilot <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

alphaMsPre <- alphaCpPre <- alphaIpwPre <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsPre <- betaCpPre <- betaIpwPre <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

alphaMsOpt <- alphaCpOpt <- alphaIpwOpt <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsOpt <- betaCpOpt <- betaIpwOpt <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

alphaMsOrac <- alphaCpOrac <- alphaIpwOrac <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaMsOrac <- betaCpOrac <- betaIpwOrac <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))




stratAll <- stratDesign(obsY,obsD,obsZ,tauCut)
strat <- c()
for (l in 1:length(stratAll)) {
  strat[stratAll[[l]]] <- l
}
strat <- as.factor(strat)


tableFullMC <- tableBalMC <- tableOracMC <- tableOptMC <- tableAdaptMC <- matrix(0,nrow=(6+1),ncol=(4+1))



timeMC1 <- Sys.time()

for (m in 1:M) {
  
  time1 <- Sys.time()
  
  set.seed(M*(indMC-1)+m)
  
  print(paste('MC replication ',m,sep=''))
  
  
  
  time01 <- Sys.time()
  ### completely balanced sampling
  set.seed(M*(indMC-1)+m)
  balSam0 <- balSample(obsY,obsD,obsZ,r,tauCut,pamela=TRUE)
  
  balWt0 <- balSam0$wt
  balTable0 <- balSam0$tableStrat
  indStrat0 <- balSam0$indStrat
  
  indPhaseBal <- balSam0[c(1,2)]
  
  indPhase1Bal <- indPhaseBal[[1]]
  indPhase2Bal <- indPhaseBal[[2]]
  
  print('Bal')
  print(length(indPhase2Bal)/(length(indPhase1Bal)+length(indPhase2Bal)))
  
  
  
  
  # balanced sample
  fitBal <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseBal,tauCut,tol=0.001,wt=balWt0,link=link)
  
  alphaCpBal[m,] <- fitBal$alphaCp
  alphaIpwBal[m,] <- fitBal$alphaIpw
  alphaMsBal[m,] <- fitBal$alphaMs
  alphaFull[m,] <- fitBal$alphaFull
  
  betaCpBal[m,] <- fitBal$betaCp
  betaIpwBal[m,] <- fitBal$betaIpw
  betaMsBal[m,] <- fitBal$betaMs
  betaFull[m,] <- fitBal$betaFull
  
  time02 <- Sys.time()
  time02-time01
  
  trueAlpha <- fitBal$alphaFull
  trueBeta <- fitBal$betaFull
  
  
  
  
  
  ### completely SRS
  set.seed(M*(indMC-1)+m)
  
  indPhase2Rand <- sort(sample(1:N,n))
  indPhase1Rand <- (1:N)[-indPhase2Rand]
  
  indPhaseRand <- list(indPhase1Rand,indPhase2Rand)
  
  print('Rand')
  print(length(indPhase2Rand)/(length(indPhase1Rand)+length(indPhase2Rand)))
  
  
  
  time01 <- Sys.time()
  ### completely SRS
  fitRand <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseRand,tauCut,tol=0.001,link=link)
  
  alphaCpRand[m,] <- fitRand$alphaCp
  alphaMsRand[m,] <- fitRand$alphaMs
  alphaIpwRand[m,] <- fitRand$alphaIpw
  
  betaCpRand[m,] <- fitRand$betaCp
  betaMsRand[m,] <- fitRand$betaMs
  betaIpwRand[m,] <- fitRand$betaIpw
  
  time02 <- Sys.time()
  time02-time01
  
  
  
  
  
  
  
  
  
  
  
  time01 <- Sys.time()
  ### practical implementation for the second stage sampling (pilot + adaptive)
  set.seed(M*(indMC-1)+m)
  
  
  # balanced pilot
  pilotSam <- balSample(obsY,obsD,obsZ,r0,tauCut,pamela=TRUE)
  print('Pilot - Bal')
  
  pilotWt <- pilotSam$wt
  pilotTable <- pilotSam$tableStrat
  indStratPilot <- pilotSam$indStrat
  
  indPhasePilot <- pilotSam[c(1,2)] 
  
  indPhase1Pilot <- indPhasePilot[[1]]
  indPhase2Pilot <- indPhasePilot[[2]]
  
  print(length(indPhase2Pilot)/(length(indPhase1Pilot)+length(indPhase2Pilot)))
  
  
  # pilot sample
  fitPilot <- msFitNew(obsY,obsD,obsZ,obsX,indPhasePilot,tauCut,tol=0.001,wt=pilotWt,link=link)
  
  alphaCpPilot[m,] <- fitPilot$alphaCp
  alphaIpwPilot[m,] <- fitPilot$alphaIpw
  alphaMsPilot[m,] <- fitPilot$alphaMs
  
  betaCpPilot[m,] <- fitPilot$betaCp
  betaIpwPilot[m,] <- fitPilot$betaIpw
  betaMsPilot[m,] <- fitPilot$betaMs
  
  alphaTmp <- fitPilot$alphaIpw
  betaTmp <- fitPilot$betaIpw
  infoTmp <- fitPilot$infoIpw
  
  time02 <- Sys.time()
  time02-time01
  
  
  
  
  
  time01 <- Sys.time()
  ## optimal design
  designResult <- msDesignNew(obsY,obsD,obsZ,obsX,indPhasePilot,tauCut,r,alphaTmp,betaTmp,infoTmp,pilotTable,indStratPilot,link=link)
  time02 <- Sys.time()
  time02-time01
  
  
  time01 <- Sys.time()
  # optimal
  indSample <- data.frame(designResult$indSample)
  prob <- data.frame(designResult$prob)
  n2Ayz <- designResult$n2Ayz
  n2Byz <- designResult$n2Byz
  n2Byz0 <- designResult$n2Byz0
  n2ByzOpt <- designResult$n2ByzOpt
  
  apply(indSample,2,'sum')
  round(apply(prob,2,'sum'))
  sum(n2Ayz)
  apply(n2Byz,2,'sum')
  apply(n2Byz0,2,'sum')
  apply(n2ByzOpt,2,'sum')
  
  
  indPhase2Opt <- sort(c(indPhase2Pilot,indPhase1Pilot[which(indSample[,length(trueAlpha)+optVar]==1)]))
  indPhase1Opt <- indPhase1Pilot[which(indSample[,length(trueAlpha)+optVar]==0)]
  
  
  
  wt <- c()
  for (l in 1:length(indStratPilot)) {
    
    indOptTmp <- intersect(indStratPilot[[l]],indPhase2Opt)
    p2Opt <- length(indOptTmp)/length(indStratPilot[[l]])
    wt[indStratPilot[[l]]] <- 1/p2Opt
    
  }
  
  
  indPhaseOpt <- list(indPhase1Opt,indPhase2Opt)
  
  print('Opt')
  print(length(indPhase2Opt)/(length(indPhase1Opt)+length(indPhase2Opt)))
  
  fitOpt <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseOpt,tauCut,tol=0.001,wt,link=link)
  
  alphaCpOpt[m,] <- fitOpt$alphaCp
  alphaIpwOpt[m,] <- fitOpt$alphaIpw
  alphaMsOpt[m,] <- fitOpt$alphaMs
  
  betaCpOpt[m,] <- fitOpt$betaCp
  betaIpwOpt[m,] <- fitOpt$betaIpw
  betaMsOpt[m,] <- fitOpt$betaMs
  
  time02 <- Sys.time()
  time02-time01
  
  
  
  
  
  
  
  
  time01 <- Sys.time()
  ## Oracle design
  designOrac <- msDesignOracle(obsY,obsD,obsZ,obsX,tauCut,r,5000,M*(indMC-1)+m,fitAlphaBeta=list(trueAlpha,trueBeta),link=link)
  time02 <- Sys.time()
  time02-time01
  
  
  time01 <- Sys.time()
  set.seed(M*(indMC-1)+m)
  
  indSampleOrac <- data.frame(designOrac$indSample)
  prob <- data.frame(designOrac$prob)
  n2yz <- designOrac$n2yz
  n2yz0 <- designOrac$n2yz0
  indStratObs <- designOrac$indStratObs
  
  apply(indSampleOrac,2,'sum')
  round(apply(prob,2,'sum'))
  apply(n2yz,2,'sum')
  apply(n2yz0,2,'sum')
  
  
  indPhase2Orac <- sort(which(indSampleOrac[,length(trueAlpha)+optVar]==1))
  indPhase1Orac <- sort(which(indSampleOrac[,length(trueAlpha)+optVar]==0))
  
  probOrac <- prob[,length(trueAlpha)+optVar]
  probOrac <- probOrac/sum(probOrac)*(1-r)*N
  oracWt <- c()
  oracWt <- 1/probOrac
  oracWt[which(probOrac==0)] <- 0
  
  indPhaseOrac <- list(indPhase1Orac,indPhase2Orac)
  
  print('Orac')
  print(length(indPhase2Orac)/(length(indPhase1Orac)+length(indPhase2Orac)))
  
  
  fitOrac <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseOrac,tauCut,tol=0.001,oracWt,link=link)
  fitOrac$betaIpw
  fitOrac$betaMs
  
  
  alphaCpOrac[m,] <- fitOrac$alphaCp
  alphaIpwOrac[m,] <- fitOrac$alphaIpw
  alphaMsOrac[m,] <- fitOrac$alphaMs
  
  betaCpOrac[m,] <- fitOrac$betaCp
  betaIpwOrac[m,] <- fitOrac$betaIpw
  betaMsOrac[m,] <- fitOrac$betaMs
  
  time02 <- Sys.time()
  time02-time01
  
  
  
  time2 <- Sys.time()
  
  print(time2-time1)
  
  
}
timeMC2 <- Sys.time()

print(timeMC2 - timeMC1)




GetSqrtMse <- function(betaEst) {
  sqrtMse <- sqrt(apply(t(t(betaEst) - trueBeta)^2,2,'mean')) 
  return(sqrtMse)
}


GetBias <- function(betaEst) {
  bias <- apply(betaEst,2,'mean') - trueBeta
  return(bias)
}

print('beta - Cc Rand')
round(t(data.frame(sqrtMse=GetSqrtMse(betaCpRand),
                   bias=GetBias(betaCpRand),
                   sqrtVar=sqrt(GetSqrtMse(betaCpRand)^2-GetBias(betaCpRand)^2))),3)

print('beta - Ms Rand')
round(t(data.frame(sqrtMse=GetSqrtMse(betaMsRand),
                   bias=GetBias(betaMsRand),
                   sqrtVar=sqrt(GetSqrtMse(betaMsRand)^2-GetBias(betaMsRand)^2))),3)

print('beta - Ipw Bal')
round(t(data.frame(sqrtMse=GetSqrtMse(betaIpwBal),
                   bias=GetBias(betaIpwBal),
                   sqrtVar=sqrt(GetSqrtMse(betaIpwBal)^2-GetBias(betaIpwBal)^2))),3)

print('beta - Ms Opt')
round(t(data.frame(sqrtMse=GetSqrtMse(betaMsOpt),
                   bias=GetBias(betaMsOpt),
                   sqrtVar=sqrt(GetSqrtMse(betaMsOpt)^2-GetBias(betaMsOpt)^2))),3)

print('beta - Ms Orac')
round(t(data.frame(sqrtMse=GetSqrtMse(betaMsOrac),
                   bias=GetBias(betaMsOrac),
                   sqrtVar=sqrt(GetSqrtMse(betaMsOrac)^2-GetBias(betaMsOrac)^2))),3)


