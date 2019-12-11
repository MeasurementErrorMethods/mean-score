
indMC <- 1

M <- 10
cens <- 50

link <- 'cloglog'

N <- 4000; r <- 1-400/N; r0 <- 1 - 400/N/2
n <- (1-r)*N

source('./source-functions-adaptive-github.R')
  

#### global parameters

# discrete survival times
J <- 10
tau <- 1:J
indTauCut <- 6
tauCut <- tau[1:indTauCut]

# baseline hazards
baseLambda <- exp(seq(log(0.0325),log(0.95),length.out=J))
  
# regression coefficient
beta0 <- c(log(1.5),log(0.7)) # conti
beta1 <- c(log(1.3),-log(1.3))#/2 # discret

d0 <- length(beta0)
d1 <- length(beta1)

trueAlpha <- cloglog(baseLambda[1:indTauCut])
trueBeta <- c(beta0,beta1)

# censoring rate
pCen <- 1 - 1 # no random censoring


# r = 1 - validation rate
# r0: 1 - (phase II-A rate)
# r1: 1 - (phase II-B rate where phase II-A is excluded)
r1 <- 1-(r0-r)/r0
(1-r0)+r0*(1-r1)

# AR parameter for covariate correlation
rho <- 0.3

sig <- 0.1
para <- c(2,1.5,3,3)



### m-th MC loop

alphaCox <- alphaMsRand <- alphaCpRand <- alphaIpwRand <- alphaFull <- matrix(nrow=M,ncol=length(trueAlpha))
betaCox <- betaMsRand <- betaCpRand <- betaIpwRand <- betaFull <- matrix(nrow=M,ncol=length(trueBeta))

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


optVar <- 1
timeMC1 <- Sys.time()

for (m in 1:M) {
  
  time1 <- Sys.time()
  
  set.seed(2*M*(indMC-1)+2*m)
  
  print(paste('MC replication ',m,'/',M,sep=''))
  
  ## data generation
  
  simSam <- dataGen(N,r,tau,baseLambda,pCen,beta0,beta1,d0,d1,rho,link=link,sig=sig,para=para)
  
  obsY <- simSam$obsY
  obsD <- simSam$obsD
  obsZ <- simSam$obsZ
  obsX <- simSam$obsX
  indPhaseRand <- simSam$indPhase
  
  indPhase1Rand <- indPhaseRand[[1]]
  indPhase2Rand <- indPhaseRand[[2]]
  
  print('Rand')
  # print(length(indPhase2Rand)/(length(indPhase1Rand)+length(indPhase2Rand)))

  
  
  
  
  
  time01 <- Sys.time()
  ### completely balanced sampling
  set.seed(2*M*(indMC-1)+2*m)
  balSam0 <- balSample(obsY,obsD,obsZ,r,tauCut)
  
  balWt0 <- balSam0$wt
  balTable0 <- balSam0$tableStrat
  indStrat0 <- balSam0$indStrat
  
  indPhaseBal <- balSam0[c(1,2)]
  
  indPhase1Bal <- indPhaseBal[[1]]
  indPhase2Bal <- indPhaseBal[[2]]
  
  print('Bal')
  print(length(indPhase2Bal)/(length(indPhase1Bal)+length(indPhase2Bal)))
  
  
  
  
  # balanced sample
  fitBal <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseBal,tauCut,tol=0.0001,wt=balWt0,link=link)
  
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
  
  
  
  ### completely SRS
  set.seed(2*M*(indMC-1)+2*m)
  
  indPhase2Rand <- sort(sample(1:N,n))
  indPhase1Rand <- (1:N)[-indPhase2Rand]
  
  indPhaseRand <- list(indPhase1Rand,indPhase2Rand)
  
  print('Rand')
  print(length(indPhase2Rand)/(length(indPhase1Rand)+length(indPhase2Rand)))
  
  
  time01 <- Sys.time()
  ### completely SRS
  fitRand <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseRand,tauCut,tol=0.0001,link=link)
  
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
  set.seed(2*M*(indMC-1)+2*m)
  
  
  # balanced pilot
  pilotSam <- balSample(obsY,obsD,obsZ,r0,tauCut,pamela=FALSE)
  print('Pilot - Bal')
  
  pilotWt <- pilotSam$wt
  pilotTable <- pilotSam$tableStrat
  indStratPilot <- pilotSam$indStrat
  
  indPhasePilot <- pilotSam[c(1,2)] 
  
  indPhase1Pilot <- indPhasePilot[[1]]
  indPhase2Pilot <- indPhasePilot[[2]]
  
  print(length(indPhase2Pilot)/(length(indPhase1Pilot)+length(indPhase2Pilot)))
  
  
  # pilot sample
  fitPilot <- msFitNew(obsY,obsD,obsZ,obsX,indPhasePilot,tauCut,tol=0.0001,wt=pilotWt,link=link)
  
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
  
  
  fitOpt <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseOpt,tauCut,tol=0.0001,wt,link=link)
  
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
  designOrac <- msDesignOracle(obsY,obsD,obsZ,obsX,tauCut,r,5000,M*(indMC-1)+m,fitAlphaBeta=list(trueAlpha,trueBeta),link=link,sig=sig,para=para)
  time02 <- Sys.time()
  time02-time01
  
  
  time01 <- Sys.time()
  set.seed(2*M*(indMC-1)+2*m)
  
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
  
  
  
  fitOrac <- msFitNew(obsY,obsD,obsZ,obsX,indPhaseOrac,tauCut,tol=0.0001,oracWt,link=link)
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




getSqrtMse <- function(betaEst) {
  
  sqrtMse <- sqrt(apply(t(t(betaEst) - trueBeta)^2,2,'mean')) 
  return(sqrtMse)
}


getBias <- function(betaEst) {
  
  bias <- apply(betaEst,2,'mean') - trueBeta
  return(bias)
}

getSqrtVar <- function(betaEst) {
  
  sqrtMse <- getSqrtMse(betaEst)
  bias <- getBias(betaEst)
  
  sqrtVar <- sqrt(sqrtMse^2 - bias^2)
  return(sqrtVar)
}



print('beta - Full')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaFull),
                         bias=getBias(betaFull),
                         sqrtVar=getSqrtVar(betaFull))),3))

print('beta - Cc Rand')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaCpRand),
                         bias=getBias(betaCpRand),
                         sqrtVar=getSqrtVar(betaCpRand))),3))

print('beta - Ms Rand')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaMsRand),
                         bias=getBias(betaMsRand),
                         sqrtVar=getSqrtVar(betaMsRand))),3))

print('beta - Ipw Bal')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaIpwBal),
                         bias=getBias(betaIpwBal),
                         sqrtVar=getSqrtVar(betaIpwBal))),3))

print('beta - Ms Opt')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaMsOpt),
                         bias=getBias(betaMsOpt),
                         sqrtVar=getSqrtVar(betaMsOpt))),3))

print('beta - Ms Orac')
print(round(t(data.frame(sqrtMse=getSqrtMse(betaMsOrac),
                         bias=getBias(betaMsOrac),
                         sqrtVar=getSqrtVar(betaMsOrac))),3))



