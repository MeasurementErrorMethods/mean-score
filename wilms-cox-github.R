
library(survival)

nwts <- read.table('./nwts-expanded.txt',header=TRUE)
  

stage34 <- ifelse(nwts$stage%in%c(1,2),0,1)
nwts <- data.frame(nwts,stage34=stage34)

head(nwts)

X <- data.frame(nwts$instit,nwts$histol,nwts$stage34,nwts$age,nwts$tumdiam,nwts$histol*nwts$stage34)
colnames(X) <- c('instit','histol','stage34','age','tumdiam','histol:stage34')


dim(nwts)
with(nwts, table(instit,histol))
1-sum(nwts$instit==nwts$histol)/nrow(nwts)



# cox regression (by histol+stage34+age+tumdiam)
fitTrue <- coxph(Surv(nwts$trel,nwts$relaps,type='right')~X$histol+X$stage34+X$age+X$tumdiam+X$histol*X$stage34)
summary(fitTrue)

# cox regression (by instit+stage34+age+tumdiam)
fitErr <- coxph(Surv(nwts$trel,nwts$relaps,type='right')~X$instit+X$stage34+X$age+X$tumdiam+X$instit*X$stage34)
summary(fitErr)


  
indMC <- 1

M <- 10

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

coefCoxBal <- coefCoxRand <- coefCoxOpt <- coefCoxOrac <- matrix(nrow=M,ncol=length(coef(fitTrue)))


stratAll <- stratDesign(obsY,obsD,obsZ,tauCut)
strat <- c()
for (l in 1:length(stratAll)) {
  strat[stratAll[[l]]] <- l
}
strat <- factor(strat,levels=seq(1,length(stratAll)))

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
  
  
  ## application to coxph
  inSampleBal <- c()
  inSampleBal[indPhase1Bal] <- FALSE
  inSampleBal[indPhase2Bal] <- TRUE
  
  sdfWilmsBal <- data.frame(dfWilms[inSampleBal,],wt=balWt0[inSampleBal])
  coxphBal2phase <- coxph(Surv(trel,relaps,type='right')~histol+stage34+age+tumdiam+histol*stage34,weights=wt,data=sdfWilmsBal)
  
  coefCoxBal[m,] <- coef(coxphBal2phase)
  
  
  
  
  ### completely SRS
  set.seed(M*(indMC-1)+m)
  
  indPhase2Rand <- sort(sample(1:N,n))
  indPhase1Rand <- (1:N)[-indPhase2Rand]
  
  indPhaseRand <- list(indPhase1Rand,indPhase2Rand)
  
  print('Rand')
  print(length(indPhase2Rand)/(length(indPhase1Rand)+length(indPhase2Rand)))
  
  
  ## application to coxph
  inSampleRand <- c()
  inSampleRand[indPhase1Rand] <- FALSE
  inSampleRand[indPhase2Rand] <- TRUE
  
  
  randTable0 <- balTable0
  randTable0[,4] <- as.numeric(table(inSampleRand,strat)[2,])
  
  randWt0 <- c()
  for (l in 1:nrow(randTable0)) {
    
    randWt0[stratAll[[l]]] <- randTable0[l,5]/randTable0[l,4]
    
    if (randTable0[l,5]==0) {
      randWt0[stratAll[[l]]] <- 0
    } 
    
  }
  
  sdfWilmsRand <- data.frame(dfWilms[inSampleRand,],wt=randWt0[inSampleRand])
  coxphRand2phase <- coxph(Surv(trel,relaps,type='right')~histol+stage34+age+tumdiam+histol*stage34,weights=wt,data=sdfWilmsRand)
  
  coefCoxRand[m,] <- coef(coxphRand2phase)
  
  
  
  
  
  
  
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
  
  alphaTmp <- fitPilot$alphaIpw
  betaTmp <- fitPilot$betaIpw
  infoTmp <- fitPilot$infoIpw
  
  time02 <- Sys.time()
  time02-time01
  
  trueAlpha <- fitPilot$alphaFull
  trueBeta <- fitPilot$betaFull
  
  
  
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
  
  
  ## application to coxph
  inSampleOpt <- c()
  inSampleOpt[indPhase1Opt] <- FALSE
  inSampleOpt[indPhase2Opt] <- TRUE
  
  
  optTable0 <- balTable0
  optTable0[,4] <- as.numeric(table(inSampleOpt,strat)[2,])
  
  optWt0 <- c()
  for (l in 1:nrow(optTable0)) {
    
    optWt0[stratAll[[l]]] <- optTable0[l,5]/optTable0[l,4]
    
    if (optTable0[l,5]==0) {
      optWt0[stratAll[[l]]] <- 0
    } 
    
  }
  
  
  sdfWilmsOpt <- data.frame(dfWilms[inSampleOpt,],wt=optWt0[inSampleOpt])
  coxphOpt2phase <- coxph(Surv(trel,relaps,type='right')~histol+stage34+age+tumdiam+histol*stage34,weights=wt,data=sdfWilmsOpt)
  
  coefCoxOpt[m,] <- coef(coxphOpt2phase)
  
  
  
  
  
  
  
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
  
  
  ## application to coxph
  inSampleOrac <- c()
  inSampleOrac[indPhase1Orac] <- FALSE
  inSampleOrac[indPhase2Orac] <- TRUE
  
  
  
  oracTable0 <- balTable0
  oracTable0[,4] <- as.numeric(table(inSampleOrac,strat)[2,])
  
  oracWt0 <- c()
  for (l in 1:nrow(oracTable0)) {
    
    oracWt0[stratAll[[l]]] <- oracTable0[l,5]/oracTable0[l,4]
    
    if (oracTable0[l,5]==0) {
      oracWt0[stratAll[[l]]] <- 0
    } 
    
  }
  
  sdfWilmsOrac <- data.frame(dfWilms[inSampleOrac,],wt=oracWt[inSampleOrac])
  coxphOrac2phase <- coxph(Surv(trel,relaps,type='right')~histol+stage34+age+tumdiam+histol*stage34,weights=wt,data=sdfWilmsOrac)
  
  coefCoxOrac[m,] <- coef(coxphOrac2phase)
  
  
  
  
  
  time2 <- Sys.time()
  
  print(time2-time1)
  
  
}
timeMC2 <- Sys.time()

print(timeMC2 - timeMC1)






resultRand <- t(data.frame(sqrtMse=round(sqrt(apply(t(t(coefCoxRand)-coef(fitTrue))^2,2,'mean')),3),
                           bias=round(as.numeric(apply(coefCoxRand,2,'mean')-coef(fitTrue)),3),
                           sqrtVar=round(as.numeric(sqrt(apply(t(t(coefCoxRand)-coef(fitTrue))^2,2,'mean')-(apply(coefCoxRand,2,'mean')-coef(fitTrue))^2)),3)
))

resultBal <- t(data.frame(sqrtMse=round(sqrt(apply(t(t(coefCoxBal)-coef(fitTrue))^2,2,'mean')),3),
                          bias=round(as.numeric(apply(coefCoxBal,2,'mean')-coef(fitTrue)),3),
                          sqrtVar=round(as.numeric(sqrt(apply(t(t(coefCoxBal)-coef(fitTrue))^2,2,'mean')-(apply(coefCoxBal,2,'mean')-coef(fitTrue))^2)),3)
))

resultOpt <- t(data.frame(sqrtMse=round(sqrt(apply(t(t(coefCoxOpt)-coef(fitTrue))^2,2,'mean')),3),
                          bias=round(as.numeric(apply(coefCoxOpt,2,'mean')-coef(fitTrue)),3),
                          sqrtVar=round(as.numeric(sqrt(apply(t(t(coefCoxOpt)-coef(fitTrue))^2,2,'mean')-(apply(coefCoxOpt,2,'mean')-coef(fitTrue))^2)),3)
))

resultOrac <- t(data.frame(sqrtMse=round(sqrt(apply(t(t(coefCoxOrac)-coef(fitTrue))^2,2,'mean')),3),
                           bias=round(as.numeric(apply(coefCoxOrac,2,'mean')-coef(fitTrue)),3),
                           sqrtVar=round(as.numeric(sqrt(apply(t(t(coefCoxOrac)-coef(fitTrue))^2,2,'mean')-(apply(coefCoxOrac,2,'mean')-coef(fitTrue))^2)),3)
))




fitTrue

print('two-phase cox regression - Rand')
resultRand

print('two-phase cox regression - Bal')
resultBal

print('two-phase cox regression - Opt')
resultOpt

print('two-phase cox regression - Orac')
resultOrac


