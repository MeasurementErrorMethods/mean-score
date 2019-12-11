# rm(list=ls(all=TRUE))

# library(survival)

library(MASS)


### logit transformation
logit <- function(t) log(t/(1-t))
invLogit <- function(t) exp(t)/(1+exp(t))

### cloglog transformation
cloglog <- function(t) log(-log(1-t))
invCll <- function(t) 1-exp(-exp(t))


### autoregressive covariance generator
sqrtArCov <- function(rho,d) {
  V <- matrix(nrow=d,ncol=d)
  
  for (i in 1:d) {
    for (j in 1:d) {
      V[i,j] <- rho^(abs(i-j))
    } 
  }
  
  eigenS <- eigen(V)
  
  P <- eigenS$vectors
  R <- eigenS$values
  
  sqrtV <- P%*%diag(sqrt(R))%*%t(P)
  
  return(sqrtV)
}


### conditional proportional hazards from logit hazards
condPH <- function(baseLambda,X,beta,link) {
  
  if (link=='logit') {
    alpha <- logit(baseLambda)
  } else if (link=='cloglog') {
    alpha <- cloglog(baseLambda)
  }
  betaX <- X%*%beta
  
  tmpMat <- matrix(nrow=nrow(X),ncol=length(baseLambda))
  for (j in 1:length(baseLambda)) {
    tmpMat[,j] <- alpha[j] + betaX
  }
  
  if (link=='logit') {
    lambdaMat <- invLogit(tmpMat)
  } else if (link=='cloglog') {
    lambdaMat <- invCll(tmpMat)
  }
  
  return(lambdaMat)
}


### discrete survival probability
pSurv <- function(lambda,X=NULL,beta=NULL,link) {
  
  if (is.null(X)==TRUE || is.null(beta)==TRUE) {
    tmpSurv <- c()
    
    tmpSurv[1] <- 1
    for (j in 2:length(lambda)) {
      tmpSurv[j] <- prod(1-lambda[1:(j-1)])
    }
    
    return(tmpSurv)
    
  } else {
    condLambda <- condPH(lambda,X,beta,link=link)
    
    tmpSurv <- matrix(nrow=nrow(X),ncol=length(lambda))
    for (i in 1:nrow(X)) {
      tmpSurv[i,1] <- 1
      for (j in 2:length(lambda)) {
        tmpSurv[i,j] <- prod(1-condLambda[i,1:(j-1)])
      }
    } 
    
    return(tmpSurv)
  }
}







### data generation
dataGen <- function(N,r,tau,baseLambda,pCen,beta0,beta1,d0,d1,rho,dimAux=NULL,link,sig=0.1,cutVal=c(-Inf,0.25,0.5,0.75,Inf),para=c(2,3,3,3)) {
  
  if (is.null(dimAux)) {
    dimAux <- 1
  }
  
  trueBeta <- c(beta0,beta1)
  
  sqrtV <- sqrtArCov(rho,length(trueBeta))
  W <- matrix(rnorm(N*(length(trueBeta))),nrow=N,ncol=(length(trueBeta)))%*%sqrtV
  obsX <- pnorm(W)
  
  obsX[,1] <- qbeta(obsX[,1],para[1],para[2])
  obsX[,2] <- qbeta(obsX[,2],para[3],para[4])
  
  obsA <- obsX[,1:dimAux] + matrix((rnorm(dimAux*N,0,sig)),nrow=N,ncol=dimAux)
  
  obsZ <- matrix(nrow=N,ncol=dimAux)
  cutA <- cutVal
  for (j in 1:ncol(obsZ)) {
    for (k in 1:length(cutA)) {
      ind <- which(obsA[,j]>cutA[k-1] & obsA[,j]<cutA[k])
      obsZ[ind,j] <- k-1
    }
  }
  
  obsX0 <- obsX[,1:d0]
  obsX1 <- obsX[,(d0+1):(d0+d1)]
  for (i in 1:N) {
    for (j in 1:ncol(obsX1)) {
      obsX1[i,j] <- sample(c(1,0),1,prob=c(obsX1[i,j],1-obsX1[i,j]))
    }
  }
  
  obsX <- cbind(obsX0,obsX1)
  
  # round(cor(cbind(obsX,obsZ)),2)
  
  # conditional hazards
  condLambda <- condPH(baseLambda,obsX,trueBeta,link=link)
  
  # conditional survival probability
  condSurv <- pSurv(baseLambda,obsX,trueBeta,link=link)
  baseSurv <- pSurv(baseLambda,link=link)
  
  # conditional probability mass
  # condDens <- data.frame(condLambda*condSurv)
  tmp <- condSurv[,-ncol(condSurv)]-condSurv[,-1]
  tmp <- cbind(tmp,1-apply(tmp,1,'sum'))
  condDens <- tmp
  names(condDens) <- tau
  baseDens <- baseLambda*baseSurv
  
  
  # survival time generation
  obsY <- c()
  for (i in 1:N) {
    obsY[i] <- sample(tau,1,prob=condDens[i,])
  }
  
  
  obsD <- sample(c(1,0),N,replace=TRUE,prob=c(1-pCen,pCen))
  # indOver <- which(obsY>max(tauCut))
  # obsD[indOver] <- 0
  
  indPhase1 <- sort(sample(1:length(obsY),length(obsY)*r))
  indPhase2 <- (1:N)[-indPhase1]
  indPhase <- list(indPhase1,indPhase2)
  
  return(list(obsY=obsY,
              obsD=obsD,
              obsZ=obsZ,
              obsX=obsX,
              indPhase=indPhase))
}


stratDesign <- function(obsY,obsD,obsZ,tauCut) {
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- list()
  tableStrat <- c()
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        indStrat[[l]] <- tmp
        l <- l+1
      }
    }
  }
  
  return(indStrat)
  
}



balSample <- function(obsY,obsD,obsZ,r0,tauCut,pamela=FALSE) {
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- list()
  tableStrat <- c()
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tableStrat <- rbind(tableStrat,c(i,j,k,length(tmp)))
        indStrat[[l]] <- tmp
        l <- l+1
      }
    }
  }
  
  nonZeroStrat <- which(tableStrat[,4]!=0)
  numSam <- round(length(obsY)*(1-r0))
  numSamStrat <- numSam/length(nonZeroStrat)
  
  # nFloor <- floor(numSamStrat)
  # nCeiling <- ceiling(numSamStrat)
  # pSam <- nCeiling - numSamStrat
  # tmpSam <- sample(c(nFloor,nCeiling),length(nonZeroStrat),replace=TRUE,prob=c(pSam,1-pSam))
  
  indPhase2Bal <- c()
  openStrat <- c()
  w <- c()
  for (l in 1:length(nonZeroStrat)) {
    tmp <- which(apply(obsMat,1,'identical',y=tableStrat[nonZeroStrat[l],1:3])==TRUE)
    
    # if (length(tmp)<tmpSam[l]) {
    if (length(tmp)<numSamStrat) {
      
      if (pamela==FALSE) {
        indPhase2Bal <- c(indPhase2Bal,tmp)
        w[tmp] <- 1
        openStrat[nonZeroStrat[l]] <- 0
        # w[tmp] <- length(tmp)/length(obsY)
      } else {
        if (tableStrat[nonZeroStrat[l],1]<max(tauCut) & tableStrat[nonZeroStrat[l],2]==0) {
          if (length(tmp)<=1) {
            indPhase2Bal <- c(indPhase2Bal,tmp)
          } else {
            indPhase2Bal <- c(indPhase2Bal,sample(tmp,min(numSamStrat/2,length(tmp))))
          }
          w[tmp] <- min(numSamStrat/2,length(tmp))/length(tmp)
          openStrat[nonZeroStrat[l]] <- 0
        } else {
          indPhase2Bal <- c(indPhase2Bal,tmp)
          w[tmp] <- 1
          openStrat[nonZeroStrat[l]] <- 0
        }
      }
      
    } else {
      
      if (pamela==FALSE) {
        indPhase2Bal <- c(indPhase2Bal,sample(tmp,numSamStrat))
        w[tmp] <- numSamStrat/length(tmp)
        openStrat[nonZeroStrat[l]] <- 1
      } else {
        if (tableStrat[nonZeroStrat[l],1]<max(tauCut) & tableStrat[nonZeroStrat[l],2]==0) {
          if (length(tmp)<=1) {
            indPhase2Bal <- c(indPhase2Bal,tmp)
          } else {
            indPhase2Bal <- c(indPhase2Bal,sample(tmp,numSamStrat/2))
          }
          w[tmp] <- numSamStrat/2/length(tmp)
          openStrat[nonZeroStrat[l]] <- 0
        } else {
          indPhase2Bal <- c(indPhase2Bal,sample(tmp,numSamStrat))
          w[tmp] <- numSamStrat/length(tmp)
          openStrat[nonZeroStrat[l]] <- 1
        }
      }
      # w[tmp] <- numSamStrat/length(obsY)
      
      # indPhase2Bal <- c(indPhase2Bal,sample(tmp,tmpSam[l]))
      # w[tmp] <- tmpSam[l]/length(tmp)
      
    }  
  }
  
  indPhase2Bal <- sort(indPhase2Bal)
  indPhase1Bal <- (1:length(obsY))[-indPhase2Bal]
  
  length(indPhase2Bal)
  
  indOpen <- which(openStrat==1)
  numSamOpenStrat <- (numSam-length(indPhase2Bal))/length(indOpen)
  for (l in 1:length(indOpen)) {
    
    tmp <- which(apply(obsMat[indPhase1Bal,],1,'identical',y=tableStrat[indOpen[l],1:3])==TRUE)
    
    if (length(tmp)<numSamOpenStrat) {
      
      # print(c(indOpen[l],length(tmp),length(tmp)))
      
      indPhase2Bal <- c(indPhase2Bal,indPhase1Bal[tmp])
      
    } else {
      
      # print(c(indOpen[l],length(tmp),numSamOpenStrat))
      
      indPhase2Bal <- c(indPhase2Bal,sample(indPhase1Bal[tmp],numSamOpenStrat))
      
    }
    
  }
  
  indPhase2Bal <- sort(indPhase2Bal)
  indPhase1Bal <- (1:length(obsY))[-indPhase2Bal]
  
  length(indPhase2Bal)
  
  # randomly distribute remaining slots for validation 
  if (length(indPhase2Bal)<numSam) {
    
    if (pamela==FALSE) {
      numRes <- numSam-length(indPhase2Bal)
      ind0 <- sample(indPhase1Bal,numRes,prob=w[indPhase1Bal])
      
      # ind0 <- order(w0[indPhase1Bal])[1:numRes]
      # ind0 <- indPhase1Bal[ind0]
      
      indPhase2Bal <- sort(c(indPhase2Bal,ind0))
      indPhase1Bal <- (1:length(obsY))[-indPhase2Bal]
    } else {
      tmp <- which(obsMat[indPhase1Bal,1]<max(tauCut) & obsMat[indPhase1Bal,2]==0)
      
      numRes <- numSam-length(indPhase2Bal)
      ind0 <- sample(indPhase1Bal[-tmp],numRes,prob=w[indPhase1Bal[-tmp]])
      
      # ind0 <- order(w0[indPhase1Bal])[1:numRes]
      # ind0 <- indPhase1Bal[ind0]
      
      indPhase2Bal <- sort(c(indPhase2Bal,ind0))
      indPhase1Bal <- (1:length(obsY))[-indPhase2Bal]
    }
    
  }
  length(indPhase2Bal)
  
  obsFreq <- c()
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tmp1 <- which(apply(obsMat[indPhase2Bal,],1,'identical',y=c(i,j,k))==TRUE)
        obsFreq <- c(obsFreq,length(tmp1))
        
        w[tmp] <- length(tmp1)/length(tmp)
        # w[tmp] <- length(tmp0)/length(obsY)
        
      }
    }
  }
  
  # w[which(w<0.01)] <- 0.01
  
  tableStrat <- cbind(tableStrat,obsFreq)
  colnames(tableStrat) <- c('Y','D','Z','sizeStrat','sizeSample')
  tableStrat <- tableStrat[,c(1,2,3,5,4)]
  
  return(list(indPhase1Bal=indPhase1Bal,
              indPhase2Bal=indPhase2Bal,
              wt=1/w,
              tableStrat=tableStrat,
              indStrat=indStrat))
}





randSample <- function(obsY,obsD,obsZ,r0,tauCut) {
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- list()
  tableStrat <- c()
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tableStrat <- rbind(tableStrat,c(i,j,k,length(tmp)))
        indStrat[[l]] <- tmp
        l <- l+1
      }
    }
  }
  
  nonZeroStrat <- which(tableStrat[,4]!=0)
  
  indPhase2Pilot <- sort(sample(1:N,N*(1-r0)))
  indPhase1Pilot <- (1:length(obsY))[-indPhase2Pilot]
  
  length(indPhase2Pilot)
  
  obsFreq <- c()
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tmp1 <- which(apply(obsMat[indPhase2Pilot,],1,'identical',y=c(i,j,k))==TRUE)
        obsFreq <- c(obsFreq,length(tmp1))
        
      }
    }
  }
  
  # w[which(w<0.01)] <- 0.01
  
  tableStrat <- cbind(tableStrat,obsFreq)
  colnames(tableStrat) <- c('Y','D','Z','sizeStrat','sizeSample')
  tableStrat <- tableStrat[,c(1,2,3,5,4)]
  
  return(list(indPhase1Pilot=indPhase1Pilot,
              indPhase2Pilot=indPhase2Pilot,
              wt=rep(1,N),
              tableStrat=tableStrat,
              indStrat=indStrat))
}






msFitNew <- function(obsY,obsD,obsZ,obsX,indPhase,tauCut,tol=0.001,wt=NULL,link) {
  
  N <- length(obsY)
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  if (is.null(wt)) {
    
    wt <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
    wt <- rep(wt,length(obsY))
    
    # wt[indPhase1] <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
    # wt[indPhase2] <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
  }
  
  valWt <- wt[indPhase2]
  nonvalWt <- wt[indPhase1]
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  phase1Mat <- obsMat[indPhase1,]
  phase2Mat <- obsMat[indPhase2,]
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- indStrat1 <- indStrat2 <- list()
  tableStrat <- c()
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        # tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tmp1 <- which(apply(phase1Mat,1,'identical',y=c(i,j,k))==TRUE)
        tmp2 <- which(apply(phase2Mat,1,'identical',y=c(i,j,k))==TRUE)
        
        tableStrat <- rbind(tableStrat,c(i,j,k,length(tmp1),length(tmp2)))
        
        indStrat1[[l]] <- indPhase1[tmp1]
        indStrat2[[l]] <- indPhase2[tmp2]
        # indStrat[[l]] <- tmp
        l <- l+1
      }
    }
  }
  nonZeroStrat <- which((tableStrat[,4]+tableStrat[,5])!=0)
  
  # indStrat <- indStrat[nonZeroStrat]
  indStrat1 <- indStrat1[nonZeroStrat]
  indStrat2 <- indStrat2[nonZeroStrat]
  
  tableStrat <- tableStrat[nonZeroStrat,]
  
  
  valX <- as.matrix(obsX[indPhase2,])
  nonvalX <- as.matrix(obsX[indPhase1,])
  
  valZ <- as.matrix(obsZ[indPhase2,])
  nonvalZ <- as.matrix(obsZ[indPhase1,])
  
  valY <- ceiling(obsY[indPhase2])
  nonvalY <- ceiling(obsY[indPhase1])
  
  # pseudo outervations
  D <- matrix(0,nrow=length(obsY),ncol=length(tauCut))
  for (i in 1:length(obsY)) {
    for (j in 1:length(tauCut)) {
      if ((tmpY[i]) > (tauCut[j]-1) & (tmpY[i]) <= tauCut[j] & tmpD[i]==1) {
        D[i,j] <- 1
      }
    }
  }
  colnames(D) <- tauCut
  
  # cbind(valY,D[indPhase2,],obsY[indPhase2])[1:50,]
  
  valY[which(valY>max(tauCut))] <- max(tauCut)
  nonvalY[which(nonvalY>max(tauCut))] <- max(tauCut)
  
  valD <- D[indPhase2,]
  nonvalD <- D[indPhase1,]
  
  # cbind(valY,valD,obsY[indPhase2],obsD[indPhase2])[1:50,]
  
  combY <- c(valY,nonvalY)
  combD <- rbind(valD,nonvalD)
  combZ <- rbind(valZ,nonvalZ)
  combX <- rbind(valX,nonvalX)
  combWt <- c(valWt,nonvalWt)
  
  combY[indPhase1] <- nonvalY
  combY[indPhase2] <- valY
  
  combD[indPhase1,] <- nonvalD
  combD[indPhase2,] <- valD
  
  combZ[indPhase1,] <- nonvalZ
  combZ[indPhase2,] <- valZ
  
  combX[indPhase1,] <- nonvalX
  combX[indPhase2,] <- valX
  
  combWt[indPhase1] <- nonvalWt
  combWt[indPhase2] <- valWt
  
  
  
  
  
  ## new version
  
  # initial estimates
  msAlpha <- cpAlpha <- ipwAlpha <- fullAlpha <- rep(0,length(tauCut))
  msBeta <- cpBeta <- ipwBeta <- fullBeta <- rep(0,ncol(obsX))
  
  coefMs <- coefCp <- coefIpw <- coefFull<- c(msAlpha,msBeta)
  
  eps <- 10
  iter <- 1
  while (eps>tol) {
    eps0 <- eps
    # print(paste('Newton-Rhapson iteration ',iter,sep=''))
    
    tmp1Alpha <- coefMs[1:length(tauCut)]
    tmp1Beta <- coefMs[-(1:length(tauCut))]
    
    tmp2Alpha <- coefCp[1:length(tauCut)]
    tmp2Beta <- coefCp[-(1:length(tauCut))]
    
    tmp3Alpha <- coefFull[1:length(tauCut)]
    tmp3Beta <- coefFull[-(1:length(tauCut))]
    
    tmp4Alpha <- coefIpw[1:length(tauCut)]
    tmp4Beta <- coefIpw[-(1:length(tauCut))]
    
    
    ## score and information
    scoreA <- scoreACp <- scoreAIpw <- scoreAFull <- rep(0,length(tauCut))
    scoreB <- scoreBCp <- scoreBIpw <- scoreBFull <- rep(0,ncol(valX))
    
    hessAA <- hessAACp <- hessAAIpw <- hessAAFull <- matrix(0,nrow=length(tauCut),ncol=length(tauCut))
    hessAB <- hessABCp <- hessABIpw <- hessABFull <- matrix(0,nrow=length(tauCut),ncol=ncol(valX))
    hessBB <- hessBBCp <- hessBBIpw <- hessBBFull <- matrix(0,nrow=ncol(valX),ncol=ncol(valX))
    
    
    for (l in 1:nrow(tableStrat)) {
      
      I1 <- tableStrat[l,4]
      I2 <- tableStrat[l,5]
      
      for (i in indStrat2[[l]]) {
        
        for (j in 1:min(combY[i],max(tauCut))) {
          
          if (link=='logit') {
            
            mu1 <- invLogit(tmp1Alpha[j]+c(t(tmp1Beta)%*%combX[i,]))
            mu2 <- invLogit(tmp2Alpha[j]+c(t(tmp2Beta)%*%combX[i,]))
            mu3 <- invLogit(tmp3Alpha[j]+c(t(tmp3Beta)%*%combX[i,]))
            mu4 <- invLogit(tmp4Alpha[j]+c(t(tmp4Beta)%*%combX[i,]))
            
            # mean score (phase 2 + auxiliary phase 1)
            scoreA[j] <- scoreA[j] + (1+I1/I2)*(combD[i,j] - mu1)
            hessAA[j,j] <- hessAA[j,j] - (1+I1/I2)*mu1*(1-mu1)
            hessAB[j,] <- hessAB[j,] - (1+I1/I2)*mu1*(1-mu1)*combX[i,]
            
            scoreB <- scoreB + (1+I1/I2)*(combD[i,j] - mu1)*combX[i,]
            hessBB <- hessBB - (1+I1/I2)*mu1*(1-mu1)*combX[i,]%*%t(combX[i,])
            
            # complete (phase 2 only)
            scoreACp[j] <- scoreACp[j] + (combD[i,j] - mu2)
            hessAACp[j,j] <- hessAACp[j,j] - mu2*(1-mu2)
            hessABCp[j,] <- hessABCp[j,] - mu2*(1-mu2)*combX[i,]
            
            scoreBCp <- scoreBCp + (combD[i,j] - mu2)*combX[i,]
            hessBBCp <- hessBBCp - mu2*(1-mu2)*combX[i,]%*%t(combX[i,])
            
            # full cohorts
            scoreAFull[j] <- scoreAFull[j] + (combD[i,j] - mu3)
            hessAAFull[j,j] <- hessAAFull[j,j] - mu3*(1-mu3)
            hessABFull[j,] <- hessABFull[j,] - mu3*(1-mu3)*combX[i,]
            
            scoreBFull <- scoreBFull + (combD[i,j] - mu3)*combX[i,]
            hessBBFull <- hessBBFull - mu3*(1-mu3)*combX[i,]%*%t(combX[i,])
            
            # IPW
            scoreAIpw[j] <- scoreAIpw[j] + (combD[i,j] - mu4)*combWt[i]
            hessAAIpw[j,j] <- hessAAIpw[j,j] - mu4*(1-mu4)*combWt[i]
            hessABIpw[j,] <- hessABIpw[j,] - mu4*(1-mu4)*combX[i,]*combWt[i]
            
            scoreBIpw <- scoreBIpw + (combD[i,j] - mu4)*combX[i,]*combWt[i]
            hessBBIpw <- hessBBIpw - mu4*(1-mu4)*combX[i,]%*%t(combX[i,])*combWt[i] 
            
          } else if (link=='cloglog') {
            
            cll1 <- tmp1Alpha[j]+c(t(tmp1Beta)%*%combX[i,])
            cll2 <- tmp2Alpha[j]+c(t(tmp2Beta)%*%combX[i,])
            cll3 <- tmp3Alpha[j]+c(t(tmp3Beta)%*%combX[i,])
            cll4 <- tmp4Alpha[j]+c(t(tmp4Beta)%*%combX[i,]) 
            
            mu1 <- invCll(cll1)
            mu2 <- invCll(cll2)
            mu3 <- invCll(cll3)
            mu4 <- invCll(cll4)
            
            # mean score (phase 2 + auxiliary phase 1)
            tmp1 <- combD[i,j]*exp(cll1)/mu1 - exp(cll1)
            tmp2 <- combD[i,j]*exp(cll1)/mu1*(1-exp(cll1-exp(cll1))/mu1) - exp(cll1)
            
            scoreA[j] <- scoreA[j] + (1+I1/I2)*tmp1
            hessAA[j,j] <- hessAA[j,j] + (1+I1/I2)*tmp2
            hessAB[j,] <- hessAB[j,] + (1+I1/I2)*tmp2*combX[i,]
            
            scoreB <- scoreB + (1+I1/I2)*tmp1*combX[i,]
            hessBB <- hessBB + (1+I1/I2)*tmp2*combX[i,]%*%t(combX[i,])
            
            # complete (phase 2 only)
            tmp1 <- combD[i,j]*exp(cll2)/mu2 - exp(cll2)
            tmp2 <- combD[i,j]*exp(cll2)/mu2*(1-exp(cll2-exp(cll2))/mu2) - exp(cll2)
            
            scoreACp[j] <- scoreACp[j] + tmp1
            hessAACp[j,j] <- hessAACp[j,j] + tmp2
            hessABCp[j,] <- hessABCp[j,] + tmp2*combX[i,]
            
            scoreBCp <- scoreBCp + tmp1*combX[i,]
            hessBBCp <- hessBBCp + tmp2*combX[i,]%*%t(combX[i,])
            
            # full cohorts
            tmp1 <- combD[i,j]*exp(cll3)/mu3 - exp(cll3)
            tmp2 <- combD[i,j]*exp(cll3)/mu3*(1-exp(cll3-exp(cll3))/mu3) - exp(cll3)
            
            scoreAFull[j] <- scoreAFull[j] + tmp1
            hessAAFull[j,j] <- hessAAFull[j,j] + tmp2
            hessABFull[j,] <- hessABFull[j,] + tmp2*combX[i,]
            
            scoreBFull <- scoreBFull + tmp1*combX[i,]
            hessBBFull <- hessBBFull + tmp2*combX[i,]%*%t(combX[i,])
            
            # IPW
            tmp1 <- combD[i,j]*exp(cll4)/mu4 - exp(cll4)
            tmp2 <- combD[i,j]*exp(cll4)/mu4*(1-exp(cll4-exp(cll4))/mu4) - exp(cll4)
            
            scoreAIpw[j] <- scoreAIpw[j] + tmp1*combWt[i]
            hessAAIpw[j,j] <- hessAAIpw[j,j] + tmp2*combWt[i]
            hessABIpw[j,] <- hessABIpw[j,] + tmp2*combX[i,]*combWt[i]
            
            scoreBIpw <- scoreBIpw + tmp1*combX[i,]*combWt[i]
            hessBBIpw <- hessBBIpw + tmp2*combX[i,]%*%t(combX[i,])*combWt[i] 
            
          }
          
          
        }
      }
      
      # extra updates for full cohorts
      for (i in indStrat1[[l]]) {
        
        for (j in 1:min(combY[i],max(tauCut))) {
          
          if (link=='logit') {
            
            mu3 <- invLogit(tmp3Alpha[j]+c(t(tmp3Beta)%*%combX[i,]))
            
            # full cohorts
            scoreAFull[j] <- scoreAFull[j] + (combD[i,j] - mu3)
            hessAAFull[j,j] <- hessAAFull[j,j] - mu3*(1-mu3)
            hessABFull[j,] <- hessABFull[j,] - mu3*(1-mu3)*combX[i,]
            
            scoreBFull <- scoreBFull + (combD[i,j] - mu3)*combX[i,]
            hessBBFull <- hessBBFull - mu3*(1-mu3)*combX[i,]%*%t(combX[i,])
            
          } else if (link=='cloglog') {
            
            cll3 <- tmp3Alpha[j]+c(t(tmp3Beta)%*%combX[i,])
            mu3 <- invCll(cll3)
            
            # full cohorts
            tmp1 <- combD[i,j]*exp(cll3)/mu3 - exp(cll3)
            tmp2 <- combD[i,j]*exp(cll3)/mu3*(1-exp(cll3-exp(cll3))/mu3) - exp(cll3)
            
            scoreAFull[j] <- scoreAFull[j] + tmp1
            hessAAFull[j,j] <- hessAAFull[j,j] + tmp2
            hessABFull[j,] <- hessABFull[j,] + tmp2*combX[i,]
            
            scoreBFull <- scoreBFull + tmp1*combX[i,]
            hessBBFull <- hessBBFull + tmp2*combX[i,]%*%t(combX[i,])
            
          }
        }
      }
      
    }
    
    # mean score
    score <- c(scoreA,scoreB)
    
    hess1 <- cbind(hessAA,hessAB)
    hess2 <- cbind(t(hessAB),hessBB)
    hess <- rbind(hess1,hess2)
    
    # complete
    scoreCp <- c(scoreACp,scoreBCp)
    
    hess1Cp <- cbind(hessAACp,hessABCp)
    hess2Cp <- cbind(t(hessABCp),hessBBCp)
    hessCp <- rbind(hess1Cp,hess2Cp)
    
    # full cohorts
    scoreFull <- c(scoreAFull,scoreBFull)
    
    hess1Full <- cbind(hessAAFull,hessABFull)
    hess2Full <- cbind(t(hessABFull),hessBBFull)
    hessFull <- rbind(hess1Full,hess2Full)
    
    # IPW
    scoreIpw <- c(scoreAIpw,scoreBIpw)
    
    hess1Ipw <- cbind(hessAAIpw,hessABIpw)
    hess2Ipw <- cbind(t(hessABIpw),hessBBIpw)
    hessIpw <- rbind(hess1Ipw,hess2Ipw)
    
    
    # Newton-Raphson
    coefMs <- coefMs - ginv(hess)%*%score
    coefCp <- coefCp - ginv(hessCp)%*%scoreCp
    coefFull <- coefFull - ginv(hessFull)%*%scoreFull
    coefIpw <- coefIpw - ginv(hessIpw)%*%scoreIpw
    
    epsMs <- max(abs(coefMs-c(tmp1Alpha,tmp1Beta)))
    epsCp <- max(abs(coefCp-c(tmp2Alpha,tmp2Beta)))
    epsFull <- max(abs(coefFull-c(tmp3Alpha,tmp3Beta)))
    epsIpw <- max(abs(coefIpw-c(tmp4Alpha,tmp4Beta)))
    
    
    epsAll <- c(epsMs,epsCp,epsIpw,epsFull)
    eps <- max(epsAll)
    
    iter <- iter+1
    # print(paste('update eps: ',round(epsAll,3),sep=''))
    
    # if (eps-eps0 > 0) {
    if (iter > 20) {
      print('The algorithm may not converge')
      break()
    }
    
  }
  
  
  
  # return(list(alphaMs=tmp1Alpha, betaMs=tmp1Beta,
  #             alphaCp=tmp2Alpha, betaCp=tmp2Beta,
  #             alphaFull=tmp3Alpha, betaFull=tmp3Beta,
  #             alphaIpw=tmp4Alpha, betaIpw=tmp4Beta,
  #             infoFull=ginv(-hessFull/length(obsY)),
  #             infoMs=ginv(-hess/length(obsY)),
  #             infoIpw=ginv(-hessIpw/length(obsY)),
  #             infoCp=ginv(-hessCp/length(valY))))
  
  return(list(alphaMs=tmp1Alpha, betaMs=tmp1Beta,
              alphaCp=tmp2Alpha, betaCp=tmp2Beta,
              alphaFull=tmp3Alpha, betaFull=tmp3Beta,
              alphaIpw=tmp4Alpha, betaIpw=tmp4Beta,
              infoFull=-hessFull/length(obsY),
              infoMs  =-hess/length(obsY),
              infoIpw =-hessIpw/length(obsY),
              infoCp  =-hessCp/length(valY)))
  
  
  
  
  # ## old version
  # 
  # # initial estimates
  # msAlpha <- cpAlpha <- ipwAlpha <- fullAlpha <- rep(0,length(tauCut))
  # msBeta <- cpBeta <- ipwBeta <- fullBeta <- rep(0,ncol(obsX))
  # 
  # coefMs <- coefCp <- coefIpw <- coefFull<- c(msAlpha,msBeta)
  # 
  # eps <- 10
  # iter <- 1
  # while (eps>tol) {
  #   eps0 <- eps
  #   # print(paste('Newton-Rhapson iteration ',iter,sep=''))
  #   
  #   tmp1Alpha <- coefMs[1:length(tauCut)]
  #   tmp1Beta <- coefMs[-(1:length(tauCut))]
  #   
  #   tmp2Alpha <- coefCp[1:length(tauCut)]
  #   tmp2Beta <- coefCp[-(1:length(tauCut))]
  #   
  #   tmp3Alpha <- coefFull[1:length(tauCut)]
  #   tmp3Beta <- coefFull[-(1:length(tauCut))]
  #   
  #   tmp4Alpha <- coefIpw[1:length(tauCut)]
  #   tmp4Beta <- coefIpw[-(1:length(tauCut))]
  #   
  #   
  #   ## score and information
  #   scoreA <- scoreACp <- scoreAIpw <- scoreAFull <- rep(0,length(tauCut))
  #   scoreB <- scoreBCp <- scoreBIpw <- scoreBFull <- rep(0,ncol(valX))
  #   
  #   hessAA <- hessAACp <- hessAAIpw <- hessAAFull <- matrix(0,nrow=length(tauCut),ncol=length(tauCut))
  #   hessAB <- hessABCp <- hessABIpw <- hessABFull <- matrix(0,nrow=length(tauCut),ncol=ncol(valX))
  #   hessBB <- hessBBCp <- hessBBIpw <- hessBBFull <- matrix(0,nrow=ncol(valX),ncol=ncol(valX))
  #   
  #   
  #   for (i in 1:nrow(valX)) {
  #     
  #     # qwer1 <-c(valY[i],valZ[i,])
  #     # ind1 <- which(apply(matrix(c(as.numeric(nonvalY),as.numeric(nonvalZ)),ncol=length(qwer1)),1,'identical',x=qwer1)==TRUE)
  #     
  #     ind1Y <- which(nonvalY==valY[i])
  #     tmp1Z <- (nonvalZ - matrix(rep(valZ[i,],nrow(nonvalZ)),nrow=nrow(nonvalZ),ncol=ncol(nonvalZ),byrow=TRUE))^2
  #     tmp1Z <- c(tmp1Z%*%rep(1,ncol(nonvalZ)))
  #     ind1Z <- which(tmp1Z==0)
  #     ind1 <- intersect(ind1Y,ind1Z)
  #     
  #     # qwer2 <-c(valY[i],valZ[i,])
  #     # ind2 <- which(apply(matrix(c(as.numeric(valY),as.numeric(valZ)),ncol=length(qwer2)),1,'identical',x=qwer2)==TRUE)
  #     
  #     ind2Y <- which(valY==valY[i])
  #     tmp2Z <- (valZ - matrix(rep(valZ[i,],nrow(valZ)),nrow=nrow(valZ),ncol=ncol(valZ),byrow=TRUE))^2
  #     tmp2Z <- c(tmp2Z%*%rep(1,ncol(valZ)))
  #     ind2Z <- which(tmp2Z==0)
  #     ind2 <- intersect(ind2Y,ind2Z)
  #     
  #     for (j in 1:min(valY[i],max(tauCut))) {
  #       
  #       ind1D <- which(nonvalD[,j]==valD[i,j])
  #       # ind1D <- which(apply(nonvalD[,1:j],1,'sum')==sum(valD[i,1:j]))
  #       ind1 <- intersect(ind1,ind1D)
  #       # ind1 <- which(apply(matrix(c(as.numeric(nonvalY),as.numeric(nonvalD[,j]),as.numeric(nonvalZ)),ncol=length(qwer1)),1,'identical',x=qwer1)==TRUE)
  #       I1 <- length(ind1)
  #       
  #       ind2D <- which(valD[,j]==valD[i,j])
  #       # ind2D <- which(apply(valD[,1:j],1,'sum')==sum(valD[i,1:j]))
  #       ind2 <- intersect(ind2,ind2D)
  #       # ind2 <- which(apply(matrix(c(as.numeric(valY),as.numeric(valD[,j]),as.numeric(valZ)),ncol=length(qwer1)),1,'identical',x=qwer1)==TRUE)
  #       I2 <- length(ind2)
  #       
  #       mu1 <- invLogit(tmp1Alpha[j]+c(t(tmp1Beta)%*%valX[i,]))
  #       mu2 <- invLogit(tmp2Alpha[j]+c(t(tmp2Beta)%*%valX[i,]))
  #       mu3 <- invLogit(tmp3Alpha[j]+c(t(tmp3Beta)%*%valX[i,]))
  #       mu4 <- invLogit(tmp4Alpha[j]+c(t(tmp4Beta)%*%valX[i,]))
  #       
  #       if (length(ind2)==0) {
  #         I2 <- 1
  #         # I1 <- 0 # correction on Sep 20, 2018
  #         I1 <- -1 # correction on Jun 6, 2019
  #         # print('warning: sparse discretization')
  #       }
  #       
  #       # print(I1/I2)
  #       
  #       # mean score (phase 2 + auxiliary phase 1)
  #       scoreA[j] <- scoreA[j] + (1+I1/I2)*(valD[i,j] - mu1)
  #       hessAA[j,j] <- hessAA[j,j] - (1+I1/I2)*mu1*(1-mu1)
  #       hessAB[j,] <- hessAB[j,] - (1+I1/I2)*mu1*(1-mu1)*valX[i,]
  #       
  #       scoreB <- scoreB + (1+I1/I2)*(valD[i,j] - mu1)*valX[i,]
  #       hessBB <- hessBB - (1+I1/I2)*mu1*(1-mu1)*valX[i,]%*%t(valX[i,])
  #       
  #       # complete (phase 2 only)
  #       scoreACp[j] <- scoreACp[j] + (valD[i,j] - mu2)
  #       hessAACp[j,j] <- hessAACp[j,j] - mu2*(1-mu2)
  #       hessABCp[j,] <- hessABCp[j,] - mu2*(1-mu2)*valX[i,]
  #       
  #       scoreBCp <- scoreBCp + (valD[i,j] - mu2)*valX[i,]
  #       hessBBCp <- hessBBCp - mu2*(1-mu2)*valX[i,]%*%t(valX[i,])
  #       
  #       # full cohorts
  #       scoreAFull[j] <- scoreAFull[j] + (valD[i,j] - mu3)
  #       hessAAFull[j,j] <- hessAAFull[j,j] - mu3*(1-mu3)
  #       hessABFull[j,] <- hessABFull[j,] - mu3*(1-mu3)*valX[i,]
  #       
  #       scoreBFull <- scoreBFull + (valD[i,j] - mu3)*valX[i,]
  #       hessBBFull <- hessBBFull - mu3*(1-mu3)*valX[i,]%*%t(valX[i,])
  #       
  #       # IPW
  #       scoreAIpw[j] <- scoreAIpw[j] + (valD[i,j] - mu4)*valWt[i]
  #       hessAAIpw[j,j] <- hessAAIpw[j,j] - mu4*(1-mu4)*valWt[i]
  #       hessABIpw[j,] <- hessABIpw[j,] - mu4*(1-mu4)*valX[i,]*valWt[i]
  #       
  #       scoreBIpw <- scoreBIpw + (valD[i,j] - mu4)*valX[i,]*valWt[i]
  #       hessBBIpw <- hessBBIpw - mu4*(1-mu4)*valX[i,]%*%t(valX[i,])*valWt[i]
  #       
  #       
  #     }
  #   }
  #   
  #   # extra updates for full cohorts
  #   for (i in 1:nrow(nonvalX)) {
  #     
  #     for (j in 1:min(nonvalY[i],max(tauCut))) {
  #       
  #       mu3 <- invLogit(tmp3Alpha[j]+c(t(tmp3Beta)%*%nonvalX[i,]))
  #       
  #       # full cohorts
  #       scoreAFull[j] <- scoreAFull[j] + (nonvalD[i,j] - mu3)
  #       hessAAFull[j,j] <- hessAAFull[j,j] - mu3*(1-mu3)
  #       hessABFull[j,] <- hessABFull[j,] - mu3*(1-mu3)*nonvalX[i,]
  #       
  #       scoreBFull <- scoreBFull + (nonvalD[i,j] - mu3)*nonvalX[i,]
  #       hessBBFull <- hessBBFull - mu3*(1-mu3)*nonvalX[i,]%*%t(nonvalX[i,])
  #       
  #     }
  #   }
  #   
  #   # mean score
  #   score <- c(scoreA,scoreB)
  #   
  #   hess1 <- cbind(hessAA,hessAB)
  #   hess2 <- cbind(t(hessAB),hessBB)
  #   hess <- rbind(hess1,hess2)
  #   
  #   # complete
  #   scoreCp <- c(scoreACp,scoreBCp)
  #   
  #   hess1Cp <- cbind(hessAACp,hessABCp)
  #   hess2Cp <- cbind(t(hessABCp),hessBBCp)
  #   hessCp <- rbind(hess1Cp,hess2Cp)
  #   
  #   # full cohorts
  #   scoreFull <- c(scoreAFull,scoreBFull)
  #   
  #   hess1Full <- cbind(hessAAFull,hessABFull)
  #   hess2Full <- cbind(t(hessABFull),hessBBFull)
  #   hessFull <- rbind(hess1Full,hess2Full)
  #   
  #   # IPW
  #   scoreIpw <- c(scoreAIpw,scoreBIpw)
  #   
  #   hess1Ipw <- cbind(hessAAIpw,hessABIpw)
  #   hess2Ipw <- cbind(t(hessABIpw),hessBBIpw)
  #   hessIpw <- rbind(hess1Ipw,hess2Ipw)
  #   
  #   
  #   # Newton-Raphson
  #   coefMs <- coefMs - ginv(hess)%*%score
  #   coefCp <- coefCp - ginv(hessCp)%*%scoreCp
  #   coefFull <- coefFull - ginv(hessFull)%*%scoreFull
  #   coefIpw <- coefIpw - ginv(hessIpw)%*%scoreIpw
  #   
  #   epsMs <- max(abs(coefMs-c(tmp1Alpha,tmp1Beta)))
  #   epsCp <- max(abs(coefCp-c(tmp2Alpha,tmp2Beta)))
  #   epsFull <- max(abs(coefFull-c(tmp3Alpha,tmp3Beta)))
  #   epsIpw <- max(abs(coefIpw-c(tmp4Alpha,tmp4Beta)))
  #   
  #   
  #   epsAll <- c(epsMs,epsCp,epsIpw,epsFull)
  #   eps <- max(epsAll)
  #   
  #   iter <- iter+1
  #   # print(paste('update eps: ',round(epsAll,3),sep=''))
  #   
  #   if (eps-eps0 > 0) {
  #     print('The algorithm may not converge')
  #     break()
  #   }
  #   
  # }
  # 
  # 
  # return(list(alphaMs=tmp1Alpha, betaMs=tmp1Beta,
  #             alphaCp=tmp2Alpha, betaCp=tmp2Beta,
  #             alphaFull=tmp3Alpha, betaFull=tmp3Beta,
  #             alphaIpw=tmp4Alpha, betaIpw=tmp4Beta,
  #             infoFull=ginv(-hessFull/length(obsY)),
  #             infoMs=ginv(-hess/length(obsY)),
  #             infoIpw=ginv(-hessIpw/length(obsY)),
  #             infoCp=ginv(-hessCp/length(valY))))
  # 
  # return(list(alphaMs=tmp1Alpha, betaMs=tmp1Beta,
  #             alphaCp=tmp2Alpha, betaCp=tmp2Beta,
  #             alphaFull=tmp3Alpha, betaFull=tmp3Beta,
  #             alphaIpw=tmp4Alpha, betaIpw=tmp4Beta,
  #             infoFull=-hessFull/length(obsY),
  #             infoMs  =-hess/length(obsY),
  #             infoIpw =-hessIpw/length(obsY),
  #             infoCp  =-hessCp/length(valY)))
  
  
}







asympConf <- function(obsY,obsD,obsZ,tauCut,indPhase,alphaTmp,betaTmp,wt=NULL,link) {
  
  N <- length(obsY)
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  if (is.null(wt)) {
    
    wt <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
    wt <- rep(wt,length(obsY))
    
    # wt[indPhase1] <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
    # wt[indPhase2] <- 1/(length(indPhase2)/(length(indPhase1) + length(indPhase2)))
  }
  
  valWt <- wt[indPhase2]
  nonvalWt <- wt[indPhase1]
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  phase1Mat <- obsMat[indPhase1,]
  phase2Mat <- obsMat[indPhase2,]
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- indStrat1 <- indStrat2 <- list()
  tableStrat <- c()
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        # tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tmp1 <- which(apply(phase1Mat,1,'identical',y=c(i,j,k))==TRUE)
        tmp2 <- which(apply(phase2Mat,1,'identical',y=c(i,j,k))==TRUE)
        
        tableStrat <- rbind(tableStrat,c(i,j,k,length(tmp1),length(tmp2)))
        
        indStrat1[[l]] <- indPhase1[tmp1]
        indStrat2[[l]] <- indPhase2[tmp2]
        # indStrat[[l]] <- tmp
        l <- l+1
      }
    }
  }
  nonZeroStrat <- which((tableStrat[,4]+tableStrat[,5])!=0)
  
  # indStrat <- indStrat[nonZeroStrat]  
  indStrat1 <- indStrat1[nonZeroStrat]  
  indStrat2 <- indStrat2[nonZeroStrat] 
  
  tableStrat <- tableStrat[nonZeroStrat,] 
  
  
  
  valX <- as.matrix(obsX[indPhase2,])
  nonvalX <- as.matrix(obsX[indPhase1,])
  
  valZ <- as.matrix(obsZ[indPhase2,])
  nonvalZ <- as.matrix(obsZ[indPhase1,])
  
  valY <- ceiling(obsY[indPhase2])
  nonvalY <- ceiling(obsY[indPhase1])
  
  # pseudo outervations
  D <- matrix(0,nrow=length(obsY),ncol=length(tauCut))
  for (i in 1:length(obsY)) {
    for (j in 1:length(tauCut)) {
      if ((tmpY[i]) > (tauCut[j]-1) & (tmpY[i]) <= tauCut[j] & tmpD[i]==1) {
        D[i,j] <- 1
      }
    }
  }
  colnames(D) <- tauCut
  
  # cbind(valY,D[indPhase2,],obsY[indPhase2])[1:50,]
  
  valY[which(valY>max(tauCut))] <- max(tauCut)
  nonvalY[which(nonvalY>max(tauCut))] <- max(tauCut)
  
  valD <- D[indPhase2,]
  nonvalD <- D[indPhase1,]
  
  # cbind(valY,valD,obsY[indPhase2],obsD[indPhase2])[1:50,]
  
  combY <- c(valY,nonvalY)
  combD <- rbind(valD,nonvalD)
  combZ <- rbind(valZ,nonvalZ)
  combX <- rbind(valX,nonvalX)
  combWt <- c(valWt,nonvalWt)
  
  combY[indPhase1] <- nonvalY
  combY[indPhase2] <- valY
  
  combD[indPhase1,] <- nonvalD
  combD[indPhase2,] <- valD
  
  combZ[indPhase1,] <- nonvalZ
  combZ[indPhase2,] <- valZ
  
  combX[indPhase1,] <- nonvalX
  combX[indPhase2,] <- valX
  
  combWt[indPhase1] <- nonvalWt
  combWt[indPhase2] <- valWt
  
  
  
  omega <- 0
  
  hessAA <- matrix(0,nrow=length(tauCut),ncol=length(tauCut))
  hessAB <- matrix(0,nrow=length(tauCut),ncol=ncol(combX))
  hessBB <- matrix(0,nrow=ncol(combX),ncol=ncol(combX))
  
  for (l in 1:nrow(tableStrat)) {
    
    # Nyz <- tableStrat[l,6]
    n1yz <- tableStrat[l,4]
    n2yz <- tableStrat[l,5]
    Nyz <- n1yz + n2yz
    
    if (n2yz!=0) {
      
      scoreStrat <- c()
      
      for (i in indStrat2[[l]]) {
        
        scoreA <-rep(0,length(tauCut))
        scoreB <-rep(0,ncol(combX))
        
        for (j in 1:min(combY[i],max(tauCut))) {
          
          if (link=='logit') {
            
            mu <- invLogit(alphaTmp[j]+c(t(betaTmp)%*%combX[i,]))
            
            scoreA[j] <- (combD[i,j] - mu)#*combWt[i]
            scoreB <- scoreB + (combD[i,j] - mu)*combX[i,]#*combWt[i]
            
            hessAA[j,j] <- hessAA[j,j] - mu*(1-mu)*(1+n1yz/n2yz)#*combWt[i]
            hessAB[j,] <- hessAB[j,] - mu*(1-mu)*combX[i,]*(1+n1yz/n2yz)#*combWt[i]
            hessBB <- hessBB - mu*(1-mu)*combX[i,]%*%t(combX[i,])*(1+n1yz/n2yz)#*combWt[i]
            
          } else if (link=='clogclog') {
            
            cll <- alphaTmp[j]+c(t(betaTmp)%*%combX[i,])
            mu <- invCll(cll)
            
            tmp1 <- combD[i,j]*exp(cll)/mu - exp(cll)
            tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll))/mu) - exp(cll)
            
            scoreA[j] <- tmp1#*combWt[i]
            scoreB <- scoreB + tmp1*combX[i,]#*combWt[i]
            
            hessAA[j,j] <- hessAA[j,j] + tmp2*(1+n1yz/n2yz)#*combWt[i]
            hessAB[j,] <- hessAB[j,] + tmp2*combX[i,]*(1+n1yz/n2yz)#*combWt[i]
            hessBB <- hessBB + tmp2*combX[i,]%*%t(combX[i,])*(1+n1yz/n2yz)#*combWt[i]
            
          }
        }
        
        score <- c(scoreA,scoreB)
        scoreStrat <- rbind(scoreStrat,score)
      }
      
      if (length(indStrat2[[l]])<2) {
        V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
      } else {
        V <- var(scoreStrat)
      }
      
      omega <- omega + Nyz/N*n1yz/n2yz*V
      
    }
    
  }
  
  hess1 <- cbind(hessAA,hessAB)
  hess2 <- cbind(t(hessAB),hessBB)
  hess <- rbind(hess1,hess2)
  
  invA <- ginv(-hess/N)
  B <- omega
  
  return(list(invA=invA,B=B))
  
}





optDist <- function(n2yz=NULL,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link) {
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  if (is.null(n2yz)) {
    
    # optimally distributed sampling design
    numer <- list()
    denom <- 0
    n2Ayz <- rep(0,nrow(pilotTable))
    for (l in 1:nrow(pilotTable)) {
      
      if (pilotTable[l,5]==0) {
        
        numer[[l]] <- rep(0,length(alphaTmp)+length(betaTmp))
        
      } else {
        
        Nyz <- pilotTable[l,5]
        Pyz <- Nyz/N
        
        n2Ayz[l] <- pilotTable[l,4]
        ind2Ayz <- intersect(indStrat[[l]],indPhase2)
        
        scoreStrat <- c()
        
        for (i in ind2Ayz) {
          
          scoreA <-rep(0,length(tauCut))
          scoreB <-rep(0,ncol(combX))
          
          for (j in 1:min(combY[i],max(tauCut))) {
            
            if (link=='logit') {
              
              mu <- invLogit(alphaTmp[j]+c(t(betaTmp)%*%combX[i,]))
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              scoreA[j] <- combD[i,j] - mu
              scoreB <- scoreB + (combD[i,j] - mu)*combX[i,]
              
            } else if (link=='cloglog') {
              
              cll <- alphaTmp[j]+c(t(betaTmp)%*%combX[i,])
              mu <- invCll(cll)
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              tmp1 <- combD[i,j]*exp(cll)/mu - exp(cll)
              # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
              
              scoreA[j] <- tmp1
              scoreB <- scoreB + tmp1*combX[i,]
              
            }
            
          }
          
          score <- c(scoreA,scoreB)
          scoreStrat <- rbind(scoreStrat,score)
          
        }
        
        if (length(ind2Ayz)<2) {
          V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
        } else {
          V <- var(scoreStrat)
        }
        
        # round(V,3)
        # round(A,3)
        
        # tmp0 <- ginv(A)%*%V%*%ginv(A)
        tmp0 <- ginv(infoTmp)%*%V%*%ginv(infoTmp)
        tmp0Svd <- svd(tmp0)
        
        svdU <- tmp0Svd$u
        svdD <- tmp0Svd$d
        svdV <- tmp0Svd$v
        
        tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
        
        numer[[l]] <- (1-r)*N*sqrt(Nyz)*sqrt(Pyz)*diag(tmp)
        denom <- denom + sqrt(Nyz)*sqrt(Pyz)*diag(tmp)
        
      }
    }
    
    n2Byz <- matrix(0,nrow=nrow(pilotTable),ncol=(length(alphaTmp)+length(betaTmp)))
    for (l in 1:nrow(pilotTable)) {
      if (pilotTable[l,5]!=0) {
        n2Byz[l,] <- numer[[l]]/denom - n2Ayz[l]
      }
    }
    apply(n2Byz,2,'sum')
    
  } else {
    
    # re-distributed design
    n2Ayz <- n2yz$A
    n2Byz0 <- n2Byz <- n2yz$B
    for (lc in 1:ncol(n2Byz)) {
      
      numer <- rep(NA,nrow(pilotTable))
      denom <- 0
      # n2Ayz <- rep(0,nrow(pilotTable))
      for (lr in 1:nrow(pilotTable)) {
        
        if (n2Byz0[lr,lc] <= 0) {
          
          # print(c(1,lr))
          
          n2Byz[lr,lc] <- 0 
          next()
          
        } else if (n2Byz0[lr,lc] > (pilotTable[lr,5] - n2Ayz[lr])) {
          
          # print(c(2,lr))
          
          n2Byz[lr,lc] <- pilotTable[lr,5] - n2Ayz[lr] 
          next()
          
        } else {
          
          # print(c(3,lr))
          
          Nyz <- pilotTable[lr,5]
          Pyz <- Nyz/N
          
          # n2Ayz[lr] <- pilotTable[lr,4]
          ind2Ayz <- intersect(indStrat[[lr]],indPhase2)
          
          scoreStrat <- c()
          
          for (i in ind2Ayz) {
            
            scoreA <-rep(0,length(tauCut))
            scoreB <-rep(0,ncol(combX))
            
            for (j in 1:min(combY[i],max(tauCut))) {
              
              if (link=='logit') {
                
                mu <- invLogit(alphaTmp[j]+c(t(betaTmp)%*%combX[i,]))
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                scoreA[j] <- combD[i,j] - mu
                scoreB <- scoreB + (combD[i,j] - mu)*combX[i,]
                
              } else if (link=='cloglog') {
                
                cll <- alphaTmp[j]+c(t(betaTmp)%*%combX[i,])
                mu <- invCll(cll)
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                tmp1 <- combD[i,j]*exp(cll)/mu - exp(cll)
                # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
                
                scoreA[j] <- tmp1
                scoreB <- scoreB + tmp1*combX[i,]
                
              }
              
            }
            
            score <- c(scoreA,scoreB)
            scoreStrat <- rbind(scoreStrat,score)
            
          }
          
          if (length(ind2Ayz)<2) {
            V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
          } else {
            V <- var(scoreStrat)
          }
          
          tmp0 <- ginv(infoTmp)%*%V%*%ginv(infoTmp)
          tmp0Svd <- svd(tmp0)
          
          svdU <- tmp0Svd$u
          svdD <- tmp0Svd$d
          svdV <- tmp0Svd$v
          
          tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
          
          ind1 <- which(n2Byz0[,lc] <= 0)
          ind2 <- which(n2Byz0[,lc] > (pilotTable[,5] - n2Ayz))
          
          nAdj <- (1-r)*N - sum(n2Ayz[ind1]) - sum(pilotTable[ind2,5])#-n2Ayz[ind2])
          
          # print(c(nAdj,sum(n2Ayz[ind1]),sum(pilotTable[ind2,5]-n2Ayz[ind2])))
          
          numer[lr] <- nAdj*sqrt(Nyz)*sqrt(Pyz)*diag(tmp)[lc]
          denom <- denom + sqrt(Nyz)*sqrt(Pyz)*diag(tmp)[lc]
          
        }
        
      }
      # sum(numer/denom)
      # sum(n2Ayz[which(n2Byz0[,lc]<=0)])
      
      
      # ind1 <- which(n2Byz0[,lc] >= 0)
      # ind2 <- which(n2Byz0[,lc] <= (pilotTable[,5] - n2Ayz))
      # ind3 <- which(pilotTable[,5]!=0)
      # 
      # ind <- unique(sort(c(ind1,ind2,ind3)))
      
      for (lr in which(is.na(numer)==FALSE)) {
        n2Byz[lr,lc] <- numer[lr]/denom - n2Ayz[lr]
      }
      # n2Byz[,lc]
      # sum(n2Byz[,lc])
      # cbind(n2Ayz,n2Byz[,lc])
      # sum(n2Ayz)
      # sum(n2Byz[,lc])
      
      asdf <- cbind(numer/denom,n2Ayz,n2Byz0[,lc],n2Byz[,lc],pilotTable[,5])
      colnames(asdf) <- c('design','n2Ayz','n2Byz0','n2Byz','total')
      asdf
      sum(n2Byz[,lc])
      
    }
    apply(n2Byz,2,'sum') + sum(n2Ayz)
    # which(n2Byz<0,arr.ind=TRUE)
    
  }
  
  return(list(n2Ayz=n2Ayz,n2Byz=n2Byz))
  
}





msDesignNew <- function(obsY,obsD,obsZ,obsX,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,pilotTable,indStrat,link){ 
  
  N <- length(obsY)
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  valX <- as.matrix(obsX[indPhase2,])
  nonvalX <- as.matrix(obsX[indPhase1,])
  
  valZ <- as.matrix(obsZ[indPhase2,])
  nonvalZ <- as.matrix(obsZ[indPhase1,])
  
  valY <- ceiling(obsY[indPhase2])
  nonvalY <- ceiling(obsY[indPhase1])
  
  # pseudo outervations
  D <- matrix(0,nrow=length(obsY),ncol=length(tauCut))
  for (i in 1:length(obsY)) {
    for (j in 1:length(tauCut)) {
      if ((tmpY[i]) > (tauCut[j]-1) & (tmpY[i]) <= tauCut[j] & tmpD[i]==1) {
        D[i,j] <- 1
      }
    }
  }
  colnames(D) <- tauCut
  
  # cbind(valY,D[indPhase2,],obsY[indPhase2])[1:50,]
  
  valY[which(valY>max(tauCut))] <- max(tauCut)
  nonvalY[which(nonvalY>max(tauCut))] <- max(tauCut)
  
  valD <- D[indPhase2,]
  nonvalD <- D[indPhase1,]
  
  # cbind(valY,valD,outY[indPhase2],outD[indPhase2])[1:50,]
  
  combY <- c(valY,nonvalY)
  combD <- rbind(valD,nonvalD)
  combZ <- rbind(valZ,nonvalZ)
  combX <- rbind(valX,nonvalX)
  
  combY[indPhase1] <- nonvalY
  combY[indPhase2] <- valY
  
  combD[indPhase1,] <- nonvalD
  combD[indPhase2,] <- valD
  
  combZ[indPhase1,] <- nonvalZ
  combZ[indPhase2,] <- valZ
  
  combX[indPhase1,] <- nonvalX
  combX[indPhase2,] <- valX
  
  
  # optimally distributed sampling design
  optTmp <- optDist(n2yz=NULL,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link=link)
  n2Ayz <- optTmp$n2Ayz
  n2Byz0 <- n2ByzOpt <- optTmp$n2Byz
  
  round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  apply(n2Byz0,2,'sum');sum(n2Ayz)
  
  # round(cbind(pilotTable[,c(1:3,5)],n2Ayz,n2Byz0),1)
  
  # ind <- which(n2Byz0<0,arr.ind=TRUE)
  # n2Byz0[ind] <- 0
  # 
  # round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  # apply(n2Byz0,2,'sum') + sum(n2Ayz)
  
  
  # re-distributed design 
  iter <- 1
  while (nrow(which(n2Byz0<0,arr.ind=TRUE))>0) {
    # print(iter)
    iter <- iter+1
    n2yz <- list(A=n2Ayz,B=n2Byz0)
    optTmp <- optDist(n2yz,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link=link)
    # n2Ayz <- optTmp$n2Ayz
    n2Byz0 <- optTmp$n2Byz
    apply(round(n2Byz0,3),2,'sum');sum(n2Ayz)
    
  }
  
  apply(round(n2Byz0,1),2,'sum') + sum(n2Ayz)
  which(n2Byz0<0,arr.ind=TRUE)
  
  # round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  
  
  n2Byz <- floor(n2Byz0)
  tmp <- n2Byz0 - n2Byz
  
  for (lc in 1:ncol(n2Byz)) {
    nRemain <- (1-r)*N - sum(n2Byz[,lc]) - sum(n2Ayz)
    if (nRemain > 0) {
      indAdd <- order(tmp[,lc],decreasing=TRUE)[1:nRemain]
      tmp[indAdd,lc] <- 1
      tmp[which(tmp[,lc]<1),lc] <- 0
    }
  }
  n2Byz <- n2Byz + tmp
  
  apply(round(n2Byz,1),2,'sum') + sum(n2Ayz)
  which(n2Byz<0,arr.ind=TRUE)
  
  round(cbind(n2Byz,n2Ayz,pilotTable[,5]),3)
  
  
  # design - McIsaac and Cook
  indSample <- prob <- matrix(0,nrow=nrow(nonvalZ),ncol=ncol(infoTmp))
  for (lc in 1:ncol(n2Byz)) {
    
    for (lr in 1:nrow(n2Byz)) {
      
      if (n2Byz[lr,lc]!=0) {
        
        ind2Byz <- intersect(indStrat[[lr]],indPhase1)
        indTmp <- sample(ind2Byz,min(n2Byz[lr,lc],length(ind2Byz)))
        
        indSample[match(indTmp,indPhase1),lc] <- 1
        
        prob[match(ind2Byz,indPhase1),lc] <- n2Byz0[lr,lc]/length(ind2Byz)
        
      }
    }
    
  }
  
  return(list(indSample=indSample,
              prob=prob,
              n2Ayz=n2Ayz,
              n2ByzOpt=n2ByzOpt,
              n2Byz0=n2Byz0,
              n2Byz=n2Byz))
  
}







preDist <- function(n2yz=NULL,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link) {
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  if (is.null(n2yz)) {
    
    # preimally distributed sampling design
    numer <- list()
    denom <- 0
    n2Ayz <- rep(0,nrow(pilotTable))
    for (l in 1:nrow(pilotTable)) {
      
      if (pilotTable[l,5]==0) {
        
        numer[[l]] <- rep(0,length(alphaTmp)+length(betaTmp))
        
      } else {
        
        Nyz <- pilotTable[l,5]
        Pyz <- Nyz/N
        
        n2Ayz[l] <- pilotTable[l,4]
        ind2Ayz <- intersect(indStrat[[l]],indPhase2)
        
        scoreStrat <- c()
        
        for (i in ind2Ayz) {
          
          scoreA <-rep(0,length(tauCut))
          scoreB <-rep(0,ncol(combX))
          
          for (j in 1:min(combY[i],max(tauCut))) {
            
            if (link=='logit') {
              
              mu <- invLogit(alphaTmp[j]+c(t(betaTmp)%*%combX[i,]))
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              scoreA[j] <- combD[i,j] - mu
              scoreB <- scoreB + (combD[i,j] - mu)*combX[i,]
              
            } else if (link=='cloglog') {
              
              cll <- alphaTmp[j]+c(t(betaTmp)%*%combX[i,])
              mu <- invCll(cll)
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              tmp1 <- combD[i,j]*exp(cll)/mu - exp(cll)
              # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
              
              scoreA[j] <- tmp1
              scoreB <- scoreB + tmp1*combX[i,]
              
            }
            
          }
          
          score <- c(scoreA,scoreB)
          scoreStrat <- rbind(scoreStrat,score)
          
        }
        
        if (length(ind2Ayz)<2) {
          V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
        } else {
          V <- var(scoreStrat)
        }
        
        # round(V,3)
        # round(A,3)
        
        # tmp0 <- ginv(A)%*%V%*%ginv(A)
        tmp0 <- ginv(infoTmp)%*%V%*%ginv(infoTmp)
        tmp0Svd <- svd(tmp0)
        
        svdU <- tmp0Svd$u
        svdD <- tmp0Svd$d
        svdV <- tmp0Svd$v
        
        tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
        
        numer[[l]] <- (r0-r)*N*Pyz*diag(tmp)
        denom <- denom + Pyz*diag(tmp)
        
      }
    }
    
    n2Byz <- matrix(0,nrow=nrow(pilotTable),ncol=(length(alphaTmp)+length(betaTmp)))
    for (l in 1:nrow(pilotTable)) {
      if (pilotTable[l,5]!=0) {
        n2Byz[l,] <- numer[[l]]/denom# - n2Ayz[l]
      }
    }
    apply(n2Byz,2,'sum')
    
  } else {
    
    # re-distributed design
    n2Ayz <- n2yz$A
    n2Byz0 <- n2Byz <- n2yz$B
    for (lc in 1:ncol(n2Byz)) {
      
      numer <- rep(NA,nrow(pilotTable))
      denom <- 0
      # n2Ayz <- rep(0,nrow(pilotTable))
      for (lr in 1:nrow(pilotTable)) {
        
        if (n2Byz0[lr,lc] <= 0) {
          
          # print(c(1,lr))
          
          n2Byz[lr,lc] <- 0 
          next()
          
        } else if (n2Byz0[lr,lc] > (pilotTable[lr,5] - n2Ayz[lr])) {
          
          # print(c(2,lr))
          
          n2Byz[lr,lc] <- pilotTable[lr,5] - n2Ayz[lr] 
          next()
          
        } else {
          
          # print(c(3,lr))
          
          Nyz <- pilotTable[lr,5]
          Pyz <- Nyz/N
          
          # n2Ayz[lr] <- pilotTable[lr,4]
          ind2Ayz <- intersect(indStrat[[lr]],indPhase2)
          
          scoreStrat <- c()
          
          for (i in ind2Ayz) {
            
            scoreA <-rep(0,length(tauCut))
            scoreB <-rep(0,ncol(combX))
            
            for (j in 1:min(combY[i],max(tauCut))) {
              
              if (link=='logit') {
                
                mu <- invLogit(alphaTmp[j]+c(t(betaTmp)%*%combX[i,]))
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                scoreA[j] <- combD[i,j] - mu
                scoreB <- scoreB + (combD[i,j] - mu)*combX[i,]
                
              } else if (link=='cloglog') {
                
                cll <- alphaTmp[j]+c(t(betaTmp)%*%combX[i,])
                mu <- invCll(cll)
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                tmp1 <- combD[i,j]*exp(cll)/mu - exp(cll)
                # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
                
                scoreA[j] <- tmp1
                scoreB <- scoreB + tmp1*combX[i,]
                
              }
              
            }
            
            score <- c(scoreA,scoreB)
            scoreStrat <- rbind(scoreStrat,score)
            
          }
          
          if (length(ind2Ayz)<2) {
            V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
          } else {
            V <- var(scoreStrat)
          }
          
          tmp0 <- ginv(infoTmp)%*%V%*%ginv(infoTmp)
          tmp0Svd <- svd(tmp0)
          
          svdU <- tmp0Svd$u
          svdD <- tmp0Svd$d
          svdV <- tmp0Svd$v
          
          tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
          
          ind1 <- which(n2Byz0[,lc] <= 0)
          ind2 <- which(n2Byz0[,lc] > (pilotTable[,5] - n2Ayz))
          
          nAdj <- (r0-r)*N - sum(n2Ayz[ind1]) - sum(pilotTable[ind2,5])#-n2Ayz[ind2])
          
          # print(c(nAdj,sum(n2Ayz[ind1]),sum(pilotTable[ind2,5]-n2Ayz[ind2])))
          
          numer[lr] <- nAdj*Pyz*diag(tmp)[lc]
          denom <- denom + Pyz*diag(tmp)[lc]
          
        }
        
      }
      # sum(numer/denom)
      # sum(n2Ayz[which(n2Byz0[,lc]<=0)])
      
      
      # ind1 <- which(n2Byz0[,lc] >= 0)
      # ind2 <- which(n2Byz0[,lc] <= (pilotTable[,5] - n2Ayz))
      # ind3 <- which(pilotTable[,5]!=0)
      # 
      # ind <- unique(sort(c(ind1,ind2,ind3)))
      
      for (lr in which(is.na(numer)==FALSE)) {
        n2Byz[lr,lc] <- numer[lr]/denom# - n2Ayz[lr]
      }
      # n2Byz[,lc]
      # sum(n2Byz[,lc])
      # cbind(n2Ayz,n2Byz[,lc])
      # sum(n2Ayz)
      # sum(n2Byz[,lc])
      
      asdf <- cbind(numer/denom,n2Ayz,n2Byz0[,lc],n2Byz[,lc],pilotTable[,5])
      colnames(asdf) <- c('design','n2Ayz','n2Byz0','n2Byz','total')
      asdf
      sum(n2Byz[,lc])
      
    }
    apply(n2Byz,2,'sum') + sum(n2Ayz)
    # which(n2Byz<0,arr.ind=TRUE)
    
  }
  
  return(list(n2Ayz=n2Ayz,n2Byz=n2Byz))
  
}






msDesignPre <- function(obsY,obsD,obsZ,obsX,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,pilotTable,indStrat,link){ 
  
  N <- length(obsY)
  
  indPhase1 <- indPhase[[1]]
  indPhase2 <- indPhase[[2]]
  
  valX <- as.matrix(obsX[indPhase2,])
  nonvalX <- as.matrix(obsX[indPhase1,])
  
  valZ <- as.matrix(obsZ[indPhase2,])
  nonvalZ <- as.matrix(obsZ[indPhase1,])
  
  valY <- ceiling(obsY[indPhase2])
  nonvalY <- ceiling(obsY[indPhase1])
  
  # pseudo outervations
  D <- matrix(0,nrow=length(obsY),ncol=length(tauCut))
  for (i in 1:length(obsY)) {
    for (j in 1:length(tauCut)) {
      if ((tmpY[i]) > (tauCut[j]-1) & (tmpY[i]) <= tauCut[j] & tmpD[i]==1) {
        D[i,j] <- 1
      }
    }
  }
  colnames(D) <- tauCut
  
  # cbind(valY,D[indPhase2,],obsY[indPhase2])[1:50,]
  
  valY[which(valY>max(tauCut))] <- max(tauCut)
  nonvalY[which(nonvalY>max(tauCut))] <- max(tauCut)
  
  valD <- D[indPhase2,]
  nonvalD <- D[indPhase1,]
  
  # cbind(valY,valD,outY[indPhase2],outD[indPhase2])[1:50,]
  
  combY <- c(valY,nonvalY)
  combD <- rbind(valD,nonvalD)
  combZ <- rbind(valZ,nonvalZ)
  combX <- rbind(valX,nonvalX)
  
  combY[indPhase1] <- nonvalY
  combY[indPhase2] <- valY
  
  combD[indPhase1,] <- nonvalD
  combD[indPhase2,] <- valD
  
  combZ[indPhase1,] <- nonvalZ
  combZ[indPhase2,] <- valZ
  
  combX[indPhase1,] <- nonvalX
  combX[indPhase2,] <- valX
  
  
  # optimally distributed sampling design
  preTmp <- preDist(n2yz=NULL,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link=link)
  n2Ayz <- preTmp$n2Ayz
  n2Byz0 <- n2ByzPre <- preTmp$n2Byz
  
  round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  apply(n2Byz0,2,'sum');sum(n2Ayz)
  
  # round(cbind(pilotTable[,c(1:3,5)],n2Ayz,n2Byz0),1)
  
  # ind <- which(n2Byz0<0,arr.ind=TRUE)
  # n2Byz0[ind] <- 0
  # 
  # round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  # apply(n2Byz0,2,'sum') + sum(n2Ayz)
  
  
  # re-distributed design 
  iter <- 1
  while (nrow(which(n2Byz0<0,arr.ind=TRUE))>0) {
    # print(iter)
    iter <- iter+1
    n2yz <- list(A=n2Ayz,B=n2Byz0)
    preTmp <- preDist(n2yz,pilotTable,indStrat,combY,combX,combD,indPhase,tauCut,r,alphaTmp,betaTmp,infoTmp,link=link)
    # n2Ayz <- optTmp$n2Ayz
    n2Byz0 <- preTmp$n2Byz
    apply(round(n2Byz0,3),2,'sum');sum(n2Ayz)
    
  }
  
  apply(round(n2Byz0,1),2,'sum') + sum(n2Ayz)
  which(n2Byz0<0,arr.ind=TRUE)
  
  # round(cbind(n2Byz0,n2Ayz,pilotTable[,5]),3)
  
  
  n2Byz <- floor(n2Byz0)
  tmp <- n2Byz0 - n2Byz
  
  for (lc in 1:ncol(n2Byz)) {
    nRemain <- (r0-r)*N - sum(n2Byz[,lc]) #- sum(n2Ayz)
    if (nRemain > 0) {
      indAdd <- order(tmp[,lc],decreasing=TRUE)[1:nRemain]
      tmp[indAdd,lc] <- 1
      tmp[which(tmp[,lc]<1),lc] <- 0
    }
  }
  n2Byz <- n2Byz + tmp
  
  apply(round(n2Byz,1),2,'sum') #+ sum(n2Ayz)
  which(n2Byz<0,arr.ind=TRUE)
  
  round(cbind(n2Byz,n2Ayz,pilotTable[,5]),3)
  
  
  # design - McIsaac and Cook
  indSample <- prob <- matrix(0,nrow=nrow(nonvalZ),ncol=ncol(infoTmp))
  for (lc in 1:ncol(n2Byz)) {
    
    for (lr in 1:nrow(n2Byz)) {
      
      if (n2Byz[lr,lc]!=0) {
        
        ind2Byz <- intersect(indStrat[[lr]],indPhase1)
        indTmp <- sample(ind2Byz,min(n2Byz[lr,lc],length(ind2Byz)))
        
        indSample[match(indTmp,indPhase1),lc] <- 1
        
        prob[match(ind2Byz,indPhase1),lc] <- n2Byz0[lr,lc]/length(ind2Byz)
        
      }
    }
    
  }
  
  return(list(indSample=indSample,
              prob=prob,
              n2Ayz=n2Ayz,
              n2ByzPre=n2ByzPre,
              n2Byz0=n2Byz0,
              n2Byz=n2Byz))
  
}


















genStrat <- function(obsY,obsD,obsZ,tauCut) {
  
  tmpY <- obsY
  tmpD <- obsD
  
  indOver <- which(tmpY>max(tauCut))
  tmpY[indOver] <- max(tauCut)
  tmpD[indOver] <- 0
  
  obsMat <- matrix(c(tmpY,tmpD,obsZ),nrow=length(tmpY),ncol=3)
  names(obsMat) <- NULL
  
  stratY <- as.numeric(names(table(tmpY)))
  stratD <- as.numeric(names(table(tmpD)))
  stratZ <- as.numeric(names(table(obsZ)))
  
  indStrat <- list()
  tableStrat <- c()
  wt <- rep(0,length(obsY))
  l <- 1
  for (i in stratY) {
    for (j in stratD) {
      for (k in stratZ) {
        
        tmp <- which(apply(obsMat,1,'identical',y=c(i,j,k))==TRUE)
        tableStrat <- rbind(tableStrat,c(i,j,k,length(tmp)))
        indStrat[[l]] <- tmp
        if (length(tmp)>0) {
          wt[tmp] <- 1/(length(tmp)/length(obsY))
        }
        l <- l+1
      }
    }
  }
  
  colnames(tableStrat) <- c('Y','D','Z','sizeStrat')
  
  return(list(wt=wt,
              tableStrat=tableStrat,
              indStrat=indStrat))
}




oracDist <- function(n2yz=NULL,obsTable,
                     oracTable,indStratOrac,combYOrac,combXOrac,combDOrac,tauCut,r,infoOrac,link) {
  
  newN <- sum(oracTable[,4])
  
  if (is.null(n2yz)) {
    
    # optimally distributed oracle sampling design
    numer <- list()
    denom <- 0
    for (l in 1:nrow(oracTable)) {
      
      if (oracTable[l,4]==0) {
        
        numer[[l]] <- rep(0,length(trueAlpha)+length(trueBeta))
        
      } else {
        
        Pyz <- oracTable[l,4]/newN
        
        ind2yzOrac <- indStratOrac[[l]]
        
        scoreStrat <- c()
        
        for (i in ind2yzOrac) {
          
          scoreA <-rep(0,length(tauCut))
          scoreB <-rep(0,ncol(combXOrac))
          
          for (j in 1:min(combYOrac[i],max(tauCut))) {
            
            if (link=='logit') {
              
              mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,]))
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              scoreA[j] <- combDOrac[i,j] - mu
              scoreB <- scoreB + (combDOrac[i,j] - mu)*combXOrac[i,]
              
            } else if (link=='cloglog') {
              
              cll <- trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,])
              mu <- invCll(cll)
              # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
              
              tmp1 <- combDOrac[i,j]*exp(cll)/mu - exp(cll)
              # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
              
              scoreA[j] <-tmp1
              scoreB <- scoreB + tmp1*combXOrac[i,]
              
            }
            
          }
          
          score <- c(scoreA,scoreB)
          scoreStrat <- rbind(scoreStrat,score)
          
        }
        
        if (length(ind2yzOrac)<2) {
          V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
        } else {
          V <- var(scoreStrat)
        }
        
        tmp0 <- ginv(infoOrac)%*%V%*%ginv(infoOrac)
        # tmp0 <- ginv(infoTmp)%*%VOrac%*%ginv(infoTmp)
        
        tmp0Svd <- svd(tmp0)
        
        svdU <- tmp0Svd$u
        svdD <- tmp0Svd$d
        svdV <- tmp0Svd$v
        
        tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
        
        numer[[l]] <- (1-r)*N*Pyz*diag(tmp)
        denom <- denom + Pyz*diag(tmp)
        
        # numer[[l]] <- Pyz*(1-r)*N*sqrt(diag(tmp0))
        # denom <- denom + Pyz*sqrt(diag(tmp0))
      }
    }
    
    n2yz <- matrix(0,nrow=nrow(oracTable),ncol=(length(trueAlpha)+length(trueBeta)))
    for (l in 1:nrow(oracTable)) {
      if (oracTable[l,4]!=0) {
        n2yz[l,] <- numer[[l]]/denom
      }
    }
    
    # for (lc in 1:ncol(n2yz)) {
    #   n2yz[,lc] <- n2yz[,lc]/sum(n2yz[,lc])*(1-r)*sum(obsTable[,4])
    # }
    # apply(n2yz,2,'sum')
    
  } else {
    
    # re-distributed design
    
    n2yz0 <- n2yz
    for (lc in 1:ncol(n2yz)) {
      
      numer <- rep(NA,nrow(oracTable))
      denom <- 0
      # n2Ayz <- rep(0,nrow(pilotTable))
      for (lr in 1:nrow(oracTable)) {
        
        if (n2yz0[lr,lc] <= 0) {
          
          # print(c(1,lr))
          
          n2yz[lr,lc] <- 0
          next()
          
        } else if (n2yz0[lr,lc] >= obsTable[lr,4]) {
          
          # print(c(2,lr))
          
          n2yz[lr,lc] <- obsTable[lr,4]
          next()
          
        } else {
          
          # print(c(3,lr))
          
          Pyz <- oracTable[lr,4]/newN
          
          ind2yzOrac <- indStratOrac[[lr]]
          
          scoreStrat <- c()
          
          for (i in ind2yzOrac) {
            
            scoreA <-rep(0,length(tauCut))
            scoreB <-rep(0,ncol(combXOrac))
            
            for (j in 1:min(combYOrac[i],max(tauCut))) {
              
              if (link=='logit') {
                
                mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,]))
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                scoreA[j] <- combDOrac[i,j] - mu
                scoreB <- scoreB + (combDOrac[i,j] - mu)*combXOrac[i,]
                
              } else if (link=='cloglog') {
                
                cll <- trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,])
                mu <- invCll(cll)
                # mu <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combX[i,]))
                
                tmp1 <- combDOrac[i,j]*exp(cll)/mu - exp(cll)
                # tmp2 <- combD[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll)))/mu - exp(cll)
                
                scoreA[j] <-tmp1
                scoreB <- scoreB + tmp1*combXOrac[i,]
                
              }
              
            }
            
            score <- c(scoreA,scoreB)
            scoreStrat <- rbind(scoreStrat,score)
            
          }
          
          if (length(ind2yzOrac)<2) {
            V <- matrix(0,nrow=(length(trueAlpha)+length(trueBeta)),ncol=(length(trueAlpha)+length(trueBeta)))
          } else {
            V <- var(scoreStrat)
          }
          
          tmp0 <- ginv(infoOrac)%*%V%*%ginv(infoOrac)
          tmp0Svd <- svd(tmp0)
          
          svdU <- tmp0Svd$u
          svdD <- tmp0Svd$d
          svdV <- tmp0Svd$v
          
          tmp <- svdU%*%diag(sqrt(svdD))%*%t(svdV)
          
          ind1 <- which(n2yz0[,lc] <= 0)
          ind2 <- which(n2yz0[,lc] >= obsTable[,4])
          
          nAdj <- (1-r)*N - sum(obsTable[ind2,4])
          
          # print(c(nAdj,sum(obsTable[ind2,4])))
          
          numer[lr] <- nAdj*Pyz*diag(tmp)[lc]
          denom <- denom + Pyz*diag(tmp)[lc]
          
        }
      }
      
      for (lr in which(is.na(numer)==FALSE)) {
        n2yz[lr,lc] <- numer[lr]/denom
      }
      
      apply(round(n2yz,1),2,'sum')
      which(n2yz<0,arr.ind=TRUE)
      
      tmp <- matrix(rep(obsTable[,4],ncol(n2yz)),nrow=nrow(n2yz),ncol=ncol(n2yz))
      which(n2yz>tmp,arr.ind=TRUE)
      #
      # n2yz[,lc] <- n2yz[,lc]/sum(n2yz[,lc])*(1-r)*sum(obsTable[,4])
      
    }
    which(n2yz>tmp,arr.ind=TRUE)
    
  }
  
  return(n2yz)
  
}













msDesignOracle <- function(obsY,obsD,obsZ,obsX,tauCut,r,newN=10000,seedNum=100,fitAlphaBeta=NULL,link,sig=0.1,para=c(2,3,3,3)) {
  
  set.seed(seedNum)
  
  N <- length(obsY)
  
  # generate an independent large sample
  
  if (is.null(fitAlphaBeta)) {
    
    oracSam <- dataGen(newN,1,tau,baseLambda,pCen,beta0,beta1,d0,d1,rho,link=link,sig,para)
    
    obsYOrac <- oracSam$obsY
    obsDOrac <- oracSam$obsD
    obsZOrac <- oracSam$obsZ
    obsXOrac <- oracSam$obsX
    indPhaseOrac <- oracSam$indPhase
    
    indPhase2Orac <- 1:newN
    indPhase1Orac <- (1:newN)[-indPhase2Orac]
    
  } else {
    
    newN <- N
    
    obsYOrac <- obsY
    obsDOrac <- obsD
    obsZOrac <- obsZ
    obsXOrac <- obsX
    
    indPhase2Orac <- 1:newN
    indPhase1Orac <- (1:newN)[-indPhase2Orac]
    
    indPhaseOrac <- list(indPhase1Orac,indPhase2Orac)
    
    trueAlpha <- fitAlphaBeta[[1]]
    trueBeta <- fitAlphaBeta[[2]]
    
  }
  
  tmpYOrac <- obsYOrac
  tmpDOrac <- obsDOrac
  
  indOverOrac <- which(tmpYOrac>max(tauCut))
  tmpYOrac[indOverOrac] <- max(tauCut)
  tmpDOrac[indOverOrac] <- 0
  
  valXOrac <- as.matrix(obsXOrac[indPhase2Orac,])
  nonvalXOrac <- as.matrix(obsXOrac[indPhase1Orac,])
  
  valZOrac <- as.matrix(obsZOrac[indPhase2Orac,])
  nonvalZOrac <- as.matrix(obsZOrac[indPhase1Orac,])
  
  valYOrac <- ceiling(obsYOrac[indPhase2Orac])
  nonvalYOrac <- ceiling(obsYOrac[indPhase1Orac])
  
  # pseudo outervations
  DOrac <- matrix(0,nrow=length(obsYOrac),ncol=length(tauCut))
  for (i in 1:length(obsYOrac)) {
    for (j in 1:length(tauCut)) {
      if ((tmpYOrac[i]) > (tauCut[j]-1) & (tmpYOrac[i]) <= tauCut[j] & tmpDOrac[i]==1) {
        DOrac[i,j] <- 1
      }
    }
  }
  colnames(DOrac) <- tauCut
  
  valYOrac[which(valYOrac>max(tauCut))] <- max(tauCut)
  nonvalYOrac[which(nonvalYOrac>max(tauCut))] <- max(tauCut)
  
  valDOrac <- DOrac[indPhase2Orac,]
  nonvalDOrac <- DOrac[indPhase1Orac,]
  
  combYOrac <- c(valYOrac,nonvalYOrac)
  combDOrac <- rbind(valDOrac,nonvalDOrac)
  combZOrac <- rbind(valZOrac,nonvalZOrac)
  combXOrac <- rbind(valXOrac,nonvalXOrac)
  
  combYOrac[indPhase1Orac] <- nonvalYOrac
  combYOrac[indPhase2Orac] <- valYOrac
  
  combDOrac[indPhase1Orac,] <- nonvalDOrac
  combDOrac[indPhase2Orac,] <- valDOrac
  
  combZOrac[indPhase1Orac,] <- nonvalZOrac
  combZOrac[indPhase2Orac,] <- valZOrac
  
  combXOrac[indPhase1Orac,] <- nonvalXOrac
  combXOrac[indPhase2Orac,] <- valXOrac
  
  
  
  # oracle information
  hessAAOrac <- matrix(0,nrow=length(tauCut),ncol=length(tauCut))
  hessABOrac <- matrix(0,nrow=length(tauCut),ncol=ncol(valXOrac))
  hessBBOrac <- matrix(0,nrow=ncol(valXOrac),ncol=ncol(valXOrac))
  
  for (i in 1:nrow(combXOrac)) {
    
    for (j in 1:min(combYOrac[i],max(tauCut))) {
      
      if (link=='logit') {
        
        muOrac <- invLogit(trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,]))
        
        # full cohorts
        hessAAOrac[j,j] <- hessAAOrac[j,j] - muOrac*(1-muOrac)
        hessABOrac[j,] <- hessABOrac[j,] - muOrac*(1-muOrac)*combXOrac[i,]
        hessBBOrac <- hessBBOrac - muOrac*(1-muOrac)*combXOrac[i,]%*%t(combXOrac[i,])
        
      } else if (link=='cloglog') {
        
        cll <- trueAlpha[j]+c(t(trueBeta)%*%combXOrac[i,])
        mu <- invCll(cll)
        
        # tmp1 <- combDOrac[i,j]*exp(cll)/mu - exp(cll)
        tmp2 <- combDOrac[i,j]*exp(cll)/mu*(1-exp(cll-exp(cll))/mu) - exp(cll)
        
        # full cohorts
        hessAAOrac[j,j] <- hessAAOrac[j,j] + tmp2
        hessABOrac[j,] <- hessABOrac[j,] + tmp2*combXOrac[i,]
        hessBBOrac <- hessBBOrac + tmp2*combXOrac[i,]%*%t(combXOrac[i,])
        
      }
      
    }
  }
  
  # full cohorts
  hess1Orac <- cbind(hessAAOrac,hessABOrac)
  hess2Orac <- cbind(t(hessABOrac),hessBBOrac)
  hessOrac <- rbind(hess1Orac,hess2Orac)
  
  infoOrac <- -hessOrac/newN
  
  
  
  #### stratification
  # obs table
  obsTmp <- genStrat(obsY,obsD,obsZ,tauCut)
  
  obsTable <- obsTmp$tableStrat
  indStratObs <- obsTmp$indStrat
  
  # oracle table
  oracTmp <- genStrat(obsYOrac,obsDOrac,obsZOrac,tauCut)
  
  oracTable <- oracTmp$tableStrat
  indStratOrac <- oracTmp$indStrat
  
  
  
  n2yz0 <- oracDist(n2yz=NULL,obsTable,
                    oracTable,indStratOrac,combYOrac,combXOrac,combDOrac,tauCut,r,infoOrac,link=link)
  
  round(cbind(obsTable[,c(1:4)],n2yz0),1)
  apply(n2yz0,2,'sum')
  
  # round(cbind(n2yz0[,length(trueAlpha)+1],n2Byz[,length(trueAlpha)+1]+n2Ayz,n2ByzOpt[,length(trueAlpha)+1]+n2Ayz),3)
  
  
  
  
  # re-distributed design
  iter <- 1
  tmp <- matrix(rep(obsTable[,4],ncol(n2yz0)),nrow=nrow(n2yz0),ncol=ncol(n2yz0))
  while (nrow(which(n2yz0>tmp,arr.ind=TRUE))>0) {
    # print(iter)
    iter <- iter+1
    n2yz <- n2yz0
    n2yz0 <- oracDist(n2yz,obsTable,
                      oracTable,indStratOrac,combYOrac,combXOrac,combDOrac,tauCut,r,infoOrac,link=link)
    # n2Ayz <- optTmp$n2Ayz
    apply(round(n2yz0,1),2,'sum')
    which(n2yz0>tmp,arr.ind=TRUE)
  }
  
  apply(round(n2yz0,1),2,'sum')
  which(n2yz0<0,arr.ind=TRUE)
  
  round(cbind(n2yz0,oracTable[,4]),3)
  
  
  
  n2yz <- floor(n2yz0)
  tmp <- n2yz0 - n2yz
  
  for (lc in 1:ncol(n2yz)) {
    nRemain <- (1-r)*N - sum(n2yz[,lc]) #- sum(n2Ayz)
    if (nRemain > 0) {
      indAdd <- order(tmp[,lc],decreasing=TRUE)[1:nRemain]
      tmp[indAdd,lc] <- 1
      tmp[which(tmp[,lc]<1),lc] <- 0
    }
  }
  n2yz <- n2yz + tmp
  
  round(cbind(n2yz,oracTable[,4]),3)
  apply(n2yz,2,'sum')
  
  
  indSample <- prob <- matrix(0,nrow=N,ncol=(length(trueAlpha)+length(trueBeta)))
  set.seed(seedNum)
  for (lc in 1:ncol(n2yz0)) {
    
    for (lr in 1:nrow(n2yz0)) {
      
      if (n2yz[lr,lc]!=0) {
        
        indYZ <- indStratObs[[lr]]
        indTmp <- sample(indYZ,n2yz[lr,lc])
        
        indSample[indTmp,lc] <- 1
        
        prob[indYZ,lc] <- length(indTmp)/length(indYZ)
        
      }
    }
    
  }
  
  
  # return(prob)
  return(list(indSample=indSample,
              prob=prob,
              n2yz=n2yz,
              n2yz0=n2yz0,
              indStratObs=indStratObs))
  
}


