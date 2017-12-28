sicGroupBF <- function(inData, sictest="bf", domtest="ks", plotSIC=TRUE, ...) {

  subjects <- sort(unique(inData$Subject))
  nsubjects <- length(subjects)

  conditions <- sort(unique(inData$Condition))
  nconditions <- length(conditions)

  SICnames <- c("SerialOR", "ParallelAND", "ParallelOR", "Coactive",
                "ParallelAND", "SerialAND")

  times <- sort(unique(round(inData$RT)))

  n <- 0
  nc1 <- 1
  Dom <- vector("list")

  sicAllMat <- numeric()
  subj.out <- character()
  cond.out <- character()

  if( sictest=="ks") {
    KS <- numeric()
    KSp <- numeric()
    micstat <- numeric()
    micp <- numeric()
    KSwin <- character()
    N <- numeric()
  } else if(sictest=="bf") {
    BF <- numeric()
    BFwin <- numeric() 
    BFnames <- c("SerialOR", "ParallelAND", "ParallelOR", "Coactive",
                 "ParallelAND", "SerialAND")
  } else {
    cat("Only KS-SIC test is currently implemented.\n")
    return(NA)
  }
  
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    for ( sn in 1:nsubjects ) {
      if (is.factor(subjects)) {subj <- levels(subjects)[sn]} else {subj <- subjects[sn] }
      if (plotSIC & ( sn %% 9 == 1) ) {
        dev.new()
        par(mfrow=c(3,3))
      }

      cond.out <- c(cond.out, cond)
      subj.out <- c(subj.out, subj)

      n <- n+1
      HH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==2] )
      HL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==1] )
      LH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==2] )
      LL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==1] )

      if ( min( length(HH), length(HL), length(LH), length(LL)) > 10 ) {
        sicn <- sic(HH=HH, HL=HL, LH=LH, LL=LL)
        sicAllMat <- rbind(sicAllMat, sicn$SIC(times))

        if (sictest=="ks") {
          N <- rbind(N, sicn$N)
          KS  <- rbind(KS,  sicn$Dvals[,1])
          KSp <- rbind(KSp, sicn$Dvals[,2])
          micstat <- c(micstat, sicn$MIC[[1]])
          micp <- c(micp, sicn$MIC[[2]])
          Dom[[n]] <- sicn$Dominance
          if (sicn$Dvals[1,2] < .05) {
            if (sicn$Dvals[2,2] < .05) {
              if (sicn$MIC[[2]] < .05) {
                KSwin <- c(KSwin, "Coactive")
              } else {
                KSwin <- c(KSwin, "SerialAND")
              }
            } else {
              KSwin <- c(KSwin, "ParallelOR")
            }
          } else {
            if (sicn$Dvals[2,2] < .05) {
              KSwin <- c(KSwin, "ParallelAND")
            } else {
              KSwin <- c(KSwin, "SerialOR")
            }
          }
        } else if(sictest=="bf") {
          dat <- list(HH, HL, LH, LL) 
          BF <- rbind(BF, BFsic(dat,maxn=2e4,priorinfluence=1,verbose=FALSE,
                          tolSIC=.1,tolMIC=.3,fasttest2=T)$BF)
          BFwin <- c(BFwin, BFnames[BF[n,] == max(BF[n,])])

        }
        if(plotSIC) {
          plot(times, sicn$SIC(times), type='l',
            main=paste(cond, " Condition\nParticipant ", subj, sep=""), 
            xlab="Time",ylab="SIC(t)",...)
        }
      } else { 
        #if(sictest=="ks") {BF[n,] == NA; BFwin[n] == NA }
        #if(sictest=="dp") {BF[n,] == NA; BFwin[n] == NA }
      }
    }

    if(plotSIC) {
      dev.new()
      matplot(times, t(sicAllMat[cond.out==cond,]),type='l',lty=1,
        main=paste(cond, " Condition", sep=""), 
        xlab="Time",ylab="SIC(t)",...)
    }
    nc1 <- n+1
  }

  if(sictest=="ks") {
    colnames(KS) <- c("D+", "D-")
    colnames(KSp) <- c("D+", "D-")
    statistic <- as.data.frame(list(Subject=subj.out, Condition=cond.out, N=N,
        Dpositive= KS[,1], p.val.positive=KSp[,1], Dnegative=KS[,2], p.val.negative=KSp[,2], MIC=micstat, p.val.mic=micp, Model=KSwin))
  } else if (sictest=="bf") {
    statistic <- as.data.frame(list(Subject=subj.out, Condition=cond.out, 
        BF=BF, Model=BFwin))
  }
  return(list(statistic=statistic, SIC=sicAllMat, Dominance=Dom, times=times))
}



sictestBayes <- function(HH, HL, LH, LL, method=c("DP", "IG"), model=NULL) { 
  DNAME <- paste("\nHH:", deparse(substitute(HH)), "\tHL:", deparse(substitute(HL)), 
                 "\nLH:", deparse(substitute(LH)), "\tLL:", deparse(substitute(LL)) )
  if ( method == "DP") { 
    METHOD <- "Nonparametric Bayesian SIC test"
    statistic <- sicDPtest(list(HH, HL, LH, LL))$BF
    statistic <- statistic[c(1,6,3,2,4)]
    names(statistic) <- c("Zero", "NegPos.MIC0", "Positive", "Negative", "NegPos.MICpos")

  } else if ( method == "IG") { 
    METHOD <- "Parametric Bayesian SIC test"

    dat <- list(NHH=length(HH), NHL=length(HL), NLH=length(LH), NLL=length(LL), 
                yHH=HH, yHL=HL, yLH=LH, yLL=LL)
    if (is.null(model)) { 
      coactive.model <- stan_model(file="coactive_waic.stan")
      paralleland.model <- stan_model(file="parallel-and_waic.stan")
      parallelor.model <- stan_model(file="parallel-or_waic.stan")
      serialand.model <- stan_model(file="serial-and_waic.stan")
      serialor.model <- stan_model(file="serial-or_waic.stan")
    }
    coactive.posterior <- sampling(fit=coactive.model, data=dat, chains=3, iter=80000, warmup=8000)
    coactive.posterior <- extract(coactive.posterior, c("summandsHH", "summandsHL", "summandsLH", "summandsLL"), permute=TRUE)
    coactive.posterior <- cbind(coactive.posterior$summandsHH, coactive.posterior$summandsHL, coactive.posterior$summandsLH, coactive.posterior$summandsLL)
    coactive.waic <- waic(coactive.posterior)$waic
    rm(coactive.posterior)

    paralleland.posterior <- sampling(paralleland.model, data=dat, chains=3, iter=80000, warmup=8000)
    paralleland.posterior <- extract(paralleland.posterior, c("summandsHH", "summandsHL", "summandsLH", "summandsLL"), permute=TRUE)
    paralleland.posterior <- cbind(paralleland.posterior$summandsHH, paralleland.posterior$summandsHL, paralleland.posterior$summandsLH, paralleland.posterior$summandsLL)
    paralleland.waic <- waic(paralleland.posterior)$waic
    rm(paralleland.posterior)

    parallelor.posterior <- sampling(parallelor.model, data=dat, chains=3, iter=80000, warmup=8000)
    parallelor.posterior <- extract(parallelor.posterior, c("summandsHH", "summandsHL", "summandsLH", "summandsLL"), permute=TRUE)
    parallelor.posterior <- cbind(parallelor.posterior$summandsHH, parallelor.posterior$summandsHL, parallelor.posterior$summandsLH, parallelor.posterior$summandsLL)
    parallelor.waic <- waic(parallelor.posterior)$waic
    rm(parallelor.posterior)

    serialand.posterior <- sampling(serialand.model, data=dat, chains=3, iter=80000, warmup=8000)
    serialand.posterior <- extract(serialand.posterior, c("summandsHH", "summandsHL", "summandsLH", "summandsLL"), permute=TRUE)
    serialand.posterior <- cbind(serialand.posterior$summandsHH, serialand.posterior$summandsHL, serialand.posterior$summandsLH, serialand.posterior$summandsLL)
    serialand.waic <- waic(serialand.posterior)$waic
    rm(serialand.posterior)

    serialor.posterior <- sampling(serialor.model, data=dat, chains=3, iter=80000, warmup=8000)
    serialor.posterior <- extract(serialor.posterior, c("summandsHH", "summandsHL", "summandsLH", "summandsLL"), permute=TRUE)
    serialor.posterior <- cbind(serialor.posterior$summandsHH, serialor.posterior$summandsHL, serialor.posterior$summandsLH, serialor.posterior$summandsLL)
    serialor.waic <- waic(serialor.posterior)$waic
    rm(serialor.posterior)

    statistic <- c(serialor.waic, serialand.waic, parallelor.waic, paralleland.waic, coactive.waic)
    names(statistic) <- c("SerialOR", "SerialAND", "ParallelOR", "ParallelAND", "Coactive")
  }
  rn <- list(statistic=statistic, method=METHOD, data.name=DNAME)
  return(rn)
}


# Little function to calculate posterior variances from simulation
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function (log_lik){
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}



sicDPtest <- function(dat) { 
  RTall <- unlist(dat)
   
  nbin <- min(100, min(unlist(lapply(dat,length)))-1)
  nsamp <- 1e4

  prec<-.01
  ci<-.95
  allconverge<-F
  fasttest<-F
  sig<-.95
  fasttest2<-F
  maxn<-5e5
  tolSIC<-2e-3
  tolMIC<-1e0
  verbose<-F

  # Define Bins
  pvec <- c(0, (1:(nbin-1))/nbin, 1)
  bins<- quantile(RTall, pvec)
  bins[1] <- 0 
  bins[length(bins)] <- max(RTall)+1
  names(bins) <- NULL

  # Calculate Bin Widths
  dx <- diff(bins)

  binn <- sapply(dat, tabnodrop, bins<-bins)

  # Set prior measure for Dirichlet
  alpha <- rep(1/nbin,nbin)

  # Variable to track number of samples of each type
  Nprs <- rep(0,6); 
  names(Nprs) <- c("Z","N","P","nP","Np","np")
  Npos <- Nprs; 
  NencPrs <- Nprs; 

  N <- 0; 

  repeat {

    # Sample PDF from the prior and posterior
    priorp <- samplePriorPDF(nsamp, alpha)
    postp <- samplePostPDF(nsamp, binn, alpha)

    # Convert sample PDFs to sample SICs
    priorSIC <- apply(priorp,1,pdfs2SIC)
    postSIC  <- apply(postp,1,pdfs2SIC)

    # Count number of SIC types in prior and posterior
    Npr <- apply(apply(priorSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    Npo <- apply(apply(postSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    Npos <- Npos+Npo; 
    Nprs <- Nprs+Npr

    N <- N+nsamp

    # Samples from encompassing prior
    encPrior <- array(as.vector(rdirichlet(4*nsamp,alpha)),dim=c(nsamp,4,length(alpha)))
    encPriorSIC <- apply(encPrior,1,pdfs2SIC)
    NencPr <- apply(apply(encPriorSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    NencPrs <- NencPrs + NencPr

    # Bayes factor relative to the encompassing prior
    BF <- (Npos+1)/(NencPrs+1)
    CIs <- cis(ci,Npos,NencPrs,N)

    dcis <- apply(CIs,1,diff)

    if (fasttest) {
      BFp=BF/sum(BF)
      converge <- !(any( CIs[,1]<sig & CIs[,2]>sig))
    } else if (fasttest2) {
      BFp=BF/sum(BF)
      # Check if the largest lower end CI is larger than the second largest
      #   upper CI
      converge <- max( CIs[,1]) > max( CIs[CIs[,2]!=max(CIs[,2]), 2 ] )
    } else {
      if (allconverge) { 
        converge <- all(dcis<prec)
      } else converge <- max(dcis[BF==max(BF)])<prec
    }

    nmat=rbind(Nprs,Npos); dimnames(nmat)=list(c("Prior","Posterior"),names(Nprs))
    pmat=cbind(CIs[,1],BF/sum(BF),CIs[,2]); dimnames(pmat)[[2]]=c("Lo","BFp","Hi")

    if (converge | N>=maxn) break
  }
  list(BF=BF,BFp=pmat,N=nmat,ci=ci)
}


rdirichlet <- function (n, alpha) {
   # modified from MCMCpack to never return a zero (sets zeros min double = 3e-324)
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    out=x/as.vector(sm)
    out[out<3e-324]=3e-324
    return(out)
}

tabnodrop <- function(samp,bins) {
    # Function to count the number of observations in each bin
    #  samp         Sample data to be binned
    #  bins         Bins into which sample data are divided
    hist(samp, bins, plot=F)$counts
}
  
epdf=function(nsi) {
  cumsum(nsi/sum(nsi))[-length(nsi)]
}

pdfs2SIC=function(pmat) {
# Converts a list of pdfs to an SIC
# assumes columns in  "HH" "HL" "LH" "LL" order
  pmat=apply(pmat,1,cumsum)
  out=pmat[,3]-pmat[,4]-(pmat[,1]-pmat[,2])
  #out=out[-c(1,length(out))] # get rid of end bins
  out[abs(out)<1e-9]=0       # otherwise no zeros due to numerical issues
  return(out)
}

checkmods <- function(x,dx,tolSIC=5e-2,tolMIC=1e-2) {
# checks several ordinal conditions 
# returns boolean vector with element (i)=T if model true 

  # "Z" "N" "P" "nP" "Np" "np"

  #  Check for all negative regardless of tolerance
  if (all(x <= 0) & any(x < 0) ) return(c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE))

  #  Check for all positive  regardless of tolerance
  if (all(x >= 0) & any(x > 0) ) return(c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))

  #  Check for all zero within tolerance
  if (all(abs(x)< tolSIC) ) return(c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))

  #  Indicators of times for which the SIC is nonzero
  nonzero <- (abs(x) > tolSIC)

  #  Indicators of times for which the SIC is negative
  neg = x < -tolSIC

  #  Indicators of times for which the SIC is positive
  pos = x > tolSIC

  # Check if the SIC is all negative
  if ( all(neg | !nonzero )  ) return( c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE) )

  # Check if the SIC is all positive 
  if ( all(!neg) )  return( c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE) )

  # Check if the SIC starts positive, if so it does not fit one of the classes
  if ( (x[nonzero])[1]> 0 ) return( c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE) )

  # Determine the sign changes
  d=neg[-length(neg)]!=neg[-1] 

  # Check if there is more than one sign change, if so it does not fit one of the classes
  if (sum(d)>1) return( c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE) )

  # Calculate MIC
  sumN <- -sum(x[neg]*dx[neg])
  sumP <- sum(x[!neg]*dx[!neg])

  # Check if the MIC is zero
  if(abs(sumN-sumP)<tolMIC) return( c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE) )

  # Check if MIC is positive
  if(sumN < sumP) return( c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE) )
   
  # Check if MIC is negative 
  if(sumN < sumP) return( c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE) )
   
  # Shouldn't reach here, but in case it does, it does not fit one of the classes
  return( c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE) )
}


cis=function(ci,npo,npr,N) {
# worst case credible interval ci widths for BF probs
  nbf=length(npo)
  poNci=matrix(calcci(ci,npo,N),ncol=2)*N
  prNci=matrix(calcci(ci,npr,N),ncol=2)*N
  loBF=(poNci[,1]+1)/(prNci[,2]+1)
  hiBF=(poNci[,2]+1)/(prNci[,1]+1)
  bfp=matrix(nrow=nbf,ncol=2)
  for (i in 1:nbf) {
    bfp[i,]=c(loBF[i]/(sum(hiBF[-i])+loBF[i]),hiBF[i]/(sum(loBF[-i])+hiBF[i]))
  }
  bfp
}

calcci=function(ci,n,N) {
  c(qbeta((1-ci)/2,n+1,N-n+1),qbeta(1-(1-ci)/2,n+1,N-n+1))
}

samplePriorPDF <- function(ns, alpha) {
    return(array(as.vector(rdirichlet(4*ns,alpha)),dim=c(ns,4,length(alpha))))
}

samplePostPDF <- function(ns, binn, alpha) {
    return(array(rbind(rdirichlet(ns,binn[,1]+alpha),rdirichlet(ns,binn[,2]+alpha),
                rdirichlet(ns,binn[,3]+alpha),rdirichlet(ns,binn[,4]+alpha)),
                dim=c(ns,4,length(alpha))))
}
