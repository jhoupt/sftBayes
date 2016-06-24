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


#bfSIC <- function(inData) {
#  subjvec <- sort(unique(inData$Subject))
#  condvec <- sort(unique(inData$Condition))
#  BFnames <- c("SerialOR", "ParallelAND", "ParallelOR", "Coactive",
#                "ParallelAND", "SerialAND")
#  n <- 0
#  BF <- matrix(NA, length(subjvec)*length(condvec),6)
#  BFwin <- rep(NA, length(subjvec)*length(condvec))
#  
#  dat <- new.env()
#  for (cond in condvec) {
#    for (subj in subjvec) {
#      n <- n+1
#      dat$HH <- inData$RT[inData$Subject==subj & inData$Condition==cond &
#                          inData$Correct & 
#                          inData$Channel1==2 & inData$Channel2==2]
#      dat$HL <- inData$RT[inData$Subject==subj & inData$Condition==cond &
#                          inData$Correct & 
#                          inData$Channel1==2 & inData$Channel2==1]
#      dat$LH <- inData$RT[inData$Subject==subj & inData$Condition==cond &
#                          inData$Correct & 
#                          inData$Channel1==1 & inData$Channel2==2]
#      dat$LL <- inData$RT[inData$Subject==subj & inData$Condition==cond &
#                          inData$Correct & 
#                          inData$Channel1==1 & inData$Channel2==1]
#      dat <- as.list(dat)
#      if ( min(c(lapply(dat, length), recursive=T)) > 10 ) {
#        BF[n,] <- BFsic(dat,maxn=2e4,priorinfluence=1,verbose=FALSE,
#                        tolSIC=.1,tolMIC=.3,fasttest2=T)$BF
#        BFwin[n] <- BFnames[BF[n,] == max(BF[n,])]
#      } else { BF[n,] == NA; BFwin[n] == NA }
#
#    }
#  }
#  return(list(BF=BF, BFwin=BFwin))
#}




BFsic=function(dat,                                      
               maxnbins=100,
               nbin=min(unlist(lapply(dat,length)))-1,
               bins=NULL,
               nsamp1=1e4,nsamp=1e4,
               priorinfluence=1,alpha=NULL,
               prec=.01,ci=.95,
               allconverge=F,
               fasttest=F,sig=.95,
               fasttest2=F,
               maxn=5e5,
               tolSIC=2e-3,
               tolMIC=1e0,
               verbose=F) {
  # get histograms
  if (is.null(bins)) {
    if ( !is.na(maxnbins) ) nbin=min(maxnbins,nbin)
    pvec <- c(0, (1:(nbin-1))/nbin, 1)
    allData <- unlist(dat)
    bins= quantile(allData, pvec)
    bins[1] <- 0; bins[length(bins)] <- max(allData)+1
    names(bins) <- NULL
  }

  #dx=diff(bins[-c(1,length(bins))]) # bin widths
  dx = diff(bins)
  binn <- sapply(dat, tabnodrop, bins=bins)


  # get dirichlet samples
  if (is.null(alpha)) alpha=rep(priorinfluence/nbin,nbin)
  alphaenc=rep(priorinfluence/nbin,nbin)
  Nprs=rep(0,6); 
  names(Nprs)=c("Z","N","P","nP","Np","np")
  Npos=Nprs; 
  NencPrs=Nprs; 
  N=0; 
  ns=nsamp1

  repeat {

    priorp <- samplePriorPDF(ns, alpha)
    postp <- samplePostPDF(ns, binn, alpha)

    priorSIC=apply(priorp,1,pdfs2SIC)
    postSIC=apply(postp,1,pdfs2SIC)

    Npr=apply(apply(priorSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    Npo=apply(apply(postSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    Npos=Npos+Npo; Nprs=Nprs+Npr; N=N+ns

    # get dirichlet samples from encompassing prior
    encPrior <- array(as.vector(rdirichlet(4*ns,alphaenc)),dim=c(ns,4,length(alphaenc)))
    encPriorSIC <- apply(encPrior,1,pdfs2SIC)
    NencPr <- apply(apply(encPriorSIC,2,checkmods,dx=dx,tolSIC=tolSIC,tolMIC=tolMIC),1,sum)
    NencPrs <- NencPrs + NencPr

    BF=(Npos+1)/(NencPrs+1)
    CIs=cis(ci,Npos,NencPrs,N)

    dcis=apply(CIs,1,diff)

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
    if (verbose) {
      print(N) 
      print(nmat)
      print(round(t(pmat),3))
    }
    if (converge | N>=maxn) break
    ns=nsamp
  }
  list(BF=BF,BFp=pmat,N=nmat,ci=ci)
}


rdirichlet=function (n, alpha) {
# modified from MCMCpack to never return a zero (sets zeros min double = 3e-324)
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    out=x/as.vector(sm)
    out[out<3e-324]=3e-324
    return(out)
}

tabnodrop=function(samp,bins) {
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

  if (all(x <= 0) & any(x < 0) ) return(c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE))
  if (all(x >= 0) & any(x > 0) ) return(c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))

  if (all(abs(x)< tolSIC) ) out=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE) else {
    nonzero <- (abs(x) > tolSIC)
    #x=x[nonzero] 
    #dx=dx[nonzero]
    #neg = x < 0
    neg = x < -tolSIC
    pos = x > tolSIC
    # Check if the SIC is all negative
    if ( all(neg | !nonzero ) & any(neg) ) out=c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE) else {
      # Check if the SIC is all positive 
      if ( all(!neg) ) out=c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE) else {
        # Check if the SIC starts positive
        if ( (x[nonzero])[1]> 0 ) out=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE) else {
          d=neg[-length(neg)]!=neg[-1] # Determine the sign changes
          if (sum(d)==1) {
            sumN=-sum(x[neg]*dx[neg]); sumP=sum(x[!neg]*dx[!neg])
            # Check if the MIC is zero
            if(abs(sumN-sumP)<tolMIC) { 
                out=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
            } else out=c(FALSE,FALSE,FALSE,sumN<sumP,sumN>sumP,FALSE)
          } else out=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
        }  
      }  
    }
  }
  out
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


randDFP <- function(nbins) {
#  Randomly samples channel completion times under the assumption
#   that low RTs dominate high RTs

    # Only use half of the bins so that at most all bins are possible 
    #  samples for the low channel completion times
    nH <- ceiling(nbins/2)

    # Sample high channel completion times from a Dirichlet process
    H1 <- rdirichlet(1, rep(1/nH,nH))
    H2 <- rdirichlet(1, rep(1/nH,nH))

    # Sample low channel completion times by convolving high channel times
    #  with a random variable with a density sampled from a Dirichlet process
    L1 <- c(dconvolve(rdirichlet(1, rep(2/(nbins-nH),nbins-nH)), H1), 0)
    L2 <- c(dconvolve(rdirichlet(1, rep(2/(nbins-nH),nbins-nH)), H2), 0)

    # Pad the High channel bins with zero they are defined over the same 
    #  lengt vector as the low channel
    H1 <- c(H1, rep(0,nbins/2))
    H2 <- c(H2, rep(0,nbins/2))

    return(list(H1=H1, H2=H2, L1=L1, L2=L2))
}



samplePriorSIC <- function(ns, bins, priorsic) {
# Samples distributions for HH, HL, LH, LL based on the prior distribution over SICs
#  ns:      number of samples to generate
#  bins:    bins over which the distributions are defined

    nbins <- length(bins)-1
    if (sum(priorsic) > 1) {
        warning("Sum of prior density over SICs is larger than 1.\n")
    }
    priorsic <- cumsum(priorsic)

    priorlist <- vector("list",4)
    priorarray <- array(NA, dim=c(ns, 4, nbins))

    for ( i in 1:ns ) {
        model <- runif(1)
        if (model > 10) { 
            # Wiener Coactive Model
            cbins <- bins[1:(nbins-1)] + diff(bins)/2
            cbins <- cbins[1:(nbins-1)]
            rates <- sort(runif(3,min=1e-5,max=.3), decreasing=T)
            while(thresh > 0 ) {thresh <- rnorm(1, mean=44, sd=15)}
            rates <- c(rates, rates[2] + rates[3] - rates[1])
            priorarray[i,1,] <- ns*c(dinvGauss(cbins, nu=thresh/rates[[1]], lambda=thresh^2/4), 
                            ns*  1-pinvGauss(bins[length(bins)-1],nu=thresh/rates[[1]], lambda=thresh^2/4) )
            priorarray[i,2,] <- ns*c(dinvGauss(cbins, nu=thresh/rates[[2]], lambda=thresh^2/4), 
                            ns*  1-pinvGauss(bins[length(bins)-1],nu=thresh/rates[[2]], lambda=thresh^2/4) )
            priorarray[i,3,] <- ns*c(dinvGauss(cbins, nu=thresh/rates[[3]], lambda=thresh^2/4), 
                            ns*  1-pinvGauss(bins[length(bins)-1],nu=thresh/rates[[3]], lambda=thresh^2/4) )
            priorarray[i,4,] <- ns*c(dinvGauss(cbins, nu=thresh/rates[[4]], lambda=thresh^2/4), 
                                1-pinvGauss(bins[length(bins)-1],nu=thresh/rates[[4]], lambda=thresh^2/4) )
        }
        else {
            # Sample chanel completion times
            rdfp <- randDFP(nbins) 
            if (model < priorsic[1] ) { 
                # Serial OR
                p <- runif(1) # proabability that channel 1 is first
                priorarray[i,1,] <- p*rdfp$H1 + (1-p)*rdfp$H2
                priorarray[i,2,] <- p*rdfp$H1 + (1-p)*rdfp$L2
                priorarray[i,3,] <- p*rdfp$L1 + (1-p)*rdfp$H2
                priorarray[i,4,] <- p*rdfp$L1 + (1-p)*rdfp$L2
            } else if (model < priorsic[2]) { 
                # Serial AND
                priorarray[i,1,] <- rebin(dconvolve(rdfp$H1, rdfp$H2))
                priorarray[i,2,] <- rebin(dconvolve(rdfp$H1, rdfp$L2))
                priorarray[i,3,] <- rebin(dconvolve(rdfp$L1, rdfp$H2))
                priorarray[i,4,] <- rebin(dconvolve(rdfp$L1, rdfp$L2))
            } else if (model < priorsic[3]) { 
                # Parallel OR
                FH1 <- cumsum(rdfp$H1)
                FL1 <- cumsum(rdfp$L1)
                FH2 <- cumsum(rdfp$H2)
                FL2 <- cumsum(rdfp$L2)
                priorarray[i,1,] <- makePDF(rdfp$H1*(1-FH2) + (1-FH1)*rdfp$H2 + rdfp$H1*rdfp$H2)
                priorarray[i,2,] <- makePDF(rdfp$H1*(1-FL2) + (1-FH1)*rdfp$L2 + rdfp$H1*rdfp$L2)
                priorarray[i,3,] <- makePDF(rdfp$L1*(1-FH2) + (1-FL1)*rdfp$H2 + rdfp$L1*rdfp$H2)
                priorarray[i,4,] <- makePDF(rdfp$L1*(1-FL2) + (1-FL1)*rdfp$L2 + rdfp$L1*rdfp$L2)
            } else if (model < priorsic[4]) { 
                # Parallel AND
                FH1 <- cumsum(rdfp$H1)
                FL1 <- cumsum(rdfp$L1)
                FH2 <- cumsum(rdfp$H2)
                FL2 <- cumsum(rdfp$L2)
                priorarray[i,1,] <- rdfp$H1*FH2 + FH1*rdfp$H2 - rdfp$H1*rdfp$H2
                priorarray[i,2,] <- rdfp$H1*FL2 + FH1*rdfp$L2 - rdfp$H1*rdfp$L2
                priorarray[i,3,] <- rdfp$L1*FH2 + FL1*rdfp$H2 - rdfp$L1*rdfp$H2
                priorarray[i,4,] <- rdfp$L1*FL2 + FL1*rdfp$L2 - rdfp$L1*rdfp$L2
            }
        }
    }
    return(priorarray)
}



samplePostSIC <- function(ns, binn, priorp, priorinfluence) {
# Samples from the posterior distributions for HH, HL, LH, LL 
#  ns:      number of samples to generate
#  binn:    counts in each bin in the data
#  priorp:  samples from the prior distribution
#  priorinfluence:  weight of the prior distribution
    dimprior <- dim(priorp)
    ns <- dimprior[1]
    out <- array(NA, dimprior)
    priorparr <- priorinfluence*priorp + aperm(array(rep(c(binn[,'HH'], binn[,'HL'], binn[,'LH'], binn[,'LL']),ns), c(dimprior[3], 4,ns)))
    for ( i in 1:ns ) {
        for ( j in 1:4) {
            out[i, j, ] <- rdirichlet(1, priorparr[i,j,])
        }
    }
    return( out )
}


dconvolve <- function(x,y) {
#  Function to compute a discrete convolution
    n <- length(x)+ length(y) -1
    if (length(x) < n ) x <- c(x, rep(0, n-length(x)))
    if (length(y) < n ) y <- c(y, rep(0, n-length(y)))

    out <- rep(0,n) 
    for( i in 1:n ) {
        out[i] =  sum(x[1:i] * y[i:1])
    }
    return(out)
}


makePDF <- function(x) {
# Remove small values from the bins and re-normalize
    x[ x < 1E-10 ] <- 0
    x / sum(x)
}


rebin <- function(x) {
    even <- 2*(1:(floor(length(x)/2)))
    c(x[1], x[even]+x[even+1])
}

samplePriorPDF <- function(ns, alpha) {
    return(array(as.vector(rdirichlet(4*ns,alpha)),dim=c(ns,4,length(alpha))))
}

samplePostPDF <- function(ns, binn, alpha) {
    return(array(rbind(rdirichlet(ns,binn[,1]+alpha),rdirichlet(ns,binn[,2]+alpha),
                rdirichlet(ns,binn[,3]+alpha),rdirichlet(ns,binn[,4]+alpha)),
                dim=c(ns,4,length(alpha))))
}
