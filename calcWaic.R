args <- commandArgs(TRUE)
modelname <- args[1]
sn <- args[2]

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


# Read in and prepare the data 
datamodels <- c("SEROR")
#datamodels <- c("COA", "PARAND", "PAROR", "SERAND", "SEROR")

datamvec <- c()
snvec <- c()
waic.out <- c()
for(datam in datamodels) {
  log_lik <- c()
  for(chainid in 1:3) {
    filename <- paste(datam, "_n100_s", sn, "_c", chainid, ".txt", sep="")
    system(paste("sed -i '$ d'", filename))
    #badFile <- tryCatch( { chain <- head(read.csv(filename, comment.char='#'), -1); badFile <-FALSE }, 
    badFile <- tryCatch( 
      { 
        header <- scan(filename,comment.char='#', what=character(), sep=",", skip=31, nlines=10) 
        ncol <- length(header)
        chain <- matrix(scan(filename, comment.char='#', sep=",", skip=38), ncol=ncol, byrow=T)
        badFile <-FALSE 
      }, 
      warning = function(w) { cat(filename, "not found!\n"); return(TRUE) } ,
      error   = function(e) { cat(filename, "not found!\n"); return(TRUE) } 
    )
    if(! badFile) { 
      summandcols <- grep("summands", header)
      badrows <- which(  apply(is.na(chain[,summandcols]), 1, any) )
      if(length(badrows > 0) ) {
        log_lik <- rbind(log_lik, chain[-badrows,summandcols] )
      } else {
        log_lik <- rbind(log_lik, chain[,summandcols] )
      }
    }
  }
  datamvec <- c(datamvec, datam)
  snvec <- c(snvec, sn)
  if( !is.null(log_lik) ) {
    waic.out <- c(waic.out, waic(log_lik)$waic )
    modelvec <- rep(modelname, length(waic.out))
    outdata <- data.frame(Model=modelvec, Data=datamvec, Subject=snvec, WAIC=waic.out)
    outfile <- paste(modelvec, "_", datam, "_s", sn, ".csv", sep="")
    write.csv(outdata, file=outfile, row.names=FALSE)
  } 
  datamvec <- c()
  snvec <- c()
  waic.out <- c()
}
