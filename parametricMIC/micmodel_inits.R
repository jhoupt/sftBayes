micmodel_inits <- function(stanData, DUMP=FALSE) {
   nSubjects <- stanData$nSubjects

   modelProb_subj <- array(.25,  c(nSubjects, 3))
   modelProb_subj[,2] <- .5
   modelProb_group <- c(.25, .5, .25)
   

   A<- rep(NA, nSubjects)
   for (subj in 1:nSubjects) { 
      A[subj] <- with(stanData, -.5 * (mean(rtHL[subjectHL==subj]) -  mean(rtHH[subjectHH==subj]) + mean(rtLL[subjectLL==subj]) -  mean(rtLH[subjectLH==subj]) ) )
   }
   A[A>0] <- -1E2

   B<- rep(NA, nSubjects)
   for (subj in 1:nSubjects) { 
      B[subj] <- with(stanData, .5 * (mean(rtLH[subjectLH==subj]) -  mean(rtHH[subjectHH==subj]) + mean(rtLL[subjectLL==subj]) -  mean(rtHL[subjectHL==subj]) ) )
   }
   B[B<0] <- 1E2

   C<- rep(NA, nSubjects)
   for (subj in 1:nSubjects) { 
      C[subj] <- with(stanData, 2* mean(c(rtHH[subjectHH==subj], rtHL[subjectHL==subj], rtLH[subjectLH==subj], rtLL[subjectLL==subj])) )
   }


   good <- FALSE
   while (! good) {
      good <- TRUE
      p_mic <- abs(rnorm(nSubjects,0,1))
      mic <- 100 + 50 * p_mic

      p_A <- (A + 100) / 50 
      p_B <- (B - 100) / 50 
      p_C <- (C - 400) / 100 

      rateHH <- rgamma(nSubjects,1,1);
      rateHL <- rgamma(nSubjects,1,1);
      rateLH <- rgamma(nSubjects,1,1);
      rateLL <- rgamma(nSubjects,1,1);

      muHH_pos <-  .25 * mic + .5 * A - .5 * B + .5 * C
      muHL_pos <- -.25 * mic - .5 * A - .5 * B + .5 * C
      muLH_pos <- -.25 * mic + .5 * A + .5 * B + .5 * C
      muLL_pos <-  .25 * mic - .5 * A + .5 * B + .5 * C
                                                       
      muHH_neg <- -.25 * mic + .5 * A - .5 * B + .5 * C
      muHL_neg <-  .25 * mic - .5 * A - .5 * B + .5 * C
      muLH_neg <-  .25 * mic + .5 * A + .5 * B + .5 * C
      muLL_neg <- -.25 * mic - .5 * A + .5 * B + .5 * C
                                                                                                 
      muHH_0   <-            + .5 * A - .5 * B + .5 * C
      muHL_0   <-            - .5 * A - .5 * B + .5 * C
      muLH_0   <-            + .5 * A + .5 * B + .5 * C
      muLL_0   <-            - .5 * A + .5 * B + .5 * C

      if (with(stanData, any(muHH_pos < 0) | any(muHL_pos < 0) | any(muLH_pos <0) | any(muLL_pos < 0))) { good <- FALSE}
      if (with(stanData, any(muHH_neg < 0) | any(muHL_neg < 0) | any(muLH_neg <0) | any(muLL_neg < 0))) { good <- FALSE}
      if (with(stanData, any(muHH_0 < 0) | any(muHL_0 < 0) | any(muLH_0 <0) | any(muLL_0 < 0))) { good <- FALSE}

      #if (any(muLH_pos - muHH_pos < 0) | any(muLH_neg - muHH_neg < 0) | any(muLH_0 - muHH_0 < 0) ) { good <- FALSE }
      #if (any(muHL_pos - muHH_pos < 0) | any(muHL_neg - muHH_neg < 0) | any(muHL_0 - muHH_0 < 0) ) { good <- FALSE }
      #if (any(muLL_pos - muLH_pos < 0) | any(muLL_neg - muLH_neg < 0) | any(muLL_0 - muLH_0 < 0) ) { good <- FALSE }
      #if (any(muLL_pos - muHL_pos < 0) | any(muLL_neg - muHL_neg < 0) | any(muLL_0 - muHL_0 < 0) ) { good <- FALSE }

      if(!good) { C <- 1.1*C }    
   } 

   if (DUMP) {
      stan_rdump(c("p_mic", "p_A", "p_B", "p_C", "rateHH", "rateHL", "rateLH", "rateLL",  "modelProb_group", "modelProb_subj"), file="sftStanInits.Rdmp")
   }
   return(list(p_mic=p_mic, p_A=p_A, p_B=p_B, p_C=p_C,
               rateHH=rateHH, rateHL=rateHL, rateLH=rateLH, rateLL=rateLL, 
               modelProb_group=modelProb_group,
               modelProb_subj=modelProb_subj))
}


mic_shiftedwald_inits <- function(stanData, DUMP=FALSE) {
   nSubjects <- stanData$nSubjects

   modelProb_subj <- array(.25,  c(nSubjects, 3))
   modelProb_subj[,2] <- .5
   modelProb_group <- c(.25, .5, .25)
   

   A<- rep(NA, nSubjects)
   B<- rep(NA, nSubjects)
   C<- rep(NA, nSubjects)
   theta <- rep(NA, nSubjects)

   attach(stanData)
   for (subj in 1:nSubjects) { 
      theta[subj] <- min(c(rtHH[subjectHH==subj], rtHL[subjectHL==subj], 
                           rtLH[subjectLH==subj], rtLL[subjectLL==subj]))
      theta[subj] <- .95 * theta[subj]
      A[subj] <- -.5 * (mean(rtHL[subjectHL==subj]) 
                        - mean(rtHH[subjectHH==subj]) 
                        + mean(rtLL[subjectLL==subj]) 
                        - mean(rtLH[subjectLH==subj]))
      B[subj] <-  .5 * (mean(rtLH[subjectLH==subj]) 
                        - mean(rtHH[subjectHH==subj]) 
                        + mean(rtLL[subjectLL==subj]) 
                        - mean(rtHL[subjectHL==subj]))
      C[subj] <- 2 * (mean(c(rtHH[subjectHH==subj], rtHL[subjectHL==subj], 
                            rtLH[subjectLH==subj], rtLL[subjectLL==subj]))
                      - theta[subj])
   }
   detach(stanData)
   A[A>0] <- -1E2
   B[B<0] <- 1E2

   good <- FALSE
   while (! good) {
      good <- TRUE
      p_mic <- abs(rnorm(nSubjects,0,1))
      mic <- 100 + 50 * p_mic

      p_A <- (A + 100) / 50 
      p_B <- (B - 100) / 50 
      p_C <- (C - 400) / 100 

      gammaHH <- rgamma(nSubjects,1,1);
      gammaHL <- rgamma(nSubjects,1,1);
      gammaLH <- rgamma(nSubjects,1,1);
      gammaLL <- rgamma(nSubjects,1,1);

      muHH_pos <-  .25 * mic + .5 * A - .5 * B + .5 * C
      muHL_pos <- -.25 * mic - .5 * A - .5 * B + .5 * C
      muLH_pos <- -.25 * mic + .5 * A + .5 * B + .5 * C
      muLL_pos <-  .25 * mic - .5 * A + .5 * B + .5 * C
                                                       
      muHH_neg <- -.25 * mic + .5 * A - .5 * B + .5 * C
      muHL_neg <-  .25 * mic - .5 * A - .5 * B + .5 * C
      muLH_neg <-  .25 * mic + .5 * A + .5 * B + .5 * C
      muLL_neg <- -.25 * mic - .5 * A + .5 * B + .5 * C
                                                                                                 
      muHH_0   <-            + .5 * A - .5 * B + .5 * C
      muHL_0   <-            - .5 * A - .5 * B + .5 * C
      muLH_0   <-            + .5 * A + .5 * B + .5 * C
      muLL_0   <-            - .5 * A + .5 * B + .5 * C

      if (with(stanData, any(muHH_pos < 0) | any(muHL_pos < 0) | any(muLH_pos <0) | any(muLL_pos < 0))) { good <- FALSE}
      if (with(stanData, any(muHH_neg < 0) | any(muHL_neg < 0) | any(muLH_neg <0) | any(muLL_neg < 0))) { good <- FALSE}
      if (with(stanData, any(muHH_0 < 0) | any(muHL_0 < 0) | any(muLH_0 <0) | any(muLL_0 < 0))) { good <- FALSE}

      #if (any(muLH_pos - muHH_pos < 0) | any(muLH_neg - muHH_neg < 0) | any(muLH_0 - muHH_0 < 0) ) { good <- FALSE }
      #if (any(muHL_pos - muHH_pos < 0) | any(muHL_neg - muHH_neg < 0) | any(muHL_0 - muHH_0 < 0) ) { good <- FALSE }
      #if (any(muLL_pos - muLH_pos < 0) | any(muLL_neg - muLH_neg < 0) | any(muLL_0 - muLH_0 < 0) ) { good <- FALSE }
      #if (any(muLL_pos - muHL_pos < 0) | any(muLL_neg - muHL_neg < 0) | any(muLL_0 - muHL_0 < 0) ) { good <- FALSE }

      if(!good) { C <- 1.1*C }    
   } 

   if (DUMP) {
      stan_rdump(c("p_mic", "p_A", "p_B", "p_C", "gammaHH", "gammaHL", "gammaLH", "gammaLL",  "modelProb_group", "modelProb_subj"), file="sftStanInits.Rdmp")
   }
   return(list(p_mic=p_mic, p_A=p_A, p_B=p_B, p_C=p_C,
               theta=theta,
               gammaHH=gammaHH, gammaHL=gammaHL, 
               gammaLH=gammaLH, gammaLL=gammaLL, 
               modelProb_group=modelProb_group,
               modelProb_subj=modelProb_subj))
}
