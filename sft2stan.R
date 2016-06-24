sft2stan <- function(sftData, DUMP=FALSE) {
   require(rstan)
   nSubjects <- length(unique(sftData$Subject))
   sftData$Subject <- as.numeric(factor(sftData$Subject))

   subjectHH <- with(sftData, Subject[Channel1==2 & Channel2==2 & Correct] )
   subjectHL <- with(sftData, Subject[Channel1==2 & Channel2==1 & Correct] )
   subjectLH <- with(sftData, Subject[Channel1==1 & Channel2==2 & Correct] )
   subjectLL <- with(sftData, Subject[Channel1==1 & Channel2==1 & Correct] )

   rtHH <- with(sftData, RT[Channel1==2 & Channel2==2 & Correct] )
   rtHL <- with(sftData, RT[Channel1==2 & Channel2==1 & Correct] )
   rtLH <- with(sftData, RT[Channel1==1 & Channel2==2 & Correct] )
   rtLL <- with(sftData, RT[Channel1==1 & Channel2==1 & Correct] )

   nHH <- length(rtHH)
   nHL <- length(rtHL)
   nLH <- length(rtLH)
   nLL <- length(rtLL)
  
   modelPrior <- array(c(.25, .50, .25))

   if( DUMP ) {
      stan_rdump(c("nSubjects", "nHH", "nHL", "nLH", "nLL", "subjectHH", "subjectHL", "subjectLH", "subjectLL", "rtHH", "rtHL", "rtLH", "rtLL", "modelPrior"), file="sftStanData.Rdmp")
   }
   return( list(nSubjects=nSubjects, nHH=nHH, nHL=nHL, nLH=nLH, nLL=nLL, subjectHH=subjectHH, subjectHL=subjectHL, subjectLH=subjectLH, subjectLL=subjectLL, rtHH=rtHH, rtHL=rtHL, rtLH=rtLH, rtLL=rtLL, modelPrior=modelPrior))
}
