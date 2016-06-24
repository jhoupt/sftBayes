data {
   int<lower=1> nSubjects;
   int<lower=1> nHH;
   int<lower=1> nHL;
   int<lower=1> nLH;
   int<lower=1> nLL;
   int subjectHH[nHH];
   int subjectHL[nHL];
   int subjectLH[nLH];
   int subjectLL[nLL];
   real<lower=0> rtHH[nHH];
   real<lower=0> rtHL[nHL];
   real<lower=0> rtLH[nLH];
   real<lower=0> rtLL[nLL];
   simplex[3] modelPrior;

}
parameters {
   simplex[3] modelProb_group;
   simplex[3] modelProb_subj[nSubjects];

   vector[nSubjects] p_mic;

   vector[nSubjects] p_A;
   vector[nSubjects] p_B;
   vector[nSubjects] p_C;

   vector<lower=0>[nSubjects] rateHH;
   vector<lower=0>[nSubjects] rateHL;
   vector<lower=0>[nSubjects] rateLH;
   vector<lower=0>[nSubjects] rateLL;

}
transformed parameters { 
   vector<lower=0>[nSubjects] mic;
   vector<upper=0>[nSubjects] A;
   vector<lower=0>[nSubjects] B;
   vector<lower=0>[nSubjects] C;

   vector<lower=0>[nSubjects] muHH_pos;
   vector<lower=0>[nSubjects] muHL_pos;
   vector<lower=0>[nSubjects] muLH_pos;
   vector<lower=0>[nSubjects] muLL_pos;

   vector<lower=0>[nSubjects] muHH_neg;
   vector<lower=0>[nSubjects] muHL_neg;
   vector<lower=0>[nSubjects] muLH_neg;
   vector<lower=0>[nSubjects] muLL_neg;

   vector<lower=0>[nSubjects] muHH_0;
   vector<lower=0>[nSubjects] muHL_0;
   vector<lower=0>[nSubjects] muLH_0;
   vector<lower=0>[nSubjects] muLL_0;

   vector<lower=0>[nSubjects] shapeHH_pos;
   vector<lower=0>[nSubjects] shapeHL_pos;
   vector<lower=0>[nSubjects] shapeLH_pos;
   vector<lower=0>[nSubjects] shapeLL_pos;

   vector<lower=0>[nSubjects] shapeHH_neg;
   vector<lower=0>[nSubjects] shapeHL_neg;
   vector<lower=0>[nSubjects] shapeLH_neg;
   vector<lower=0>[nSubjects] shapeLL_neg;
   
   vector<lower=0>[nSubjects] shapeHH_0;
   vector<lower=0>[nSubjects] shapeHL_0;
   vector<lower=0>[nSubjects] shapeLH_0;
   vector<lower=0>[nSubjects] shapeLL_0;
                                         
   // FOR GAMMA PRIOR
   //mic <- 200 * p_mic;
   //A <- -200 * p_A;
   //B <- 200 * p_B;
   //C <- 400 * p_C;

   // FOR NORMAL PRIOR
   mic <-  100 + 50 * p_mic;
   A   <- -100 + 50 * p_A;
   B   <-  100 + 50 * p_B;
   C   <-  400 + 100 * p_C;



   muHH_pos <-  .25 * mic + .5 * A - .5 * B + .5 * C;
   muHL_pos <- -.25 * mic - .5 * A - .5 * B + .5 * C;
   muLH_pos <- -.25 * mic + .5 * A + .5 * B + .5 * C;
   muLL_pos <-  .25 * mic - .5 * A + .5 * B + .5 * C;

   muHH_neg <- -.25 * mic + .5 * A - .5 * B + .5 * C;
   muHL_neg <-  .25 * mic - .5 * A - .5 * B + .5 * C;
   muLH_neg <-  .25 * mic + .5 * A + .5 * B + .5 * C;
   muLL_neg <- -.25 * mic - .5 * A + .5 * B + .5 * C;

   muHH_0   <-                + .5 * A - .5 * B + .5 * C;
   muHL_0   <-                - .5 * A - .5 * B + .5 * C;
   muLH_0   <-                + .5 * A + .5 * B + .5 * C;
   muLL_0   <-                - .5 * A + .5 * B + .5 * C;

   shapeHH_pos <- muHH_pos .* rateHH;
   shapeHL_pos <- muHL_pos .* rateHL;
   shapeLH_pos <- muLH_pos .* rateLH;
   shapeLL_pos <- muLL_pos .* rateLL;

   shapeHH_neg <- muHH_neg .* rateHH;
   shapeHL_neg <- muHL_neg .* rateHL;
   shapeLH_neg <- muLH_neg .* rateLH;
   shapeLL_neg <- muLL_neg .* rateLL;

   shapeHH_0 <- muHH_0 .* rateHH;
   shapeHL_0 <- muHL_0 .* rateHL;
   shapeLH_0 <- muLH_0 .* rateLH;
   shapeLL_0 <- muLL_0 .* rateLL;
}
model {
   real logpr[3];
   for ( tr in 1:nHH) {
      //increment_log_prob( gamma_log(rtHH[tr], shapeHH_0[subjectHH[tr]], rateHH[subjectHH[tr]]) );
      logpr[1] <- log(modelProb_subj[subjectHH[tr], 1]) + gamma_log(rtHH[tr], shapeHH_pos[subjectHH[tr]], rateHH[subjectHH[tr]] );
      logpr[2] <- log(modelProb_subj[subjectHH[tr], 3]) + gamma_log(rtHH[tr], shapeHH_neg[subjectHH[tr]], rateHH[subjectHH[tr]] );
      logpr[3] <- log(modelProb_subj[subjectHH[tr], 2]) + gamma_log(rtHH[tr], shapeHH_0[subjectHH[tr]],   rateHH[subjectHH[tr]] );
      increment_log_prob( log_sum_exp(logpr) ) ;
   }
   for ( tr in 1:nHL) {
      //increment_log_prob( gamma_log(rtHL[tr], shapeHL_0[subjectHL[tr]], rateHL[subjectHL[tr]]) ) ;
      logpr[1] <- log(modelProb_subj[subjectHL[tr], 1]) + gamma_log(rtHL[tr], shapeHL_pos[subjectHL[tr]], rateHL[subjectHL[tr]]);
      logpr[2] <- log(modelProb_subj[subjectHL[tr], 3]) + gamma_log(rtHL[tr], shapeHL_neg[subjectHL[tr]], rateHL[subjectHL[tr]]);
      logpr[3] <- log(modelProb_subj[subjectHL[tr], 2]) + gamma_log(rtHL[tr], shapeHL_0[subjectHL[tr]],   rateHL[subjectHL[tr]]);
      increment_log_prob( log_sum_exp(logpr) ) ;
   }
   for ( tr in 1:nLH) {
      //increment_log_prob( gamma_log(rtLH[tr], shapeLH_0[subjectLH[tr]], rateLH[subjectLH[tr]]) ) ;
      logpr[1] <- log(modelProb_subj[subjectLH[tr], 1]) + gamma_log(rtLH[tr], shapeLH_pos[subjectLH[tr]], rateLH[subjectLH[tr]]);
      logpr[2] <- log(modelProb_subj[subjectLH[tr], 3]) + gamma_log(rtLH[tr], shapeLH_neg[subjectLH[tr]], rateLH[subjectLH[tr]]);
      logpr[3] <- log(modelProb_subj[subjectLH[tr], 2]) + gamma_log(rtLH[tr], shapeLH_0[subjectLH[tr]],   rateLH[subjectLH[tr]]);
      increment_log_prob( log_sum_exp(logpr) ) ;
   }
   for ( tr in 1:nLL) {
      //increment_log_prob( gamma_log(rtLL[tr], shapeLL_0[subjectLL[tr]], rateLL[subjectLL[tr]]) ) ;
      logpr[1] <- log(modelProb_subj[subjectLL[tr], 1]) + gamma_log(rtLL[tr], shapeLL_pos[subjectLL[tr]], rateLL[subjectLL[tr]]);
      logpr[2] <- log(modelProb_subj[subjectLL[tr], 3]) + gamma_log(rtLL[tr], shapeLL_neg[subjectLL[tr]], rateLL[subjectLL[tr]]);
      logpr[3] <- log(modelProb_subj[subjectLL[tr], 2]) + gamma_log(rtLL[tr], shapeLL_0[subjectLL[tr]],   rateLL[subjectLL[tr]]);
      increment_log_prob( log_sum_exp(logpr) ) ;
   }

   //p_mic ~ gamma(2,4);
   //p_A ~ gamma(2,4);
   //p_B ~ gamma(2,4);
   //p_C ~ gamma(2,4);

   p_mic ~ normal(0,1);
   p_A ~ normal(0,1);
   p_B ~ normal(0,1);
   p_C ~ normal(0,1);

   //rateHH ~ gamma(1,1);
   //rateHL ~ gamma(1,1);
   //rateLH ~ gamma(1,1);
   //rateLL ~ gamma(1,1);

   modelProb_group ~ dirichlet(modelPrior);
   for (sj in 1:nSubjects) {
      modelProb_subj[sj] ~ dirichlet(modelProb_group);
   }
}
