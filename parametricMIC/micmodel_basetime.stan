functions {
  real shifted_wald_lpdf(real x, real gamma, real alpha, real theta){
    real tmp1;
    real tmp2;
    tmp1 = log(alpha) - .5 * (log2() + log(pi()) + 3 * log(x - theta));
    tmp2 = -square(alpha - gamma * (x - theta))/(2 * (x - theta));
    return tmp1 + tmp2;
  }
}
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

  vector[nSubjects] p_mic_neg;
  vector[nSubjects] p_mic_pos;

  vector[nSubjects] p_A;
  vector[nSubjects] p_B;
  vector[nSubjects] p_C;

  vector<lower=0>[nSubjects] gammaHH;
  vector<lower=0>[nSubjects] gammaHL;
  vector<lower=0>[nSubjects] gammaLH;
  vector<lower=0>[nSubjects] gammaLL;

  vector<lower=0,upper=500>[nSubjects] theta;
}
transformed parameters { 
  vector<lower=0>[nSubjects] mic_magnitude_pos;
  vector<lower=0>[nSubjects] mic_magnitude_neg;
  vector<upper=0>[nSubjects] A;
  vector<lower=0>[nSubjects] B;
  vector<lower=0>[nSubjects] C;

                                        
  // FOR GAMMA PRIOR
  //mic_magnitude = 200 * p_mic;
  //A = -200 * p_A;
  //B = 200 * p_B;
  //C = 400 * p_C;

  // FOR NORMAL PRIOR
  mic_magnitude_pos =  100 + 50 * p_mic_pos;
  mic_magnitude_neg =  100 + 50 * p_mic_neg;
  A   = -100 + 50 * p_A;
  B   =  100 + 50 * p_B;
  C   =  1000 + 100 * p_C;

}
model {
  real logpr[3];

  vector[nSubjects] muHH_pos;
  vector[nSubjects] muHL_pos;
  vector[nSubjects] muLH_pos;
  vector[nSubjects] muLL_pos;

  vector[nSubjects] muHH_neg;
  vector[nSubjects] muHL_neg;
  vector[nSubjects] muLH_neg;
  vector[nSubjects] muLL_neg;

  vector[nSubjects] muHH_0;
  vector[nSubjects] muHL_0;
  vector[nSubjects] muLH_0;
  vector[nSubjects] muLL_0;

  vector[nSubjects] alphaHH_pos;
  vector[nSubjects] alphaHL_pos;
  vector[nSubjects] alphaLH_pos;
  vector[nSubjects] alphaLL_pos;

  vector[nSubjects] alphaHH_neg;
  vector[nSubjects] alphaHL_neg;
  vector[nSubjects] alphaLH_neg;
  vector[nSubjects] alphaLL_neg;
  
  vector[nSubjects] alphaHH_0;
  vector[nSubjects] alphaHL_0;
  vector[nSubjects] alphaLH_0;
  vector[nSubjects] alphaLL_0;

  muHH_pos =  .25 * mic_magnitude_pos + .5 * A - .5 * B + .5 * C;
  muHL_pos = -.25 * mic_magnitude_pos - .5 * A - .5 * B + .5 * C;
  muLH_pos = -.25 * mic_magnitude_pos + .5 * A + .5 * B + .5 * C;
  muLL_pos =  .25 * mic_magnitude_pos - .5 * A + .5 * B + .5 * C;

  muHH_neg = -.25 * mic_magnitude_neg + .5 * A - .5 * B + .5 * C;
  muHL_neg =  .25 * mic_magnitude_neg - .5 * A - .5 * B + .5 * C;
  muLH_neg =  .25 * mic_magnitude_neg + .5 * A + .5 * B + .5 * C;
  muLL_neg = -.25 * mic_magnitude_neg - .5 * A + .5 * B + .5 * C;

  muHH_0   =            + .5 * A - .5 * B + .5 * C;
  muHL_0   =            - .5 * A - .5 * B + .5 * C;
  muLH_0   =            + .5 * A + .5 * B + .5 * C;
  muLL_0   =            - .5 * A + .5 * B + .5 * C;

  alphaHH_pos = muHH_pos .* gammaHH;
  alphaHL_pos = muHL_pos .* gammaHL;
  alphaLH_pos = muLH_pos .* gammaLH;
  alphaLL_pos = muLL_pos .* gammaLL;

  alphaHH_neg = muHH_neg .* gammaHH;
  alphaHL_neg = muHL_neg .* gammaHL;
  alphaLH_neg = muLH_neg .* gammaLH;
  alphaLL_neg = muLL_neg .* gammaLL;

  alphaHH_0 = muHH_0 .* gammaHH;
  alphaHL_0 = muHL_0 .* gammaHL;
  alphaLH_0 = muLH_0 .* gammaLH;
  alphaLL_0 = muLL_0 .* gammaLL;





  for ( tr in 1:nHH) {
     logpr[1] = log(modelProb_subj[subjectHH[tr], 1]) 
                + shifted_wald_lpdf(rtHH[tr] | gammaHH[subjectHH[tr]],
                                               alphaHH_pos[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     logpr[2] = log(modelProb_subj[subjectHH[tr], 3]) 
                + shifted_wald_lpdf(rtHH[tr] | gammaHH[subjectHH[tr]],
                                               alphaHH_neg[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     logpr[3] = log(modelProb_subj[subjectHH[tr], 2]) 
                + shifted_wald_lpdf(rtHH[tr] | gammaHH[subjectHH[tr]],
                                               alphaHH_0[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     target += log_sum_exp(logpr);
  }
  for ( tr in 1:nHL) {
     logpr[1] = log(modelProb_subj[subjectHL[tr], 1]) 
                + shifted_wald_lpdf(rtHL[tr] | gammaHL[subjectHL[tr]], 
                                               alphaHL_pos[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     logpr[2] = log(modelProb_subj[subjectHL[tr], 3]) 
                + shifted_wald_lpdf(rtHL[tr] | gammaHL[subjectHL[tr]],
                                               alphaHL_neg[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     logpr[3] = log(modelProb_subj[subjectHL[tr], 2]) 
                + shifted_wald_lpdf(rtHL[tr] | gammaHL[subjectHL[tr]],
                                               alphaHL_0[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     target += log_sum_exp(logpr);
  }
  for ( tr in 1:nLH) {
     logpr[1] = log(modelProb_subj[subjectLH[tr], 1]) 
                + shifted_wald_lpdf(rtLH[tr] | gammaLH[subjectLH[tr]], 
                                               alphaLH_pos[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     logpr[2] = log(modelProb_subj[subjectLH[tr], 3]) 
                + shifted_wald_lpdf(rtLH[tr] | gammaLH[subjectLH[tr]],
                                               alphaLH_neg[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     logpr[3] = log(modelProb_subj[subjectLH[tr], 2]) 
                + shifted_wald_lpdf(rtLH[tr] | gammaLH[subjectLH[tr]],
                                               alphaLH_0[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     target += log_sum_exp(logpr);
  }
  for ( tr in 1:nLL) {
     logpr[1] = log(modelProb_subj[subjectLL[tr], 1]) 
                + shifted_wald_lpdf(rtLL[tr] | gammaLL[subjectLL[tr]], 
                                               alphaLL_pos[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     logpr[2] = log(modelProb_subj[subjectLL[tr], 3]) 
                + shifted_wald_lpdf(rtLL[tr] | gammaLL[subjectLL[tr]],
                                               alphaLL_neg[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     logpr[3] = log(modelProb_subj[subjectLL[tr], 2]) 
                + shifted_wald_lpdf(rtLL[tr] | gammaLL[subjectLL[tr]],
                                               alphaLL_0[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     target += log_sum_exp(logpr);
  }

  //p_mic ~ gamma(2,4);
  //p_A ~ gamma(2,4);
  //p_B ~ gamma(2,4);
  //p_C ~ gamma(2,4);

  p_mic_neg ~ normal(0,1);
  p_mic_pos ~ normal(0,1);
  p_A ~ normal(0,1);
  p_B ~ normal(0,1);
  p_C ~ normal(0,1);

  gammaHH ~ gamma(1,1);
  gammaHL ~ gamma(1,1);
  gammaLH ~ gamma(1,1);
  gammaLL ~ gamma(1,1);

  modelProb_group ~ dirichlet(modelPrior);
  for (sj in 1:nSubjects) {
     modelProb_subj[sj] ~ dirichlet(modelProb_group);
  }
}


