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
  vector[nSubjects] p_mic_neg;

  vector[nSubjects] p_A;
  vector[nSubjects] p_B;
  vector[nSubjects] p_C;

  vector<lower=0>[nSubjects] gammaHH_neg;
  vector<lower=0>[nSubjects] gammaHL_neg;
  vector<lower=0>[nSubjects] gammaLH_neg;
  vector<lower=0>[nSubjects] gammaLL_neg;

  vector<lower=0,upper=500>[nSubjects] theta;
}
transformed parameters { 
  vector<lower=0>[nSubjects] mic_magnitude_neg;
  vector<upper=0>[nSubjects] A;
  vector<lower=0>[nSubjects] B;
  vector<lower=0>[nSubjects] C;

  vector<lower=0>[nSubjects] muHH_neg;
  vector<lower=0>[nSubjects] muHL_neg;
  vector<lower=0>[nSubjects] muLH_neg;
  vector<lower=0>[nSubjects] muLL_neg;

  vector<lower=0>[nSubjects] alphaHH_neg;
  vector<lower=0>[nSubjects] alphaHL_neg;
  vector<lower=0>[nSubjects] alphaLH_neg;
  vector<lower=0>[nSubjects] alphaLL_neg;
  
  // FOR NORMAL PRIOR
  mic_magnitude_neg =  100 + 50 * p_mic_neg;
  A   = -100 + 50 * p_A;
  B   =  100 + 50 * p_B;
  C   =  1000 + 100 * p_C;

  muHH_neg = -.25 * mic_magnitude_neg + .5 * A - .5 * B + .5 * C;
  muHL_neg =  .25 * mic_magnitude_neg - .5 * A - .5 * B + .5 * C;
  muLH_neg =  .25 * mic_magnitude_neg + .5 * A + .5 * B + .5 * C;
  muLL_neg = -.25 * mic_magnitude_neg - .5 * A + .5 * B + .5 * C;

  alphaHH_neg = muHH_neg .* gammaHH_neg;
  alphaHL_neg = muHL_neg .* gammaHL_neg;
  alphaLH_neg = muLH_neg .* gammaLH_neg;
  alphaLL_neg = muLL_neg .* gammaLL_neg;
}
model {
  real logpr[3];
  for ( tr in 1:nHH) {
     logpr[2] = log(1)
                + shifted_wald_lpdf(rtHH[tr] | gammaHH_neg[subjectHH[tr]],
                                               alphaHH_neg[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     target += logpr[2];
  }
  for ( tr in 1:nHL) {
     logpr[2] = log(1)
                + shifted_wald_lpdf(rtHL[tr] | gammaHL_neg[subjectHL[tr]],
                                               alphaHL_neg[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     target += logpr[2];
  }
  for ( tr in 1:nLH) {
     logpr[2] = log(1)
                + shifted_wald_lpdf(rtLH[tr] | gammaLH_neg[subjectLH[tr]],
                                               alphaLH_neg[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     target += logpr[2];
  }
  for ( tr in 1:nLL) {
     logpr[2] = log(1)
                + shifted_wald_lpdf(rtLL[tr] | gammaLL_neg[subjectLL[tr]],
                                               alphaLL_neg[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     target += logpr[2];
  }

  p_mic_neg ~ normal(0,1);
  p_A ~ normal(0,1);
  p_B ~ normal(0,1);
  p_C ~ normal(0,1);

  gammaHH_neg ~ gamma(1,1);
  gammaHL_neg ~ gamma(1,1);
  gammaLH_neg ~ gamma(1,1);
  gammaLL_neg ~ gamma(1,1);

}


