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
  vector[nSubjects] p_mic_pos;

  vector[nSubjects] p_A;
  vector[nSubjects] p_B;
  vector[nSubjects] p_C;

  vector<lower=0>[nSubjects] gammaHH_pos;
  vector<lower=0>[nSubjects] gammaHL_pos;
  vector<lower=0>[nSubjects] gammaLH_pos;
  vector<lower=0>[nSubjects] gammaLL_pos;

  vector<lower=0,upper=500>[nSubjects] theta;
}
transformed parameters { 
  vector<lower=0>[nSubjects] mic_magnitude_pos;
  vector<upper=0>[nSubjects] A;
  vector<lower=0>[nSubjects] B;
  vector<lower=0>[nSubjects] C;

  vector<lower=0>[nSubjects] muHH_pos;
  vector<lower=0>[nSubjects] muHL_pos;
  vector<lower=0>[nSubjects] muLH_pos;
  vector<lower=0>[nSubjects] muLL_pos;


  vector<lower=0>[nSubjects] alphaHH_pos;
  vector<lower=0>[nSubjects] alphaHL_pos;
  vector<lower=0>[nSubjects] alphaLH_pos;
  vector<lower=0>[nSubjects] alphaLL_pos;

  // FOR NORMAL PRIOR
  mic_magnitude_pos =  100 + 50 * p_mic_pos;
  A   = -100 + 50 * p_A;
  B   =  100 + 50 * p_B;
  C   =  1000 + 100 * p_C;

  muHH_pos =  .25 * mic_magnitude_pos + .5 * A - .5 * B + .5 * C;
  muHL_pos = -.25 * mic_magnitude_pos - .5 * A - .5 * B + .5 * C;
  muLH_pos = -.25 * mic_magnitude_pos + .5 * A + .5 * B + .5 * C;
  muLL_pos =  .25 * mic_magnitude_pos - .5 * A + .5 * B + .5 * C;

  alphaHH_pos = muHH_pos .* gammaHH_pos;
  alphaHL_pos = muHL_pos .* gammaHL_pos;
  alphaLH_pos = muLH_pos .* gammaLH_pos;
  alphaLL_pos = muLL_pos .* gammaLL_pos;

}
model {
  real logpr[3];
  for ( tr in 1:nHH) {
     logpr[1] = log(1)
                + shifted_wald_lpdf(rtHH[tr] | gammaHH_pos[subjectHH[tr]],
                                               alphaHH_pos[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     target += logpr[1];
  }
  for ( tr in 1:nHL) {
     logpr[1] = log(1)
                + shifted_wald_lpdf(rtHL[tr] | gammaHL_pos[subjectHL[tr]], 
                                               alphaHL_pos[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     target += logpr[1];
  }
  for ( tr in 1:nLH) {
     logpr[1] = log(1)
                + shifted_wald_lpdf(rtLH[tr] | gammaLH_pos[subjectLH[tr]], 
                                               alphaLH_pos[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     target += logpr[1];
  }
  for ( tr in 1:nLL) {
     logpr[1] = log(1)
                + shifted_wald_lpdf(rtLL[tr] | gammaLL_pos[subjectLL[tr]], 
                                               alphaLL_pos[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     target += logpr[1];
  }

  p_mic_pos ~ normal(0,1);
  p_A ~ normal(0,1);
  p_B ~ normal(0,1);
  p_C ~ normal(0,1);

  gammaHH_pos ~ gamma(1,1);
  gammaHL_pos ~ gamma(1,1);
  gammaLH_pos ~ gamma(1,1);
  gammaLL_pos ~ gamma(1,1);

}


