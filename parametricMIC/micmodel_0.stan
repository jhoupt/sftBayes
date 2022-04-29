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
  vector[nSubjects] p_A;
  vector[nSubjects] p_B;
  vector[nSubjects] p_C;

  vector<lower=0>[nSubjects] gammaHH_0;
  vector<lower=0>[nSubjects] gammaHL_0;
  vector<lower=0>[nSubjects] gammaLH_0;
  vector<lower=0>[nSubjects] gammaLL_0;

  vector<lower=0,upper=500>[nSubjects] theta;
}
transformed parameters { 
  vector<upper=0>[nSubjects] A;
  vector<lower=0>[nSubjects] B;
  vector<lower=0>[nSubjects] C;

  vector<lower=0>[nSubjects] muHH_0;
  vector<lower=0>[nSubjects] muHL_0;
  vector<lower=0>[nSubjects] muLH_0;
  vector<lower=0>[nSubjects] muLL_0;

  vector<lower=0>[nSubjects] alphaHH_0;
  vector<lower=0>[nSubjects] alphaHL_0;
  vector<lower=0>[nSubjects] alphaLH_0;
  vector<lower=0>[nSubjects] alphaLL_0;
                                        
  // FOR NORMAL PRIOR
  A   = -50 + 25 * p_A;
  B   =  50 + 25 * p_B;
  C   =  1000 + 100 * p_C;

  muHH_0   =            + .5 * A - .5 * B + .5 * C;
  muHL_0   =            - .5 * A - .5 * B + .5 * C;
  muLH_0   =            + .5 * A + .5 * B + .5 * C;
  muLL_0   =            - .5 * A + .5 * B + .5 * C;

  alphaHH_0 = muHH_0 .* gammaHH_0;
  alphaHL_0 = muHL_0 .* gammaHL_0;
  alphaLH_0 = muLH_0 .* gammaLH_0;
  alphaLL_0 = muLL_0 .* gammaLL_0;
}
model {
  real logpr[3];
  for ( tr in 1:nHH) {
     logpr[3] = log(1) 
                + shifted_wald_lpdf(rtHH[tr] | gammaHH_0[subjectHH[tr]],
                                               alphaHH_0[subjectHH[tr]],
                                               theta[subjectHH[tr]]);
     target += logpr[3];
  }
  for ( tr in 1:nHL) {
     logpr[3] = log(1)
                + shifted_wald_lpdf(rtHL[tr] | gammaHL_0[subjectHL[tr]],
                                               alphaHL_0[subjectHL[tr]],
                                               theta[subjectHL[tr]]);
     target += logpr[3];
  }
  for ( tr in 1:nLH) {
     logpr[3] = log(1)
                + shifted_wald_lpdf(rtLH[tr] | gammaLH_0[subjectLH[tr]],
                                               alphaLH_0[subjectLH[tr]],
                                               theta[subjectLH[tr]]);
     target += logpr[3];
  }
  for ( tr in 1:nLL) {
     logpr[3] = log(1)
                + shifted_wald_lpdf(rtLL[tr] | gammaLL_0[subjectLL[tr]],
                                               alphaLL_0[subjectLL[tr]],
                                               theta[subjectLL[tr]]);
     target += logpr[3];
  }


  p_A ~ normal(0,1);
  p_B ~ normal(0,1);
  p_C ~ normal(0,1);

  gammaHH_0 ~ gamma(1,1);
  gammaHL_0 ~ gamma(1,1);
  gammaLH_0 ~ gamma(1,1);
  gammaLL_0 ~ gamma(1,1);

}


