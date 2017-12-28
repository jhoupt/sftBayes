data {
  int<lower=0> NHH;
  int<lower=0> NHL;
  int<lower=0> NLH;
  int<lower=0> NLL;
  real<lower=0> yHH[NHH];
  real<lower=0> yHL[NHL];
  real<lower=0> yLH[NLH];
  real<lower=0> yLL[NLL];
}
transformed data {
  real logpi;
  logpi <- log(pi());
}
parameters {
  real<lower=0> ratel;
  real<lower=0> nu;
  real<lower=0> threshold;
}
transformed parameters {
  real<lower=0> rateh;
  real logthresh;
  real summandsHH[NHH];
  real summandsHL[NHL];
  real summandsLH[NLH];
  real summandsLL[NLL];

  logthresh <- log(threshold);
  rateh <- ratel + nu;

  for (n in 1:NHH) {
    summandsHH[n] <- logthresh - log2() - .5 * (logpi +3*log(yHH[n])) - square((rateh+rateh)*yHH[n] - threshold)/(4*yHH[n]);
  }
  for (n in 1:NHL) {
    summandsHL[n] <- logthresh - log2() - .5 * (logpi +3*log(yHL[n])) - square((rateh+ratel)*yHL[n] - threshold)/(4*yHL[n]);
  }
  for (n in 1:NLH) {
    summandsLH[n] <- logthresh - log2() - .5 * (logpi +3*log(yLH[n])) - square((rateh+ratel)*yLH[n] - threshold)/(4*yLH[n]);
  }
  for (n in 1:NLL) {
    summandsLL[n] <- logthresh - log2() - .5 * (logpi +3*log(yLL[n])) - square((ratel+ratel)*yLL[n] - threshold)/(4*yLL[n]);
  }
}
model {
  threshold ~ gamma(4,.01);
  ratel ~ gamma(4,.01);
  nu ~ exponential(100);
  increment_log_prob( sum(summandsHH) + sum(summandsHL) + sum(summandsLH) + sum(summandsLL) );
}
