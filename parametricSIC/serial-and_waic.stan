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
  real<lower=0> threshold1;
  real<lower=0> threshold2;
}
transformed parameters {
  real<lower=0> rateh;
  real<lower=0> muhh;
  real<lower=0> muhl;
  real<lower=0> mulh;
  real<lower=0> mull;
  real<lower=0> varhh;
  real<lower=0> varhl;
  real<lower=0> varlh;
  real<lower=0> varll;
  real logmuhh;
  real logmuhl;
  real logmulh;
  real logmull;
  real logvarhh;
  real logvarhl;
  real logvarlh;
  real logvarll;
  real summandsHH[NHH];
  real summandsHL[NHL];
  real summandsLH[NLH];
  real summandsLL[NLL];

  rateh <- ratel + nu;
 
  muhh <- threshold1/rateh + threshold2/rateh;
  muhl <- threshold1/rateh + threshold2/ratel;
  mulh <- threshold1/ratel + threshold2/rateh;
  mull <- threshold1/ratel + threshold2/ratel;

  logmuhh <- log(muhh);
  logmuhl <- log(muhl);
  logmulh <- log(mulh);
  logmull <- log(mull);

  varhh <- threshold1/pow(rateh,3) + threshold2/pow(rateh,3);
  varhl <- threshold1/pow(rateh,3) + threshold2/pow(ratel,3);
  varlh <- threshold1/pow(ratel,3) + threshold2/pow(rateh,3);
  varll <- threshold1/pow(ratel,3) + threshold2/pow(ratel,3);

  logvarhh <- log(varhh);
  logvarhl <- log(varhl);
  logvarlh <- log(varlh);
  logvarll <- log(varll);

  for (n in 1:NHH) {
    summandsHH[n] <- 1.5*(logmuhh -log(yHH[n])) - .5 * ( log2() + logpi + logvarhh ) - muhh * square( yHH[n]-muhh) *inv(2*varhh*yHH[n]);
  }

  for (n in 1:NHL) {
    summandsHL[n] <- 1.5*(logmuhl -log(yHL[n])) - .5 * ( log2() + logpi + logvarhl ) - muhl * square( yHL[n]-muhl) *inv(2*varhl*yHL[n]);
  }

  for (n in 1:NLH) {
    summandsLH[n] <- 1.5*(logmulh -log(yLH[n])) - .5 * ( log2() + logpi + logvarlh ) - mulh * square( yLH[n]-mulh) *inv(2*varlh*yLH[n]);
  }

  for (n in 1:NLL) {
    summandsLL[n] <- 1.5*(logmull -log(yLL[n])) - .5 * ( log2() + logpi + logvarll ) - mull * square( yLL[n]-mull) *inv(2*varll*yLL[n]);
  }
}
model {
  threshold1 ~ gamma(4,.01);
  threshold2 ~ gamma(4,.01);
  ratel ~ gamma(4,.01);
  nu ~ exponential(100);
  increment_log_prob(sum(summandsHH) + sum(summandsHL) + sum(summandsLH) + sum(summandsLL));
}
