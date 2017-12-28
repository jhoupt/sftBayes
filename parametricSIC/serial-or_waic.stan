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
parameters {
  real<lower=0> ratel;
  real<lower=0> nu;
  real<lower=0> threshold1;
  real<lower=0> threshold2;
  real<lower=0,upper=1> order;
}
transformed parameters {
  real<lower=0> rateh;
  real logorder1;
  real logorder2;
  real logthresh1;
  real logthresh2;
  real logpi;
  real summandsHH[NHH];
  real summandsHL[NHL];
  real summandsLH[NLH];
  real summandsLL[NLL];

  rateh <- ratel + nu;

  logorder1 <- log(order);
  logorder2 <- log(1-order);
  logthresh1 <- log(threshold1);
  logthresh2 <- log(threshold2);
  logpi <- log(pi());

  for (n in 1:NHH) {
    summandsHH[n] <- log_sum_exp( logorder1 + logthresh1 - .5 * (log2() + logpi + 3*log(yHH[n])) - .5 * inv(yHH[n]) * square(rateh * yHH[n] - threshold1), 
				  logorder2 + logthresh2 - .5 * (log2() + logpi + 3*log(yHH[n])) - .5 * inv(yHH[n]) * square(rateh * yHH[n] - threshold2));
  }

  for (n in 1:NHL) {
    summandsHL[n] <- log_sum_exp( logorder1 + logthresh1 - .5 * (log2() + logpi + 3*log(yHL[n])) - .5 * inv(yHL[n]) * square(rateh * yHL[n] - threshold1), 
				  logorder2 + logthresh2 - .5 * (log2() + logpi + 3*log(yHL[n])) - .5 * inv(yHL[n]) * square(ratel * yHL[n] - threshold2));
  }

  for (n in 1:NLH) {
    summandsLH[n] <- log_sum_exp( logorder1 + logthresh1 - .5 * (log2() + logpi + 3*log(yLH[n])) - .5 * inv(yLH[n]) * square(ratel * yLH[n] - threshold1),
				  logorder2 + logthresh2 - .5 * (log2() + logpi + 3*log(yLH[n])) - .5 * inv(yLH[n]) * square(rateh * yLH[n] - threshold2));
  }

  for (n in 1:NLL) {
    summandsLL[n] <- log_sum_exp( logorder1 + logthresh1 - .5 * (log2() + logpi + 3*log(yLL[n])) - .5 * inv(yLL[n]) * square(ratel * yLL[n] - threshold1),
				  logorder2 + logthresh2 - .5 * (log2() + logpi + 3*log(yLL[n])) - .5 * inv(yLL[n]) * square(ratel * yLL[n] - threshold2));
  }
}
model {
  order ~ uniform(0,1);
  threshold1 ~ gamma(4,.01);
  threshold2 ~ gamma(4,.01);
  ratel ~ gamma(4,.01);
  nu ~ exponential(100);
  increment_log_prob( sum(summandsHH) + sum(summandsHL) + sum(summandsLH) + sum(summandsLL));
}
