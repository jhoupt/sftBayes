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
  real logthresh1;
  real logthresh2;
  real summandsHH[NHH];
  real summandsHL[NHL];
  real summandsLH[NLH];
  real summandsLL[NLL];

  logthresh1 <- log(threshold1);
  logthresh2 <- log(threshold2);
  rateh <- ratel + nu;

  for (n in 1:NHH) {
    summandsHH[n] <- log_sum_exp( 
        logthresh1 - .5 * (log2() + logpi +3*log(yHH[n])) - square(rateh*yHH[n] - threshold1)/(2*yHH[n]) + 
	log(Phi( (rateh*yHH[n]-threshold2)/sqrt(yHH[n]) ) + exp(2*threshold2*rateh)* Phi(-(rateh*yHH[n]+threshold2)/sqrt(yHH[n]) ) ),
        logthresh2 - .5 * (log2() + logpi +3*log(yHH[n])) - square(rateh*yHH[n] - threshold2)/(2*yHH[n]) + 
	log(Phi( (rateh*yHH[n]-threshold1)/sqrt(yHH[n]) ) + exp(2*threshold1*rateh)* Phi(-(rateh*yHH[n]+threshold1)/sqrt(yHH[n]) ) ) );
  }
  for (n in 1:NHL) {
    summandsHL[n] <- log_sum_exp( 
        logthresh1 - .5 * (log2() + logpi +3*log(yHL[n])) - square(rateh*yHL[n] - threshold1)/(2*yHL[n]) + 
	log(Phi( (ratel*yHL[n]-threshold2)/sqrt(yHL[n]) ) + exp(2*threshold2*ratel)* Phi(-(ratel*yHL[n]+threshold2)/sqrt(yHL[n]) ) ),
        logthresh2 - .5 * (log2() + logpi +3*log(yHL[n])) - square(ratel*yHL[n] - threshold2)/(2*yHL[n]) + 
	log(Phi( (rateh*yHL[n]-threshold1)/sqrt(yHL[n]) ) + exp(2*threshold1*rateh)* Phi(-(rateh*yHL[n]+threshold1)/sqrt(yHL[n]) ) ) );
  }
  for (n in 1:NLH) {
    summandsLH[n] <- log_sum_exp( 
        logthresh1 - .5 * (log2() + logpi +3*log(yLH[n])) - square(ratel*yLH[n] - threshold1)/(2*yLH[n]) + 
	log(Phi( (rateh*yLH[n]-threshold2)/sqrt(yLH[n]) ) + exp(2*threshold2*rateh)* Phi(-(rateh*yLH[n]+threshold2)/sqrt(yLH[n]) ) ),
        logthresh2 - .5 * (log2() + logpi +3*log(yLH[n])) - square(rateh*yLH[n] - threshold2)/(2*yLH[n]) + 
	log(Phi( (ratel*yLH[n]-threshold1)/sqrt(yLH[n]) ) + exp(2*threshold1*ratel)* Phi(-(ratel*yLH[n]+threshold1)/sqrt(yLH[n]) ) ) );
  }
  for (n in 1:NLL) {
    summandsLL[n] <- log_sum_exp( 
        logthresh1 - .5 * (log2() + logpi +3*log(yLL[n])) - square(ratel*yLL[n] - threshold1)/(2*yLL[n]) + 
	log(Phi( (ratel*yLL[n]-threshold2)/sqrt(yLL[n]) ) + exp(2*threshold2*ratel)* Phi(-(ratel*yLL[n]+threshold2)/sqrt(yLL[n]) ) ),
        logthresh2 - .5 * (log2() + logpi +3*log(yLL[n])) - square(ratel*yLL[n] - threshold2)/(2*yLL[n]) + 
	log(Phi( (ratel*yLL[n]-threshold1)/sqrt(yLL[n]) ) + exp(2*threshold1*ratel)* Phi(-(ratel*yLL[n]+threshold1)/sqrt(yLL[n]) ) ) );
  }
}
model {
  threshold1 ~ gamma(4,.01);
  threshold2 ~ gamma(4,.01);
  ratel ~ gamma(4,.01);
  nu ~ exponential(100);
  increment_log_prob( sum(summandsHH) + sum(summandsHL) + sum(summandsLH) + sum(summandsLL) );
}
