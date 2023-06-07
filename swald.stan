functions{
    // shifted Wald (log) PDF
    real shiftedWald_lpdf(real X, real gamma, real alpha, real theta){
        real tmp1;
        real tmp2;
        tmp1 = alpha / (sqrt(2 * pi() * (pow((X - theta), 3))));
        tmp2 = exp(-1 * (pow((alpha - gamma * (X-theta)),2)/(2*(X-theta))));
        return log(tmp1*tmp2) ;
    }
    // shifted Wald (log) survivor function
    real shiftedWald_lS(real X, real gamma, real alpha, real theta){
      real tmp1;
      real tmp2;
      real tmp3;
      // compute log of the survival function 1-F for swald
      tmp1 = Phi((gamma*(X-theta) - alpha)/sqrt(X-theta));
      tmp2 = exp(2*alpha*gamma);
      tmp3 = Phi(-1*(gamma*(X-theta) + alpha)/sqrt(X-theta));
      return log(1-tmp1-tmp2*tmp3);
    }
}

data{
    int ns; // number of subjects
    int nrt; // number of trials per subject 
    real rt[ns,nrt]; // matrix of response times
    real D[ns,nrt]; // censoring matrix (0 = correct RT, 1 = incorrect RT)
}

parameters{
    // noncentered parameters
    vector[ns] G_raw;
    vector[ns] A_raw;
    vector[ns] H_raw;
    // define uniform hyperpriors from Matzke & Wagenmakers (2009)
    real<lower=0.85,upper=7.43> g; 
    real<lower=0.67,upper=2.35> a;
    real<lower=0,upper=0.82> h;
    real<lower=0,upper=1.899> gS;
    real<lower=0,upper=0.485> aS;
    real<lower=0,upper=0.237> hS;
 }

transformed parameters{
    // noncentered parameterization (Matt trick) (Betancourt & Girolami, 2003)
    vector[ns] G;
    vector[ns] A;
    vector[ns] H;
    
    for(i in 1:ns){
      G[i] = g + gS*G_raw[i];
      A[i] = a + aS*A_raw[i];
      H[i] = h + hS*H_raw[i];
    }
}

model{
    for(i in 1:ns){
      // group level prior
      G_raw[i] ~ normal(0,1);
      A_raw[i] ~ normal(0,1);
      H_raw[i] ~ normal(0,1);  
    
      for(j in 1:nrt){
        target += (1-D[i,j])*shiftedWald_lpdf(rt[i,j]| G[i], A[i], H[i]) + D[i,j]*shiftedWald_lS(rt[i,j], G[i], A[i], H[i]);
      }
  }
}
