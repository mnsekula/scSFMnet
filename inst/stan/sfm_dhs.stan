data {                          // Horseshoe prior on alpha0 and alpha1, tau varying with factors
  int<lower=0> G;               // number of genes
  int<lower=0> N;               // number of cells
  int<lower=0> Fac;             // number of factors
  int<lower=0> Y[G,N];          // expression matrix
  int<lower=0> t_i[N];          // indicator vector of treatment/control
}// end data

parameters {
  matrix[Fac,N] lambda;     // Poisson distribution with treatment groups   
  matrix[G,Fac] alpha_0;        // alphas for control: t_i=0   
  matrix[G,Fac] alpha_1;        // alphas for treatment: t_i=1 
  matrix[G,Fac] alpha_star;     // mean of alpha_0 and alpha_1
  real<lower=0> zeta;    // horseshoe global of alpha_*
  matrix<lower=0>[G,Fac] omega; // horseshoe local of alpha_*
  vector[G] beta;
  vector[G] delta;              //mean adjustment term for treatment t_i=1
  real<lower=0> sigma_b;
  real<lower=0> sigma_d;
  
  vector<lower=0>[Fac] tau; //horseshoe global
  matrix<lower=0>[G,Fac] eta0;
  matrix<lower=0>[G,Fac] eta1;
}// end parameters

transformed parameters {
  matrix[G,Fac] sq_alpha_0;
  matrix[G,Fac] sq_alpha_1;
  
  for(g in 1:G){
    for(f in 1:Fac){
      sq_alpha_0[g,f] = square(alpha_0[g,f]);
      sq_alpha_1[g,f] = square(alpha_1[g,f]);
    }
  }
  
} // end tp

model {
  matrix[G,N] a0_l;
  matrix[G,N] a1_l;
  vector[G] l_adjust0;
  vector[G] l_adjust1;
  a0_l = alpha_0*lambda;
  a1_l = alpha_1*lambda;

  for(g in 1:G){
    l_adjust0[g] = sum(sq_alpha_0[g,])/2;
    l_adjust1[g] = sum(sq_alpha_1[g,])/2;
    for(i in 1:N){
      if(t_i[i] == 0){
        //control group
        Y[g,i] ~ poisson_log(beta[g]+a0_l[g,i]-l_adjust0[g]);
      }else{
        //treatment group
        Y[g,i] ~ poisson_log(delta[g]+beta[g]+a1_l[g,i]-l_adjust1[g]);
      }//end if/else 
    }// end cell loop
  }// end gene loop
    
  
  for(f in 1:Fac){
    //horseshoe prior, factor dependent
    for(g in 1:G){
      alpha_star[g,f] ~ normal(0,omega[g,f]*zeta);
      alpha_0[g,f] ~ normal(alpha_star[g,f],eta0[g,f]*tau[f]);
      alpha_1[g,f] ~ normal(alpha_star[g,f],eta1[g,f]*tau[f]);
    }
  }
  to_vector(lambda) ~ normal(0,1);
  
  zeta ~ cauchy(0,1);
  to_vector(omega) ~ cauchy(0,1);

  beta ~ normal(0,sigma_b);
  delta ~ normal(0,sigma_d);
  sigma_b ~ cauchy(0,1);
  sigma_d ~ cauchy(0,1);
  
// Horseshoe prior?
  tau ~ cauchy(0,1);
  to_vector(eta0) ~ cauchy(0,1);
  to_vector(eta1) ~ cauchy(0,1);

  
}// end model

generated quantities {
  matrix[G,G] cov_mu_mu0;
  matrix[G,G] corr_mu_mu0;
  matrix[G,G] cov_mu_mu1;
  matrix[G,G] corr_mu_mu1;
  
  for(g in 1:G){
    for(h in 1:G){
      cov_mu_mu0[g,h] = dot_product(alpha_0[g,],alpha_0[h,]);
      corr_mu_mu0[g,h] = cov_mu_mu0[g,h]/sqrt(sum(sq_alpha_0[g,])*sum(sq_alpha_0[h,]));
      cov_mu_mu1[g,h] = dot_product(alpha_1[g,],alpha_1[h,]);
      corr_mu_mu1[g,h] = cov_mu_mu1[g,h]/sqrt(sum(sq_alpha_1[g,])*sum(sq_alpha_1[h,]));
    }
  }
  
  
}// end gq

