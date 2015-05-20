library(rstan)
library(dplyr)

#These data contain complete records for 99 species. 
#The data have previously been converted to metric units (dbh (cm), mass (kg))
data <- read.csv("foliage_legacy_data.csv")
data <- subset(data, DW_FOL>0 & DBH>0)   #Removing a half dozen records with zeros, since log-transforming these will send Stan haywire. 


N <- nrow(data)   #Number of samples in legacy data
J <- 10 #Jenkins et al. species groups
agb <- log(data$BIOMASS)   #Observed foliar biomass (kg)
fr <- log(data$DW_FOL/data$BIOMASS)
dbh <- log(data$DBH)      #Observed diameter at breast height (cm)  
grp <- data$JENKINS_SPGRPCD #1.....N index of Jenkins species groups

#Collect data into list object to pass to Stan
dat = list(agb=agb, dbh=dbh, fr=fr, grp=grp, N=N, J=J)

#Specify the model
m.code="
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=1, upper=J> grp[N];
  vector[N] agb;
  vector[N] fr;
  vector[N] dbh;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_fr;
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  real<lower=0> fr_sigma_b0;
  real<lower=0> fr_sigma_b1;

  vector[J] beta0;
  vector[J] beta1;
  vector[J] fr_beta0;
  vector[J] fr_beta1;
  real theta0;
  real theta1;
  real fr_theta0;
  real fr_theta1;
}
transformed parameters{
   vector[N] y_hat; 
   vector[N] fr_hat;           

   for(i in 1:N)
     y_hat[i] <- beta0[grp[i]] + beta1[grp[i]]*dbh[i];
   for(i in 1:N)
     fr_hat[i] <- fr_beta0[grp[i]] + fr_beta1[grp[i]]*dbh[i];    
}
model{
   sigma ~ cauchy(0, 20);
   sigma_fr ~ cauchy(0, 20);
   sigma_b0 ~ cauchy(0, 20);
   sigma_b1 ~ cauchy(0, 20);
   fr_sigma_b0 ~ cauchy(0, 20);
   fr_sigma_b1 ~ cauchy(0, 20);

   theta0 ~ normal(0, 25);
   theta1 ~ normal(0, 25);
   fr_theta0 ~ normal(0, 25);
   fr_theta1 ~ normal(0, 25);

   for(j in 1:J)
         beta0[j] ~ normal(theta0, sigma_b0);
   for(j in 1:J)
         beta1[j] ~ normal(theta1, sigma_b1);
   for(j in 1:J)
         fr_beta0[j] ~ normal(fr_theta0, fr_sigma_b0);
   for(j in 1:J)
         fr_beta1[j] ~ normal(fr_theta1, fr_sigma_b1);

   agb ~ normal(y_hat, sigma);
   fr ~ normal(fr_hat, sigma_fr);
}
generated quantities{
   vector[N] exp_agb;
   vector[N] exp_fr;
   vector[N] folb;
   for(i in 1:N)
       exp_agb[i] <- exp(agb[i]);
   for(i in 1:N)
       exp_fr[i] <- exp(fr[i]);
   for(i in 1:N)
       folb[i] <- exp_agb[i] * exp_fr[i];
}
"
#stan.fit <- stan(model_code=m.code, pars=c("beta0", "beta1", "fr_beta0", "fr_beta1", "sigma", "sigma_fr"), data=dat, iter=1000)   #The call to Stan to fit the model

foliage <- stan(model_code=m.code, pars=c("folb"), data=dat, iter=500) 

fol_pred <- extract(foliage, "folb")$folb       #Extract the predicted values from Stan
