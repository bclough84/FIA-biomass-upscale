library(rstan)
library(dplyr)

#These data contain complete records >5 DBH, for almost all species.
#There were a few species that only had 1 or 2 records in the database and I removed these.
#This also contains both actual values and numerical indices (i.e., 1....N) for Jenkins species groups,
#Chojnacky species groups, and species for all observations
#The data have previously been converted to metric units (dbh (cm), mass (kg))
#These data contain complete records for 99 species. 
#The data have previously been converted to metric units (dbh (cm), mass (kg))
data <- read.csv("foliage_legacy_data.csv")
data <- subset(data, DW_FOL>0 & DBH>0)   #Removing a half dozen records with zeros, since log-transforming these will send Stan haywire. 
fia.data <- read.csv("MN_2009_2013_jenkins.csv")

fia.data <- fia.data[complete.cases(fia.data),] #Stan can't handle missing data



N <- nrow(data)   #Number of samples in legacy data
J <- 10 #Jenkins et al. species groups
y <- log(data$BIOMASS)   #Observed foliar biomass (kg)
fr <- log(data$DW_FOL/data$BIOMASS)
dbh <- log(data$DBH)      #Observed diameter at breast height (cm)  
grp <- data$JENKINS_SPGRPCD #1.....N index of Jenkins species groups

N.pred <- nrow(fia.data) #Number of prediction observations in FIA
fia.dbh <- log(fia.data$DIA) #DBH measurements from FIA
fia.grp <- fia.data$JENKINS_SPGRPCD #Jenkins group assignments for each tree


#Collect data into list object to pass to Stan
dat = list(y=y, dbh=dbh, fr=fr, grp=grp, N=N, J=J, fia_dbh=fia.dbh, fia_grp=fia.grp, N_new=N.pred)

#Specify the model
m.code="
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> N_new;
  int<lower=1, upper=J> grp[N];
  vector[N] y;
  vector[N] fr;
  vector[N] dbh;
  vector[N_new] fia_dbh;
  int<lower=1, upper=J> fia_grp[N_new];
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

   y ~ normal(y_hat, sigma);
   fr ~ normal(fr_hat, sigma_fr);
}
generated quantities{
   vector[N_new] y_pred; 
   vector[N_new] fr_pred;
   vector[N_new] fol_pred;            //Create new vector y_pred and predict
   for(n in 1:N_new)
      y_pred[n] <- normal_rng(beta0[fia_grp[n]] + beta1[fia_grp[n]]*fia_dbh[n], sigma);
   for(n in 1:N_new)
      fr_pred[n] <- normal_rng(fr_beta0[fia_grp[n]] + fr_beta1[fia_grp[n]]*fia_dbh[n], sigma_fr);
   for(n in 1:N_new)   
      fol_pred[n] <- exp(y_pred[n]) * exp(fr_pred[n]);
}
"
stan.fit <- stan(model_code=m.code, pars=c("fol_pred"), data=dat, chains=1, iter=1000, thin=2)   #The call to Stan to fit the model

fol_pred <- extract(stan.fit, "fol_pred")$fol_pred       #Extract the predicted values from Stan

fol_t <- t(fol_pred) #Transpose the matrix so nrow=N and ncol=# simulations

tpa <- fia.data$TPA_UNADJ

fol_adj <- fol_t*tpa #Return foliage biomass to original units and multiply by TPA adjustment factor

fol_adj_dat <- data.frame(fol_adj)

fol.ind <- cbind(fia.data$COUNTYCD, fia.data$PLOT, fol_adj_dat)

plot.level<- fol.ind %>%
                      group_by(fia.data$COUNTYCD, fia.data$PLOT) %>%
                      summarise_each(funs(sum), X1:X250)

plot.means <- rowMeans(plot.level[3:52])
plot.rng <- apply(plot.level[3:52], 1, function(x) { quantile(x, c(0.975)) - quantile(x, c(0.025)) } )

png("mean_uncertainty_hist_MN_plots.png", height=800, width=400)
par(mfrow = c(2,1))
hist(plot.means, xlab="Mean foliage biomass (kg/ha)",main=" ", ylab=" ", xlim=c(0, 25000), ylim=c(0,2500))
hist(plot.rng, xlab="95% uncertainty interval range (kg/ha)",main=" ", ylab=" ", xlim=c(0, 25000), ylim=c(0,2500))
