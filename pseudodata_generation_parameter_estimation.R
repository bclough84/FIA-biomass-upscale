library(rstan)
library(truncnorm)

data <- read.csv("MN_2009_2013_jenkins.csv")
data$DIA <- data$DIA*2.54

#############################################
#STEP ONE: GENERATE EMPIRICAL BIOMASS ESTIMATES FOR PSEUDOPOPULATIONS
#############################################

#Divide tree table by Jenkins species group, fit a Weibull distribution to each,
#then use the fitted parameters to generate theoretical diameter distributions for
#each species group in the sample. NOTE: Species groups 2 & 10 not represented in Minnesota
grp1 <- subset(data, JENKINS_SPGRPCD==1)
fit1 <- fitdistr(grp1$DIA, "weibull")
grp3 <- subset(data, JENKINS_SPGRPCD==3)
fit3 <- fitdistr(grp3$DIA, "weibull")
grp4 <- subset(data, JENKINS_SPGRPCD==4)
fit4 <- fitdistr(grp4$DIA, "weibull")
grp5 <- subset(data, JENKINS_SPGRPCD==5)
fit5 <- fitdistr(grp5$DIA, "weibull")
grp6 <- subset(data, JENKINS_SPGRPCD==6)
fit6 <- fitdistr(grp6$DIA, "weibull")
grp7 <- subset(data, JENKINS_SPGRPCD==7)
fit7 <- fitdistr(grp7$DIA, "weibull")
grp8 <- subset(data, JENKINS_SPGRPCD==8)
fit8 <- fitdistr(grp8$DIA, "weibull")
grp9 <- subset(data, JENKINS_SPGRPCD==9)
fit9 <- fitdistr(grp9$DIA, "weibull")

sim1 <- rweibull(1000, shape=fit1$estimate[1], scale=fit1$estimate[2])
sim3 <- rweibull(1000, shape=fit3$estimate[1], scale=fit3$estimate[2])
sim4 <- rweibull(1000, shape=fit4$estimate[1], scale=fit4$estimate[2])
sim5 <- rweibull(1000, shape=fit5$estimate[1], scale=fit5$estimate[2])
sim6 <- rweibull(1000, shape=fit6$estimate[1], scale=fit6$estimate[2])
sim7 <- rweibull(1000, shape=fit7$estimate[1], scale=fit7$estimate[2])
sim8 <- rweibull(1000, shape=fit8$estimate[1], scale=fit8$estimate[2])
sim9 <- rweibull(1000, shape=fit9$estimate[1], scale=fit9$estimate[2])

#Create some new vectors to work with
sim.dbh <- c(sim1, sim3, sim4, sim5, sim6, sim7, sim8, sim9)
jenkins.ind <- c(rep(1, 1000),rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(6,1000), rep(7,1000), rep(8,1000), rep(9,1000))

#This part of the code generates empirical estimates of total biomass and foliar ratio
#based on the parameter estimates from Jenkins et al. 2003
N <- length(sim.dbh)
bm.tot <- rep(NA, N)
fol.rat <- rep(NA, N)

for(i in 1:N){
if(jenkins.ind[i]==1){
       bm.tot[i] <- exp(-2.0336 + 2.2592*log(sim.dbh[i]))
       fol.rat[i] <- exp(-2.9584 + (4.4766/sim.dbh[i]))
   } else if(jenkins.ind[i]==3){
          bm.tot[i] <- exp(-2.5384 + 2.4814*log(sim.dbh[i]))
          fol.rat[i] <- exp(-2.9584 + (4.4766/sim.dbh[i]))
   } else if(jenkins.ind[i]==4){
             bm.tot[i] <- exp(-2.5356 + 2.4449*log(sim.dbh[i]))
             fol.rat[i] <- exp(-2.9584 + (4.4766/sim.dbh[i]))
   } else if(jenkins.ind[i]==5){
             bm.tot[i] <- exp(-2.0773 + 2.3323*log(sim.dbh[i]))
             fol.rat[i] <- exp(-2.9584 + (4.4766/sim.dbh[i]))
   } else if(jenkins.ind[i]==6){
             bm.tot[i] <- exp(-2.2094 + 2.3867*log(sim.dbh[i]))
             fol.rat[i] <- exp(-4.0813 + (5.8816/sim.dbh[i]))
   } else if(jenkins.ind[i]==7){
             bm.tot[i] <- exp(-1.9123 + 2.3651*log(sim.dbh[i]))
             fol.rat[i] <- exp(-4.0813 + (5.8816/sim.dbh[i]))
   } else if(jenkins.ind[i]==8){
             bm.tot[i] <- exp(-2.4800 + 2.4835*log(sim.dbh[i]))
             fol.rat[i] <- exp(-4.0813 + (5.8816/sim.dbh[i]))
   } else if(jenkins.ind[i]==9){
             bm.tot[i] <- exp(-2.0127 + 2.4342*log(sim.dbh[i]))
             fol.rat[i] <- exp(-4.0813 + (5.8816/sim.dbh[i]))
   }
}

fol.bm <- bm.tot * fol.rat
pseudodat <- as.data.frame(cbind(jenkins.ind, sim.dbh, fol.bm, bm.tot, fol.rat))
#Write this file in case we need it again later
write.csv(pseudodat, "simulated_biomass_data_jenkins.csv")

#########################################################
#STEP TWO: FUZZ THE BIOMASS DATA
#########################################################

stdev <- 15 #Based on fit of Jenkins model to the legacy data
fol.bm.fuzz <- rtruncnorm(n=nrow(pseudodat), a=0, b=Inf, mean=pseudodat$fol.bm, sd=stdev)  #Truncated normal to avoid negative foliage biomass estimates

fuzz.data <- as.data.frame(cbind(jenkins.ind, sim.dbh, fol.bm.fuzz))

write.csv(fuzz.data, "pseudodata_fuzzed.csv")

##########################################################
#STEP THREE: FIT THE MODEL TO THE FUZZED DATA TO GET PARAMETER ESTIMATES
##########################################################

jenkins_ind <- fuzz.data$jenkins.ind
sim_dbh <- log(fuzz.data$sim.dbh)
fol_bm <- log(fuzz.data$fol.bm.fuzz)

N <- nrow(fuzz.data)
J <- 10

dat <- list(jenkins_ind=jenkins_ind, sim_dbh=sim_dbh, fol_bm=fol_bm, N=N, J=J)

m.code="
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=1, upper=J> jenkins_ind[N];
  vector[N] fol_bm;
  vector[N] sim_dbh;

}
parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  vector[J] beta0;
  vector[J] beta1;
  real mu_b0;
  real mu_b1;
}
transformed parameters{
   vector[N] y_hat;

   for(i in 1:N)
     y_hat[i] <- beta0[jenkins_ind[i]] + beta1[jenkins_ind[i]]*sim_dbh[i];
}
model{
   sigma ~ cauchy(0, 10);      //weakly informative priors on the mean/variance of a & b
   sigma_b0 ~ cauchy(0, 10);
   sigma_b1 ~ cauchy(0, 10);

   mu_b0 ~ normal(0, 25);
   mu_b1 ~ normal(0, 25);

   for(i in 1:10)
       beta0[i] ~ normal(mu_b0, sigma_b0);
   for(i in 1:10)
       beta1[i] ~ normal(mu_b1, sigma_b1);

   fol_bm ~ normal(y_hat, sigma);
}
"

##Call to Stan to evaluate model parameters
stan.fit <- stan(model_code=m.code, pars=c("beta0", "beta1"), data=dat,iter=1000, chains=2)


