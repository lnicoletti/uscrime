library(rjags)

#let's start by setting our directory and importing relevant packages
setwd("C:/Users/Leonardo/Documents/TU_Delft/EPA_Year1/EPA1315/FinalProject")

#load dataset
SF_LA_PH_DC_prcnts_vars <- read.csv('Data/AllCities/Tables/SF_LA_PH_DC_BO_prcnts_vars.csv')
#panel regression with fixed effects. Do not normalize variables, 
#instead we set all variables on logarithmic scale
crime <- log1p(SF_LA_PH_DC_prcnts_vars$crimecount)
trees <- log1p(SF_LA_PH_DC_prcnts_vars$treecount)
lights <- log1p(SF_LA_PH_DC_prcnts_vars$lightcount)
trash <- log1p(SF_LA_PH_DC_prcnts_vars$trashcount)
SF <- log1p(SF_LA_PH_DC_prcnts_vars$CITY_San.Francisco)
LA <- log1p(SF_LA_PH_DC_prcnts_vars$CITY_Los.Angeles)
PH <- log1p(SF_LA_PH_DC_prcnts_vars$CITY_Philadelphia)
DC <- log1p(SF_LA_PH_DC_prcnts_vars$CITY_Washington.DC)
BO <- log1p(SF_LA_PH_DC_prcnts_vars$CITY_Boston)
n <- nrow(SF_LA_PH_DC_prcnts_vars)

model_string1 <- "model{

#likelihood
for (i in 1:n){
  crime[i] ~ dnorm(mu[i], inv.var) #insert function that represents distribution of crime 
  mu[i] <- beta0 + beta1*trees[i] + beta2*lights[i] + beta3*trash[i] + beta4*SF[i] + beta5*LA[i] + beta6*PH[i] + beta7*DC[i]
}

#one prior for each beta
beta0 ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)
beta3 ~ dnorm(0, 0.001)
beta4 ~ dnorm(0, 0.001)
beta5 ~ dnorm(0, 0.001)
beta6 ~ dnorm(0, 0.001)
beta7 ~ dnorm(0, 0.001)


#prior for the inverse variance
inv.var ~ dgamma(0.01, 0.0001)
sigma <- 1/sqrt(inv.var)
}" 

#initialize model
model <- jags.model(textConnection(model_string1), 
                    data = list(crime=crime,trees=trees,
                                lights=lights,trash=trash,SF=SF,
                                LA=LA,PH=PH,DC=DC, n=n), n.chains = 2)

update(model, 10000, progress.bar="none"); # Burnin for 10000 samples
samp <- coda.samples(model,variable.names=c("beta0","beta1","beta2","beta3","beta4", 
                                            "beta5","beta6", "beta7", "sigma"), 
                                            n.iter=20000, progress.bar="text")
#get summary of posterior samples
summary(samp)

traceplot(samp)
# sometimes the gelman plot won't fit on a screen
# we have to reduce the margins
par(mar=c(.4,.4,.4,.4))
gelman.plot(samp)
gelman.diag(samp)
densplot(samp)
acfplot(samp)

#visualization of the complete posterior distributions 
#with traceplots for each parameter, we first convert the model object to class mcmc
samp_mcmc <- as.mcmc(samp)
plot(samp)

# get the effective sample size
effectiveSize(samp)

# take the mean of each predicted coefficient
coeff_preds <- as.matrix(samp)
coeff_preds <- as.data.frame(coeff_preds)
beta0_pred <- mean(coeff_preds$beta0)
beta1_pred <- mean(coeff_preds$beta1)
beta2_pred <- mean(coeff_preds$beta2)
beta3_pred <- mean(coeff_preds$beta3)
beta4_pred <- mean(coeff_preds$beta4)
beta5_pred <- mean(coeff_preds$beta5)
beta6_pred <- mean(coeff_preds$beta6)
beta7_pred <- mean(coeff_preds$beta7)
sigma_pred <- mean(coeff_preds$sigma)
#compute predicted crime with model betas
crime_predictor <- function(var1, var2, var3, var4, var5, var6, var7, var8){
  # vars: trees, lights, trash, SF, LA, PH, DC, n
  #empty dataframe to be filled with predicted crime values
  crime_pred <- data.frame(matrix(ncol = 1, nrow = n))
  x <- c("predicted_crime") 
  colnames(crime_pred) <- x
    #function for crime
  crime_pred$predicted_crime <- beta0_pred + beta1_pred*var1 + beta2_pred*var2 + 
    beta3_pred*var3 + beta4_pred*var4 + beta5_pred*var5 + 
    beta6_pred*var6 + beta7_pred*var7
  return(crime_pred)
}

#for (i in 1:n){
#crime_pred$predicted_crime[i] <- beta0_pred + beta1_pred*var1[i] + beta2_pred*var2[i] + 
#  beta3_pred*var3[i] + beta4_pred*var4[i] + beta5_pred*var5[i] + 
#  beta6_pred*var6[i] + beta7_pred*var7[i]

crime_predictions <- crime_predictor(trees, lights, trash, SF, LA, PH, DC, n)
crime_obs_pred <- cbind(crime_predictions, crime)

library('ggplot2')
# plot predicted crime versus observed crime
ggplot(crime_obs_pred, aes(x=crime, y=predicted_crime)) +
  geom_point() +
  geom_smooth(method=lm)

plot(crime_obs_pred$crime, crime_obs_pred$predicted_crime, ylab='Predicted Crime', xlab='Observed Crime')
lines(seq(1,3500),seq(1,3500), col='red', lwd=2, cex=3)


model_output<-as.matrix(samp)
saveRDS(samp,"ModelResults/modelrun1.RDS")
write.csv(model_output,"modelruns2.csv",row.names = FALSE)

