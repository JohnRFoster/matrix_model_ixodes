###-------------------------------------------------------###
### Survival HB sub model for Cary ticks                  ### 
### Estimates daily survival for each site                ###
###                                                       ###
### Subset: Flat Nymphs                                   ###
### Drivers: Max Temp                                     ###
###-------------------------------------------------------###

library(ecoforecastR)
source("Functions/surv_type_iButton.R")
source("Functions/RunMCMC_Model.R")

## file path to output folder
out.folder <- "../FinalOut/SurvivalModels/Flat_Larva/"

## set the subset of ticks to estimate survival
subset <- "Flat Nymph"
subset.out <- gsub(pattern = " ", replacement = "_", x = subset)

## set season (only need for larvae)
# season <- "winter"

# set driver
driver.vec <- c("temp", "rh", "vpd")
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
driver <- driver.vec[xx]
position <- "Above"
out.name <- paste(subset.out,driver,position, "AllTicks_NoRandomNoDiff",sep = "_")


out.path <- paste(out.folder,out.name, sep = "")

# read data
data <- surv_type_iButton(subset, position)
data$b_mu <- c(0, 0)
data$b_prec <- diag(0.001, 2, 2)

if(driver == "temp"){
  data$met <- data$temp
} else if(driver == "rh"){
  data$met <- data$rh
} else {
  data$met <- data$vpd
} 
data[c("temp", "rh", "vpd")] <- NULL

cat("Number of cores:", nrow(data$y), "\n")

survival.model <- " model {

## priors
beta ~ dmnorm(b_mu, b_prec)

## loop for daily survival for each row
for(i in 1:N){
  for(t in 2:(N_days[i])){
    logit(lambda[i,t]) <- beta[1] + beta[2]*met[t,i] 
  }
}


## loop for aggrigated survival across time in each row
for(t in 1:N){
  lambda[t,1] ~ dbeta(20,1)
  log(phi[t]) <- sum(log(lambda[t,1:(N_days[t]-1)]))
}

## Binomial process for ticks surviving each day
for(i in 1:N){
  y[i,2] ~ dbinom(phi[i], y[i,1])
}

}" # end model

inits <- function(){list(beta = rnorm(2,c(3.5,0),0.001),
                         tau = rgamma(1,0.01,0.01))}

monitor <- "beta"

load.module("glm")
j.model <- jags.model(file = textConnection(survival.model),
                      data = data,
                      inits = inits,
                      n.adapt = 7500,
                      n.chains = 3)

jags.out <- runMCMC_Model(j.model = j.model, variableNames = monitor)

save(jags.out, j.model, file = paste(out.path, ".RData", sep = ""))


print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")


