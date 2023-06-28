library(tidyverse)
library(LaplacesDemon)
library(nimble)
library(rjags)
library(scoringRules)
library(lubridate)

source("Scripts/modelOutPathsDormantStages.R")

analysis.dir <- "FinalOut/A_Correct/Analysis"
save.dir <- file.path(analysis.dir, "dormantStages/nimbleModels/oneStepPredictions")

sites <- c("Green", "Henry", "Tea")
ua.exp <- c(
  "ic",
  "parameter",
  "process",
  "ic_parameter",
  "ic_process",
  "parameter_process",
  "ic_parameter_process",
  "mice",
  "mice_ic",
  "mice_parameter",
  "mice_process",
  "mice_ic_parameter",
  "mice_ic_process",
  "mice_parameter_process",
  "mice_ic_parameter_process"
)

models <- names(dormant.models)
models <- models[models != "DormantNymphWithWeatherAndMiceSiteSpecific"]
n.models <- length(models)
n.sites <- length(sites)
n.ua <- length(ua.exp)

total.jobs <- n.sites^2 * n.models * n.ua

all.jobs <- expand_grid(
  model = models,
  paramsFrom = rep(sites, each = n.sites),
  ticksFrom = rep(sites, n.sites),
  ua = ua.exp
) %>%
  distinct()

nrow(all.jobs) == total.jobs

# not all combinations are necessary (mice in ua but not in model)
all.jobs <- all.jobs %>%
  filter(!(grepl("mice", ua) & !grepl("Mice", model)))

# all.jobs <- all.jobs %>%
#   filter(model == "DormantNymphWithMNAMice",
#          !grepl("mice", ua))

array.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# if(is.na(array.num)) array.num <- sample.int(nrow(all.jobs), 1)
# array.num=63
Nmc <- 5000

job <- all.jobs %>% slice(array.num)
model <- job %>% pull(model)
paramsFrom <- job %>% pull(paramsFrom)
ticksFrom <- job %>% pull(ticksFrom)
ua <- job %>% pull(ua)


met.site <- "global"

message(paste("Model:", model))
message(paste("Parameters From:", paramsFrom))
message(paste("Ticks From:", ticksFrom))
message(paste("UA:", ua))

## get model specifics
model.dir <- pluck(dormant.models, model, "model.dir") # where model output is stored
process.met <- pluck(dormant.models, model, "process.driver") # is weather in the model
use.mice <- pluck(dormant.models, model, "mice") # are mice in the model
mice.ua <- grepl("mice", ua) # is mouse uncertainty in the simulation
mice.site <- NA
mna <- FALSE
if(use.mice){ # if we are using mice, which site do they come from
  mice.site <- ticksFrom # the same site as ticks
  mna <- pluck(dormant.models, model, "mna")
}

if(!use.mice & mice.ua) stop("Mice in UA when mice are not in model")

# function to load posterior samples
load_samples <- function(samp.dir, params, same.site = FALSE, draws = NA){
  out.files <- list.files(samp.dir)
  if("NIMBLE_samples.RData" %in% out.files){
    fileName <- "NIMBLE_samples.RData"
    load(file.path(samp.dir, fileName))
    samples <- save.ls$samples$samples
    if(params){
      mat <- samples[,-grep("x[", colnames(samples), fixed = TRUE)]  
    } else {
      mat <- samples[,grep("x[", colnames(samples), fixed = TRUE)]  
    }
  } else if("NIMBLE_check.RData" %in% out.files){
    fileName <- "NIMBLE_check.RData"
    load(file.path(samp.dir, fileName))
    if(params){
      mat <- as.matrix(save.ls$params)  
    } else {
      mat <- as.matrix(save.ls$predict)
    }
  } else {
    stop(paste("No files in", samp.dir))
  }
  
  if(!same.site){
    draws <- sample.int(nrow(mat), Nmc, replace = TRUE)  
  }
  
  post <- mat %>% 
    as_tibble() %>% 
    slice(draws)
  
  return(list(post = post,
              draws = draws))
}

param.dir <- file.path(top.dir, model.dir, paramsFrom)
tick.dir <- file.path(top.dir, model.dir, ticksFrom)

## load parameters
message("Loading Parameters")
post.params <- load_samples(param.dir, TRUE) 
params <- post.params$post

# want the same draws of params and ticks are from the same site
same.site <- if_else(paramsFrom == ticksFrom, TRUE, FALSE)

## load ticks
message("Loading Ticks")
post.ticks <- load_samples(tick.dir, FALSE, same.site, post.params$draws)
ticks <- post.ticks$post

## separate model parameters and process error
process.error <- params %>% select(starts_with("sig"))
model.params <- params %>% select(!starts_with("sig"))

## get dates
tick.raw <- read_csv("../Data/Cary_ticks.csv")
tick.actual <- tick.raw %>% 
  filter(Grid == paste(ticksFrom, "Control")) %>% 
  select(DATE, n_larvae, n_nymphs, n_adults)

source("Functions/cary_tick_met_JAGS.R")
data <- cary_ticks_met_JAGS()

source("Functions/site_data_met.R")
data <- site_data_met(
  site = paste(ticksFrom, "Control"),
  met.variable = process.met,
  data = data
)

d.index <- 2
q.index <- c(1,3,4)
y <- matrix(0, 4, data$N_est)
y[1,] <- data$y[1,]
y[2,] <- rep(NA, data$N_est)
y[3,] <- data$y[2,]
y[4,] <- data$y[3,]
data$y <- y
data$ns <- nrow(y)

## manipulate for UA experiment
if(!grepl("ic", ua)){
  ticks.mu <- apply(ticks, 2, mean)
  for(c in 1:ncol(ticks)){
    ticks[,c] <- ticks.mu[c]
  }
}

if(!grepl("parameter", ua)){
  dem.params.sub <- apply(model.params, 2, mean)
  for(c in 1:ncol(model.params)){
    model.params[,c] <- dem.params.sub[c]
  }
}

if(!grepl("process", ua)){
  process.error[process.error != 0] <- 0
}

# combine all parameters
params <- bind_cols(model.params, process.error)

if(process.met){
  source("Functions/nimble_functions.R")
  nim.data <- nimble_data(paste(ticksFrom, "Control"))
  met <- nim.data$met
  
  X <- met %>% 
    select(MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC)
  
} else {
  X <- met.site <- NA
}

message("Running simulation...")
source("predict_one.R")
out <- predict_one(
  params = params,
  ticks = ticks,
  driver = X,
  data = data,
  mice.site = mice.site,
  mna = mna,
  mice.ua = mice.ua,
  met.site = met.site
)

# for time series we just need summary stats
get_stats <- function(out){
  all.stats <- tibble()
  life.stages <- c("larvae", "dormant", "nymphs", "adults")
  for(ls in 1:4){
    for(u in 1:2){
      
      observation <- tick.actual %>% 
        mutate(DATE = as.character(DATE),
               n_dormant = NA) %>% 
        select(DATE, contains(life.stages[ls]))
      
      df.ls <- t(out[ls,,u,]) 
      colnames(df.ls) <- as.character(tick.actual$DATE)
      
      df.tb <- df.ls %>% 
        as_tibble()
      
      if(u == 1){
        type <- "latent"
      } else {
        type <- "observation"
        df.tb <- df.tb %>% select(-1)
      }
      
      df.tb <- df.tb %>% 
        pivot_longer(cols = everything(),
                     names_to = "DATE",
                     values_to = "value") %>% 
        group_by(DATE) %>% 
        summarise(mean = mean(value),
                  var = var(value),
                  conf_2.5 = quantile(value, 0.025),
                  conf_5 = quantile(value, 0.05),
                  conf_10 = quantile(value, 0.10),
                  conf_20 = quantile(value, 0.20),
                  conf_25 = quantile(value, 0.25),
                  conf_30 = quantile(value, 0.30),
                  conf_40 = quantile(value, 0.40),
                  conf_50 = quantile(value, 0.50),
                  conf_60 = quantile(value, 0.60),
                  conf_70 = quantile(value, 0.70),
                  conf_75 = quantile(value, 0.75),
                  conf_80 = quantile(value, 0.80),
                  conf_90 = quantile(value, 0.90),
                  conf_95 = quantile(value, 0.95),
                  conf_97.5 = quantile(value, 0.975)) %>% 
        pivot_longer(cols = -DATE,
                     names_to = "statistic",
                     values_to = "value") %>% 
        mutate(lifeStage = life.stages[ls],
               type = type,
               model = model,
               paramsFrom = paramsFrom,
               ticksFrom = ticksFrom,
               ua = ua)
      
      with.obs <- left_join(df.tb, observation, by = "DATE") 
      with.obs.name <- with.obs %>% 
        rename("obs" = starts_with("n_"))
      
      all.stats <- bind_rows(all.stats, with.obs.name)
      
    }
  }
  return(all.stats)
}

out.stats <- get_stats(out)
out.name <- paste("timeSeriesQuantiles", model.dir, "paramsFrom", paramsFrom, "ticksFrom", ticksFrom, ua, sep = "_")
csv.name <- paste0(out.name, ".csv")
outFile <- file.path(save.dir, csv.name)
write_csv(out.stats, path = outFile)


out.name <- paste("oneStepEnsembles", model.dir, "paramsFrom", paramsFrom, "ticksFrom", ticksFrom, ua, sep = "_")
rdata.name <- paste0(out.name, ".RData")
outFile <- file.path(save.dir, rdata.name)

message("Saving simulation")
save(list = c("out", "tick.actual", "paramsFrom", "ticksFrom", "ua"), file = outFile)

message("Calculating Scores")

out.latent.larva <- out[q.index[1],-1,1,]
out.latent.nymph <- out[q.index[2],-1,1,]
out.latent.adult <- out[q.index[3],-1,1,]
out.obs.larva <- out[q.index[1],-1,2,]
out.obs.nymph <- out[q.index[2],-1,2,]
out.obs.adult <- out[q.index[3],-1,2,]

obs <- data$y[,-1]

tick.actual <- tick.actual %>% slice(-1)

crps.tibble <- tibble(
  time = tick.actual$DATE,
  model = model,
  paramsFrom = paramsFrom,
  ticksFrom = ticksFrom,
  uaExperiment = ua,
  obsLarvae = tick.actual$n_larvae,
  obsNymphs = tick.actual$n_nymphs,
  obsAdults = tick.actual$n_adults,
  crps.latent.larva = crps_sample(obs[q.index[1],], out.latent.larva),
  crps.latent.nymph = crps_sample(obs[q.index[2],], out.latent.nymph),
  crps.latent.adult = crps_sample(obs[q.index[3],], out.latent.adult),
  crps.obs.larva = crps_sample(obs[q.index[1],], out.obs.larva),
  crps.obs.nymph = crps_sample(obs[q.index[2],], out.obs.nymph),
  crps.obs.adult = crps_sample(obs[q.index[3],], out.obs.adult)
)

message("Saving scores")
out.name <- paste("crps", model.dir, "paramsFrom", paramsFrom, "ticksFrom", ticksFrom, ua, sep = "_")
csv.name <- paste0(out.name, ".csv")
outFile <- file.path(save.dir, csv.name)
write_csv(crps.tibble, path = outFile)



message(" --- DONE ---")




