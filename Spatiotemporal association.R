######################################
###                                ###
###   SPATIOTEMPORAL ASSOCIATION   ###
###                                ###
######################################

#==========================================================
# Load Data for analyses
#==========================================================
dir <- "C:/Users/Owner/Documents/"
setwd(dir)
library(tidyverse)
library(zoo)
library(ReIns)
library(MASS)
library(progress)

# Load camera trap captures
# Original dataset contains one row per photo
# Reformatting so each row is an 'event' (independent detection) with start/stop times
# Time measured as seconds since camera period began
Captures <- read.csv("Captures.csv") %>%
  group_by(species,Event_ID) %>%
  mutate(time_start = min(sec_since_cam),
         time_end = max(sec_since_cam)) %>%
  slice(1) %>%
  dplyr::select(-X,-sec_since_cam) %>%
  dplyr::filter(!species %in% c("human","ALL_BLACK","ALL_WHITE")) # Filter out useless photos

# Load camera trap periods
# Each camera trap period has a unique identifier (Period_unique)
# Logged (1) and unlogged (0) camera sites are identified
# Mast (1) and non-mast (0) years are identified
Periods_all <- read.csv("CT_periods.csv")

#==========================================================
# Single test of spatiotemporal association
#==========================================================

# Set responder species and inducer species
# Analyze detections of responder after detections of inducer 
resp <- "sambar" # the 'responder' species
indu <- "bearded_pig" # the'inducer' species

# Set which subset of data to anlayze
Periods <- Periods_all #%>% dplyr::filter(mast == 1, logged==1)

#----------------------------------------------------------
# # Observed sightings of responder after inducer
#----------------------------------------------------------

# Create dataframe of all time differences of responder after inducer
obs_df <- Captures %>%
  dplyr::filter(species %in% c(resp,indu)) %>% # exclude other species
  dplyr::inner_join(Periods, by = "Period_unique") %>%
  dplyr::arrange(Period_unique, time_start) %>%
  dplyr::group_by(Period_unique) %>%
  dplyr::mutate(indu_time = replace(time_end, species == resp, NA), # Identify end time of each inducer detection
         most_recent_indu = zoo::na.locf(indu_time, na.rm = F), # Identify most recent inducer detection
         sec_since_indu = time_start - most_recent_indu, # Calculate number of seconds since most recent inducer detection
         day_since_indu = sec_since_indu/(24*3600)) %>% # Convert seconds to days
  dplyr::filter(species == resp, # Only include responder detections that occurred after an inducer detection
                sec_since_indu > 0)

# Save list of observed time differences
obs_all <- obs_df$day_since_indu
# Exclude time differences less than one minute
# (one detection of one species with multiple photos sometimes misidentified as multiple species when multiple photos IDed separately)
obs_all <- obs_all[obs_all > (1/60/24)]

#----------------------------------------------------------
# # Expected (simulated) sightings of responder after inducer
#----------------------------------------------------------

# Set number of simulated observations of responder
N_sim <- 50000

# Create list of all possible times since inducer that responder could have been detected
# "sec_since_indu" is duration of time (seconds) from inducer detection to either next inducer detection
#     or to end of period (if no other inducer detections before end of period)
indu_times <- Captures %>%
  dplyr::arrange(Period_unique, time_start) %>%
  dplyr::filter(species == indu) %>%
  dplyr::inner_join(Periods, by = "Period_unique") %>%
  dplyr::group_by(Period_unique) %>%
  dplyr::mutate(next_indu = lead(time_start), # Calculate sec_since_indu in each period 
                sec_since_indu = ifelse(is.na(next_indu),
                                        Duration_sec - time_end, # If last inducer detection, calcuclate number of seconds until end of period
                                        next_indu - time_end)) %>% # If not last inducer detection, calcuclate number of seconds until next inducer
  ungroup() %>%
  mutate(sumall = sum(sec_since_indu),
         indu_time_prop = sec_since_indu/sumall) # For each span of time after inducer sightings, calculate duration as proportion of total seconds

# Simulate responder detections after inducer
# First randomly select "sec_since_indu" (row number), weighted by duration of "sec_since_indu"
random_event_indu <- sample(1:nrow(indu_times), size = N_sim, replace = TRUE,
                       prob = indu_times$indu_time_prop)
# Then randomly select time (number of seconds) within each "sec_since_indu"
sim_all <- c()
for(i in 1:N_sim){
  sim_all <- c(sim_all, sample(1:indu_times$sec_since_indu[random_event_indu[i]], size = 1))
}
sim_all <- sim_all/(24*3600) # Convert seconds to days
sim_all <- sim_all[sim_all > (1/60/24)] # Exclude time differences less than one minute (same as observed dataset)

#----------------------------------------------------------
# # Comparing observed and expected/simulated spatiotemporal patterns
#----------------------------------------------------------

# Create functions for truncated Weibull distributions
# Log-likelihood
tweib_logl <- function(x, par, endpoint){
  # par = c(shape, scale)
  # when shape = 1, scale = lambda of exponential dist
  B <- exp(par)
  tmp <- ReIns::dtweibull(x, B[1], B[2], endpoint)
  logL <- sum(log(tmp))
  return(logL)
}
# Fit Weibull distribution
fit_tw <- function(x, endpoint){
  out <- optim(par = log(c(1, mean(x, na.rm = T))),
               endpoint = endpoint,
               tweib_logl,
               x = x,
               control = list(fnscale = -1),
               hessian = T)
  return(out)
}

# Comparing two distributions, following advice from Ben Bolker
# # https://r.789695.n4.nabble.com/Comparison-of-two-weibull-distributions-td4679632.html
# P-value of the difference, according to likelihood ratio test
p_val <- function(fit_fn, dat1, dat2, endpoint){
  # 1. Fit model to combined data.
  tog_f <- fit_fn(c(dat1, dat2), endpoint)
  # 2. Fit separate models to each data set on its own
  dat1_f <- fit_fn(dat1, endpoint)
  dat2_f <- fit_fn(dat2, endpoint)
  # 3. Compare log-likelihood of pooled model to sum of log-lik #2
  logl_sum <- dat1_f$value + dat2_f$value
  logl_pooled <- tog_f$value
  # P-value of the difference, according to likelihood ratio test
  # (as n -> Inf, becomes Chi2)
  p <- pchisq(2 * (logl_sum - logl_pooled),
              df = 2,
              lower.tail = FALSE)
  return(p)
}

# Set truncation point for analyses
trunc <- 14 # 14 days
# Exclude detections past truncation point
obs <- obs_all[obs_all<=trunc]
sim <- sim_all[sim_all<=trunc]

# Calculate parameters for obs and sim and p-value for comparison of obs vs. sim
est_obs <- fit_tw(obs, trunc)
est_sim <- fit_tw(sim, trunc)
vc_obs <- -1 * MASS::ginv(est_obs$hessian)
vc_sim <- -1 * MASS::ginv(est_sim$hessian)
length(obs) # Number of observed responder detections after inducer within truncation limit
p_val(fit_tw, obs, sim, trunc) # p-value for obs vs. sim comparison
exp(est_obs$par[1]) # k observed
msm::deltamethod(g = ~exp(x1), mean = est_obs$par[1],
                 cov = vc_obs[1, 1], ses = T) # SE of k observed
exp(est_obs$par[2]) # lambda observed
msm::deltamethod(g = ~exp(x1), mean = est_obs$par[2],
                 cov = vc_obs[1, 1], ses = T) # SE of lambda observed
exp(est_sim$par[1]) # k expected
msm::deltamethod(g = ~exp(x1), mean = est_sim$par[1],
                 cov = vc_sim[1, 1], ses = T) # SE of k expected
exp(est_sim$par[2]) # lambda expected
msm::deltamethod(g = ~exp(x1), mean = est_sim$par[2],
                 cov = vc_sim[1, 1], ses = T) # SE of lambda expected
(exp(est_obs$par[1])-exp(est_sim$par[1]))/exp(est_sim$par[1])*100 # % change k

# Plot both curves
# Expected curve
curve(dweibull(x, shape=exp(fit_tw(sim, trunc)$par[1]),
               scale=exp(fit_tw(sim, trunc)$par[2])),
      from=0.1, to=trunc, ylab="Density",
      lty=3, lwd=2, cex.axis=1.3, las=1, cex.lab=1.5,
      xlab="Days since inducer", ylim=c(0,.8))
# Observed curve
curve(dweibull(x, shape=exp(fit_tw(obs, trunc)$par[1]),
               scale=exp(fit_tw(obs, trunc)$par[2])),
      from=0.1, to=trunc, lwd=2, add=TRUE)

#==========================================================
# Pairwise comparisions of spatiotemporal association for list of species
#==========================================================

# Set number of simulated observations of responder
N_sim <- 50000

# Set truncation number of days
trunc <- 14

# Set functions needed for analysis (same as in "Single test of spatiotemporal association')
# Log-likelihood
tweib_logl <- function(x, par, endpoint){
  # par = c(shape, scale)
  # when shape = 1, scale = lambda of exponential dist
  B <- exp(par)
  tmp <- ReIns::dtweibull(x, B[1], B[2], endpoint)
  logL <- sum(log(tmp))
  return(logL)
}
# Fit Weibull distribution
fit_tw <- function(x, endpoint){
  out <- optim(par = log(c(1, mean(x, na.rm = T))),
               endpoint = endpoint,
               tweib_logl,
               x = x,
               control = list(fnscale = -1),
               hessian = T)
  return(out)
}
# P-value comparing two distributions (likelihood ratio test)
p_val <- function(fit_fn, dat1, dat2, endpoint){
  # 1. Fit model to combined data.
  tog_f <- fit_fn(c(dat1, dat2), endpoint)
  # 2. Fit separate models to each data set on its own
  dat1_f <- fit_fn(dat1, endpoint)
  dat2_f <- fit_fn(dat2, endpoint)
  # 3. Compare log-likelihood of pooled model to sum of log-lik #2
  logl_sum <- dat1_f$value + dat2_f$value
  logl_pooled <- tog_f$value
  # P-value of the difference, according to likelihood ratio test
  # (as n -> Inf, becomes Chi2)
  p <- pchisq(2 * (logl_sum - logl_pooled),
              df = 2,
              lower.tail = FALSE)
  return(p)
}

# Set which subset of data to anlayze
Periods_subset <- Periods_all #%>% dplyr::filter(logged == 1)
# If inducer or responder is argus, remove camera/year with lekking
Periods_argus <- Periods_subset %>% dplyr::filter(!Period_unique %in% c("AG16.19.1","AG16.19.2"))

# Make dataframes to fill in
# Row names are responder, column names are inducer
species_list <- c("bearded_pig","mousedeer","yellow_muntjac", "sambar",
                  "argus","pig_tailed_macaque","fireback","malay_civet",
                  "banded_civet")
results_mat <- matrix(nrow=length(species_list),ncol=(length(species_list)),
                      dimnames = list(species_list,species_list))
N_obs <- results_mat
P_vals <- results_mat
k_obs <- results_mat
k_obs_se <- results_mat
lambda_obs <- results_mat
lambda_obs_se <- results_mat
k_sim <- results_mat
k_sim_se <- results_mat
lambda_sim <- results_mat
lambda_sim_se <- results_mat
change_k <- results_mat

#----------------------------------------------------------
# # Run for loop to compile the results
#----------------------------------------------------------

# Run through all possible pairwise comparisons of inducer/responder
pb <- progress::progress_bar$new(
  format = " running [:bar] :percent in :elapsed eta: :eta",
  total = prod(dim(N_obs)), clear = FALSE, width = 100) # Progress bar
for(i in 1:nrow(results_mat)){
  resp <- rownames(results_mat)[i] # the responder
  for(j in 1:ncol(results_mat)){
    indu <- colnames(results_mat)[j] # the inducer
    Periods <- if("argus" %in% c(resp,indu)) Periods_argus else Periods_subset
    # If inducer and responder are same species, fill results with NA
    if(resp==indu) {
      N_obs[i,j] <- NA
      P_vals[i,j] <- NA
      k_obs[i,j] <- NA
      k_obs_se[i,j] <- NA
      lambda_obs[i,j] <- NA
      lambda_obs_se[i,j] <- NA
      k_sim[i,j] <- NA
      k_sim_se[i,j] <- NA
      lambda_sim[i,j] <- NA
      lambda_sim_se[i,j] <- NA
    }
    else {
      #-----------------------------------------------
      # Observed sightings of responder after inducer
      #-----------------------------------------------
      
      # Create dataframe of all time differences of responder after inducer
      obs_df <- Captures %>%
        dplyr::filter(species %in% c(resp,indu)) %>% # exclude other species
        dplyr::inner_join(Periods, by = "Period_unique") %>%
        dplyr::arrange(Period_unique, time_start) %>%
        dplyr::group_by(Period_unique) %>%
        dplyr::mutate(indu_time = replace(time_end, species == resp, NA), # Identify end time of each inducer detection
                      most_recent_indu = zoo::na.locf(indu_time, na.rm = F), # Identify most recent inducer detection
                      sec_since_indu = time_start - most_recent_indu, # Calculate number of seconds since most recent inducer detection
                      day_since_indu = sec_since_indu/(24*3600)) %>% # Convert seconds to days
        dplyr::filter(species == resp, # Only include responder detections that occurred after an inducer detection
                      sec_since_indu > 0)
      
      # Save list of observed time differences
      obs_all <- obs_df$day_since_indu
      # Exclude time differences less than one minute
      # (one detection of one species with multiple photos sometimes misidentified as multiple species when multiple photos IDed separately)
      obs_all <- obs_all[obs_all > (1/60/24)]
      # Exclude detections past truncation point
      obs <- obs_all[obs_all <= trunc]
      
      #-----------------------------------------------
      # Expected (simulated) sightings of responder after inducer
      #-----------------------------------------------
      
      # Create list of all possible times since inducer that responder could have been detected
      # "sec_since_indu" is duration of time (seconds) from inducer detection to either next inducer detection
      #     or to end of period (if no other inducer detections before end of period)
      indu_times <- Captures %>%
        dplyr::arrange(Period_unique, time_start) %>%
        dplyr::filter(species == indu) %>%
        dplyr::inner_join(Periods, by = "Period_unique") %>%
        dplyr::group_by(Period_unique) %>%
        dplyr::mutate(next_indu = lead(time_start), # Calculate sec_since_indu in each period 
                      sec_since_indu = ifelse(is.na(next_indu),
                                              Duration_sec - time_end, # If last inducer detection, calcuclate number of seconds until end of period
                                              next_indu - time_end)) %>% # If not last inducer detection, calcuclate number of seconds until next inducer
        ungroup() %>%
        mutate(sumall = sum(sec_since_indu),
               indu_time_prop = sec_since_indu/sumall) # For each span of time after inducer sightings, calculate duration as proportion of total seconds
      
      # Simulate responder detections after inducer
      # First randomly select "sec_since_indu" (row number), weighted by duration of "sec_since_indu"
      random_event_indu <- sample(1:nrow(indu_times), size = N_sim, replace = TRUE,
                                  prob = indu_times$indu_time_prop)
      # Then randomly select time (number of seconds) within each "sec_since_indu"
      sim_all <- c()
      for(x in 1:N_sim){
        sim_all <- c(sim_all, sample(1:indu_times$sec_since_indu[random_event_indu[x]], size = 1))
      }
      sim_all <- sim_all/(24*3600) # Convert seconds to days
      sim_all <- sim_all[sim_all > (1/60/24)] # Exclude time differences less than one minute (same as observed dataset)
      # Exclude detections past truncation point
      sim <- sim_all[sim_all <= trunc]
      
      #-----------------------------------------------
      # Expected (simulated) sightings of responder after inducer
      #-----------------------------------------------
      
      est_obs <- fit_tw(obs, trunc)
      est_sim <- fit_tw(sim, trunc)
      vc_obs <- -1 * MASS::ginv(est_obs$hessian)
      vc_sim <- -1 * MASS::ginv(est_sim$hessian)
      N_obs[i,j] <- length(obs)
      P_vals[i,j] <- p_val(fit_tw, obs, sim, trunc)
      k_obs[i,j] <- exp(est_obs$par[1]) # Estimate of k
      k_obs_se[i,j] <- msm::deltamethod(g = ~exp(x1), 
                                        mean = est_obs$par[1],
                                        cov = vc_obs[1, 1],
                                        ses = T) # Delta method for standard error of k
      lambda_obs[i,j] <- exp(est_obs$par[2]) # Estimate of lambda
      lambda_obs_se[i,j] <- msm::deltamethod(g = ~exp(x1), 
                                             mean = est_obs$par[2],
                                             cov = vc_obs[1, 1],
                                             ses = T) # Delta method for standard error of lambda
      k_sim[i,j] <- exp(est_sim$par[1]) # Estimate of k
      k_sim_se[i,j] <- msm::deltamethod(g = ~exp(x1), 
                                        mean = est_sim$par[1],
                                        cov = vc_sim[1, 1],
                                        ses = T) # Delta method for standard error of k
      lambda_sim[i,j] <- exp(est_sim$par[2]) # Estimate of lambda
      lambda_sim_se[i,j] <- msm::deltamethod(g = ~exp(x1), 
                                             mean = est_sim$par[2],
                                             cov = vc_sim[1, 1],
                                             ses = T) # Delta method for standard error of lambda
      change_k[i,j] <- 100*((k_obs[i,j]/k_sim[i,j])-1)
    }
    pb$tick()
  }
}

#----------------------------------------------------------
# # View results
#----------------------------------------------------------

View(N_obs) # number of observed responder detections after inducer within truncation limit
View(P_vals) # p-value for obs vs. sim comparison
View(k_obs) # k observed
View(k_obs_se) # SE of k observed
View(lambda_obs) # lambda observed
View(lambda_obs_se) # SE of lambda observed
View(k_sim) # k expected
View(k_sim_se) # SE of k expected
View(lambda_sim) # lambda expected
View(lambda_sim_se) # SE of lambda expected
View(change_k) # % change k (observed vs. expected)
