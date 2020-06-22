# Script to generate simulation output for parameter sets generated via latin hypercube sampling (LHS) using the nlrx-package to access NetLogo with R


library(sp)
library(doSNOW)
library(parallel)
library(foreach)
library(here)
library(nlrx)
library(tidyverse)

#read some helper functions (shift, get.overlap.df ...)
source(here("code_memory_model", "rscripts", "analysis_tools.R"))

# read or create new LHS design
if(file.exists(here("LHS_design.csv"))){
  params <- sample_q <- read.csv(here("LHS_design.csv"), row.names = 1)
}else{
  params <- sample_q <-  lhs::improvedLHS(n = 1000, k = 5)
  write.csv(params, here("LHS_design.csv"))
}

# fill LHS design with parameter values by defining the lower and upper limit of the parameter space
# ninds, ITV, patchiness, species-1-mean, resource cover
distros <- data.frame("upper" = c(1000, .3, 100, 1, 70),
                      "lower" = c(100, .05, 0, 0.45, 15))

for(i in seq(nrow(distros)))  params[,i] <- qunif(sample_q[,i], min = distros[i,2], max = distros[i,1])

# write LHS design to simulations folder
params <- data.frame(params)
params$ID <- seq(nrow(params))
data.table::fwrite(params, here("simulations", today, experiment.dir, "LHS.csv"))


#define simulation setup (all parameter combinations will be tested)
experiment_type <- "HR_analysis"

#today <- substr(stringi::stri_replace_all_regex(as.character(Sys.Date()), "-", ""), 3, 9)
today <- "200514"
dir.create(here("Simulations", today))

# generate new experiment folder to avoid overwriting
experiment_count <- 1
while(file.exists(here("Simulations", today, paste0(experiment_type, experiment_count)))) experiment_count <- experiment_count + 1

experiment.dir <- paste0(experiment_type, experiment_count)
outpath <- here("Simulations", today, experiment.dir)
dir.create(outpath)

#define simulation name and date
simulation_day <- today
simulation_name <- experiment.dir
simulation_name_proc <- paste0("processed_", simulation_name)

#define simulation path
sim.path <- (here("simulations", simulation_day, simulation_name))
sim.path.proc <- here("simulations", simulation_day, simulation_name_proc)

# define path to model and NetLogo 
netlogopath <- paste0(dirname(here()), "/NetLogo 6.1.0")
modelpath <- file.path(here("code_memory_model", "model", "S2-Model.nlogo"))
eval.seq = c(1000, 2000, 50)


### TEST RUN
j <- 1 # first parameter set

run_sim(netlogopath = netlogopath, modelpath = modelpath, outpath = here("simulations", today, experiment.dir), values = params[j,], experiment_type = experiment_type)

sim_output <- readr::read_csv(paste0(outpath, "/", j, ".csv"))
head(sim_output) # does everything look fine?



### EXECUTE SIMULATIONS IN PARALLEL (adjust number of cores)!
{
  cl <- makeCluster(16)
  registerDoSNOW(cl)
  sims_out <-  foreach(.combine = rbind, i = seq(nrow(params)), .packages = c("here", "nlrx", "sp", "raster"))%dopar%{
    run_sim(netlogopath = netlogopath, modelpath = modelpath, outpath = here("simulations", today, experiment.dir), values = params[i,], experiment_type = experiment_type)
  }
  stopCluster(cl) 
}

