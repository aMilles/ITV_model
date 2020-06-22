library(here)
library(raster)
library(doSNOW)
library(foreach)
library(parallel)

### THIS SCRIPT CALCULATES MOVEMENT METRICS 
### OVERLAP, DISTANCE MOVED, HOME RANGE SIZE
### FOR SIMULATION RUNS OF SCRIPT NO 5

source(here("code_memory_model", "rscripts", "analysis_tools.R"))


### TEST RUN

run_output <- read.csv(here("simulations", "190917", "HR_analysis_test_metrics.csv"))
coordinates(run_output) <- ~ XCOR + YCOR

# THIS PATH NEEDS TO BE ADJUSTED (HAS TO SET WITHOUT HERE() AS SIMULATIONS WERE IN ANOTHER FOLDER)
in.path <- "Y:/Gruppen/oesa/milles/proj_MemoryModel/simulations/200514/HR_analysis4"

# path where movement metrics should be stored
out.path <- here("simulations", "200514", "HRanalysis_processed_1000")
dir.create(out.path)

# check if calculation of movement metrics works for one simulation
test <- read_csv(list.files(in.path, full.names = T)[1])

# run function to calculate movement metrics
run_process(input.folder = in.path, output.folder = out.path, eval.seq = c(1000, 2000, 50), kernel = "bivnorm", grid.size = 250)


### DO IT PARALLEL WHILE SIMULATIONS ARE PROCESSED (process_SimulationInR is running)
read.csv(list.files(out.path, full.names = T)[2])
while(T){
  {
    cl <- makeCluster(9)
    registerDoSNOW(cl)
    sims_out <-  foreach(i = seq(9), .packages = c("here", "sp", "raster"))%dopar%{
      run_process(input.folder = in.path, output.folder = out.path, eval.seq = c(1000, 2000, 50), kernel = "bivnorm", grid.size = 250)
    }
    stopCluster(cl) 
  }
  print('all done for now')
  Sys.sleep(500)
}
  

proc.files <- list.files(out.path, full.names = T)
in.files <- list.files(in.path, full.names = T)


# create an additional file where processed output is stored with parameter settings

out.path2 <- here("simulations", "200514", "HRanalysis_processed_1000_fit")
dir.create(out.path2)

proc.file <- proc.files[1]
for(proc.file in proc.files){
  in.file <- in.files[basename(in.files) == basename(proc.file)]
  
  in.data <- read_csv(in.file, col_types = cols())
  proc.data <- read_csv(proc.file, col_types = cols())
  
  proc.data$trait_mean <- in.data$TRAITMEAN[1]
  proc.data$resource_cover <- in.data$RESOURCECO[1]
  proc.data$ITV <- in.data$ITV[1]
  proc.data$patchiness <- in.data$PATCHINESS[1]
  proc.data$n_inds <- in.data$NINDS[1]
 
  write_csv(proc.data,  here("simulations", "200514", "HRanalysis_processed_1000_fit", basename(proc.file))) 
}

