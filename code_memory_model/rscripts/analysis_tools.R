###############################
########### SIMULATION ###########
###############################

# setup the nl object and run the model

# i: ID of the simulation
# netlogopath: path to the netlogo.exe
# modelpath: path to the model.nlogo
# outpath: optional path to save output
# experiment_type: name of the experiment
# values: behaviour space given by latin hypercube sampling
# eval.seq: time steps with observation

run_sim <- function(netlogopath, modelpath, outpath, experiment_type, values = params, eval.seq = c(1000, 2000, 50)){
  if(!paste0(i, ".csv") %in% list.files(outpath)){
    # setup the nl object
    nl <- nl(nlversion = "6.1.0", nlpath = netlogopath, modelpath = modelpath, jvmmem = 1024)
    
    # setup the experiment
    nl@experiment <- experiment(expname = experiment_type, 
                                outpath = outpath,
                                repetition = 1,
                                tickmetrics = "true",
                                idsetup = "setup",
                                idgo = "go",
                                runtime = eval.seq[2],
                                evalticks = round(seq(eval.seq[1], eval.seq[2], length.out = eval.seq[3])),
                                metrics.turtles = list("turtles" = c("who", "pxcor", "pycor", "alpha", "mem-resources")),
                                constants = list("resource-cover" = values[1,5],
                                                 "n-inds" = values[1,1],
                                                 "max-ticks" = eval.seq[2] + 1,
                                                 "patchiness" = values[1,3],
                                                 "species-1-mean" = values[1,4],
                                                 "max-output" = F,
                                                 "n-species" = 1,
                                                 "set-this-seed" = 0,
                                                 "population-dynamics" = "false",
                                                 "ITV" = values[1,2]))
    
    # define the simulation design (simple = only constants required)
    nl@simdesign <- simdesign_simple(nl = nl, nseeds = 1)
    
    # run the nl object
    results <- run_nl_all(nl)
    
    # get the output, convert the output to points
    setsim(nl, "simoutput") <- results
    results.sf <- nlrx::nl_to_points(nl, coords = "px")
    
    
    # format the output
    run_output <- do.call(rbind, results.sf$spatial.turtles)
    run_output <- cbind(data.frame(run_output)[, c(2, 3,4, 6)], sf::st_coordinates(run_output))
    names(run_output) <- c("WHO", "ALPHA", "RESOURCEG",  "DATETIME", "XCOR", "YCOR")
    
    # do all individuals have an equal number of points:
    inds.with.n.of.points <- summary(as.factor(run_output$WHO), maxsum = 1000) == eval.seq[3]
    
    if(!all(inds.with.n.of.points)){
      warning("unequal number of points!")
      run_output <- run_output[!run_output$WHO %in% (which(!inds.with.n.of.points) - 1), ]
    }
    
    # assign the id to the output and return it
    run_output$ID <- values[1, 6]
    
    run_output <- cbind(run_output , (data.frame(values[1, 1:5])))
    names(run_output)[8:12] <- c("NINDS", "ITV", "PATCHINESS", "TRAITMEAN", "RESOURCECO")
    
    data.table::fwrite(run_output, paste0(outpath, "/", values[1, 6], ".csv"))
    
    return(run_output)
    
  }
}





#################################################
########### PROCESS SIMULATION OUTPUT ###########
#################################################

# read output from run sim to calculate home range sizes, home range overlap and the total distance moved

# run_output: spdf of relocations 
# eval.seq: the sequence of time steps with observations
# kernel: "epa" or "bvnorm" as estimator for the kernel parameter h of the utilization distribution
# grid.size: resolution of the utilization distribution
run_process <- function(input.folder, output.folder,  eval.seq = c(1000, 2000, 50), kernel, grid.size){
  
  # calculate home range metrics for different files in random order (random order to allow for parallelization with minimal inference) as loop to 
  for(input.file in sample(list.files(input.folder, full.names = T)[list.files(input.folder) != "LHS.csv"])){
    
    if(!basename(input.file) %in% list.files(output.folder)){
      
      # create a placeholder so no other core starts with processing the same simulation
      data.table::fwrite(list("placeholder"), paste0(output.folder, "/", basename(input.file)))
      
      # read the simulation file
      run_output <- read.csv(input.file)  
      coordinates(run_output) <- ~XCOR + YCOR  
      
      
      # process the simulation output on a toroid landscape to calculate true overlaps and home range sizes
      processed <- process.simulation(file = run_output, track.dist = T, landscape.size = 82, save = F, overlap.analysis = T, output.dir = here())
      
      # select the movement points in the centre of 3 x 3 replicated landscape, the overlap of individuals in the centre with all other will be calculated
      middle_inds <- processed[processed$shift == 6, ]
      
      #### get the individuals the have a overlap with the merged 95-% MCP of individuals in the centre and keep only these for the further anaylsis
      
      # 95-% MCP of individuals in the centred landscape
      mcp_middle_inds <- adehabitatHR::mcp(middle_inds[,"new_id"], percent = 95)
      
      # 95-% MCP of individuals in all landscapes
      mcp_all <- adehabitatHR::mcp(processed[,"new_id"], percent = 95)
      
      # individuals that intersect in their 95-% MCP with the 95-% MCP of individuals in the centred landscape
      keep.these <- as.character(mcp_all$id[rgeos::gIntersects(mcp_all, rgeos::gUnaryUnion(mcp_middle_inds), byid = T)])
      
      # 95-% MCP of individuals in all landscapes that intersects 95-% MCP of individuals in the centred landscape
      mcp.ss <- mcp_all[mcp_all$id %in% keep.these, ]
      
      # relocations of individuals in all landscapes that intersects 95-% MCP of individuals in the centred landscape
      ss.processed <- processed[processed$new_id %in% keep.these, ]
      
      # for each individual, get the individuals that overlap its 95-% MCP as a neighborhood-matrix
      nb.mat <- rgeos::gIntersects(mcp.ss, byid = T)
      
      # intialize output data.frame
      ol <- data.frame(matrix(nc = 6, nr = length(mcp_middle_inds$id)))
      
      # calculate the utilization distribution (UD) for all individuals that intersect the 95-% MCP of individuals in the centred landscape
      #UDs <- adehabitatHR::kernelUD(ss.processed[, "new_id"], grid = 250, same4all = T)
      
      #image(UDs.ss)
      
      for(j in seq(length(as.character(mcp_middle_inds$id)))){
        
        # select the focal_inds from the middel_inds
        focal_ind <- as.character(mcp_middle_inds$id)[j]
        
        # get the UDs of the focal individual and its neighbors
        
        # UDs.ss <- UDs[which(nb.mat[row.names(nb.mat) == focal_ind,])]
        inds <- names(which(nb.mat[row.names(nb.mat) == focal_ind,]))
        selected_relocations <- ss.processed[ss.processed$new_id %in% inds,]
        
        # calculate the utilization distribution (UD) of intersecting individuals
        UDs.ss <- adehabitatHR::kernelUD(selected_relocations[, "new_id"], grid = grid.size, same4all = T, kern = kernel)
        class(UDs.ss) <- "estUDm"
        
        # calculate home range overlaps
        if(length(UDs.ss) > 1){
          HR50 <- adehabitatHR::kerneloverlaphr(UDs.ss, method = "HR", percent = 50)
          diag(HR50) <- NA
          HR50[HR50 == 0] <- NA
          HR50 <- mean(HR50[row.names(HR50) == focal_ind, ], na.rm = T)
          
          HR95 <- adehabitatHR::kerneloverlaphr(UDs.ss, method = "HR", percent = 95)
          diag(HR95) <- NA
          HR95[HR95 == 0] <- NA
          HR95 <- mean(HR95[row.names(HR95) == focal_ind, ], na.rm = T)
        }else{
          HR50 <- 0
          HR95 <- 0
        }
        
        #select focal UD and calculate areas and get behaviour type (alpha) and distance moved and combine everything in one data.frame
        focal_UD <- UDs.ss[[focal_ind]]
        ol[j, ] <- data.frame(selected_relocations[selected_relocations$new_id == focal_ind, ]$ALPHA[1],
                              HR50,
                              HR95,
                              adehabitatHR::kernel.area(focal_UD, percent = 50),
                              adehabitatHR::kernel.area(focal_UD, percent = 95),
                              selected_relocations[selected_relocations$new_id == focal_ind, ]$dist[1]
        )
      }
      
      
      names(ol) <- c("ALPHA", "HR50", "HR95", "SIZE50", "SIZE95", "DISTMOVED")
      ol_out <- setNames(data.frame(ol, processed@data[1, c(1, 5:9)], eval.seq[2], eval.seq[1], eval.seq[3]), c(names(ol), "ID", "n_inds", "ITV", "patchiness", "trait_mean", "resource_cover", "max_ticks", "obs_start", "n_samples"))
      ol_out$ID <- as.numeric(tools::file_path_sans_ext(basename(input.file)))
      data.table::fwrite(ol_out, paste0(output.folder, "/", basename(input.file)))
    }
  }
}


###############################
########### GEOTOOL ###########
###############################


# shift the landscape in all 8 neighboring directions to work with home ranges in a torus environment

# shp: relocations of the spdf
# ls.size: size of the landscape (default is 82)

shift <- function(shp, ls.size) {
  combs <- combinat::combn(rep(seq(-ls.size, ls.size, ls.size), 2), 2)
  shifts <- combs[, !duplicated(paste0(combs[2,], combs[1,]))]
  shift.shp <- replicate(9, shp)

  for(i in seq(9)){
    shift.shp[[i]]@coords[,1] <- shift.shp[[i]]@coords[,1] + shifts[1, i]
    shift.shp[[i]]@coords[,2] <- shift.shp[[i]]@coords[,2] + shifts[2, i]
  }
  
  shift.shp <- do.call(rbind, shift.shp)
  return(shift.shp)
}



# do the shift and stitch the individuals together
# file = location of the .shp-file or a SpatialPointsDataFrame object
# landscape.size = size of the landscape in NetLogo
# keep.this.cols = vector of characters that specify columns that should be returned, if "all" then all columns are returned
# output.dir = where should the processed movement data be stored
# save = boolean, save the output or not
# track.dist = boolean, measure distance between relocations
# overlap.analysis = boolean, shift the output at the end again to calculate overlaps

process.simulation <- function(file, landscape.size = 82, keep.this.cols = "all", output.dir = sim.path.proc, save = T, track.dist = T, overlap.analysis = T){
  
  min.range <- sp::Polygon(coords = matrix(c(0, 0, 0, landscape.size - 0.01, landscape.size - 0.01, landscape.size - 0.01, landscape.size - 0.01, 0), byrow = T, nc= 2))
  min.range <- as(SpatialPolygons(list(Polygons(list(min.range), ID = 1))), "SpatialPolygonsDataFrame")
  
  if(!file.exists(output.dir)) dir.create(output.dir, showWarnings = T)
  
  # determine whether file is a SpatialPointsDataFrame or a character string as the file location
  if(class(file)[1] == "SpatialPointsDataFrame")  shp <- file
  if(class(file)[1] == "character") shp <- rgdal::readOGR(file)
  if(!class(file)[1] %in% c("character", "SpatialPointsDataFrame")) warning("Wrong file argument: Must be a character or SpatialPointsDataFrame")

  # remove unwanted columns, if any
  if(keep.this.cols != "all") shp <- shp[, match(keep.this.cols, names(shp))]
  
  # shift the landscape 2 times creating 3 x 3 identical landscapes with the original landscape in the middle
  shp.shift <- shift(shp, landscape.size)

  # initialize output files and get the sequence of time steps between observed relocations
  start.DATETIME <- min(shp$DATETIME)
  step.DATETIME <- diff(unique(shp$DATETIME))
  all.dists <- all.positions <- setNames(vector("list", length = length(unique(shp.shift$WHO))), unique(shp.shift$WHO))

  # for each individual (WHO, ID)...
  for(ID in unique(shp.shift$WHO)){
    
    # ... get its locations, start position and time in the central landscape in 3 x 3 shifted landscape
    IND <- shp.shift[shp.shift$WHO == ID, ]
    starts <- IND[which(IND$DATETIME == start.DATETIME),]
    positions <- vector("list", length = length(step.DATETIME + 1))
    dists <-  vector("list", length = length(step.DATETIME))
    start.point <- positions[[1]] <-  raster::crop(starts, min.range)
    time <- start.point$DATETIME

    # .. iterate trough the positions...
    for(pos in seq(length(positions))){
      # ... and times ...
      time <- time + step.DATETIME[pos]
      # ... and determine the next step as the nearest location among the 9 possible options that have an equal datetime attribute
      nextstep <- IND[which(IND$DATETIME == time),]
      next.dist <- rgeos::gDistance(start.point, nextstep, byid = T)
      if(track.dist) dists[[pos]] <- min(next.dist)
      start.point <- positions[[pos+1]] <- nextstep[which.min(next.dist),]
    }

    # combine all positions
    all.positions[[as.character(ID)]] <- do.call(rbind, positions)
    
    # give an output for the distance between relocations
    if(track.dist) all.dists[[as.character(ID)]] <- c(do.call(rbind, dists))
  }

  # combine all positions of all individuals
  sim.torus <- do.call(rbind, all.positions)
  if(track.dist) sim.dists <- unlist(lapply(all.dists, sum))

  # if an overlap analysis is necessary, the output needs to be shifted again to account for overlap in the periodical landscape
  if(overlap.analysis){
    shift.sim.torus <- shift(sim.torus, landscape.size)
    shift.sim.torus$shift <- rep(1:9, each = nrow(sim.torus))
    shift.sim.torus$new_id <- paste0(shift.sim.torus$WHO, "_", shift.sim.torus$shift)
  }else{
    shift.sim.torus <- sim.torus
  }
  
  # add the distance between relocations to the output
  shift.sim.torus$dist <- sim.dists[match(shift.sim.torus$WHO, names(sim.dists))]
  
  # save the output to a local file or return it to the global environment
  if(save){
    rgdal::writeOGR(shift.sim.torus, output.dir, basename(file), "ESRI Shapefile", overwrite_layer = T, check_exists = T)
  } 
  
  if(!save) return(shift.sim.torus)
}


##################################
########### FORMATTING ###########
##################################

# read many .csv output-files quickly and optionally in parallel and format them as needed

# day: day of the simulation as "YYMMDDD"
# exp.name: name of the simulation 
# nofcores: number of cores 

read.sim <- function(day, exp.name, nofcores){
  if(!paste0(exp.name, "_all.csv") %in% list.files(here("simulations", day))){
    csvs <- list.files(here("simulations", day, exp.name), full.names = T, pattern =  ".csv$")
    this.csvs <- sapply(csvs, function(x) file.info(x)$size) > 1000 & stringi::stri_count_regex(csvs, "alpha") == 0
    csvs <- list.files(here("simulations", day, exp.name), full.names = F, pattern =  ".csv$")[this.csvs]
    
    
    
    {
      cl <- makeCluster(nofcores)
      registerDoSNOW(cl)
      all.dfs <-  foreach( .combine = rbind, i = seq(length(csvs)), .packages = "here")%dopar%{
        data.frame(read.csv(file = here("simulations", day, exp.name, csvs[i])), ID = i)
      }
      
      stopCluster(cl) 
    }
    
    data.table::fwrite(file = here("simulations", day, paste0(exp.name, "_all.csv")), x = all.dfs)
  }else{
    all.dfs <- readr::read_csv(file = here("simulations", day, paste0(exp.name, "_all.csv")))
  }
  return(all.dfs)
}


##################################
########### FORMATTING ###########
##################################




### Calculate different home range metrics - function deprecated

get.overlap.df <- function(data, UD, methods = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), focal.ind){
  # apply overlap methods to the UD

  i.focal <- which(names(UD) == "focal")
  i.others <- which(names(UD) == "others")
  
  HR50 <- adehabitatHR::kerneloverlaphr(UD, method = "HR", conditional = T, percent = 50)[i.focal, i.others]
  HR95 <- adehabitatHR::kerneloverlaphr(UD, method = "HR", conditional = T, percent = 95)[i.focal, i.others]

  # combine in one data.frame
  overlap.df <- 
    data.frame(new_id = focal.ind, 
               HR.SIZE = t(adehabitatHR::kernel.area(UD$focal, percent = 95)),
               CORE.SIZE = t(adehabitatHR::kernel.area(UD$focal, percent = 50)),
               HR50 = HR50,
               HR95 = HR95)
 
  merged.overlap <- merge(data, overlap.df, by = "new_id")
  return(merged.overlap)
}