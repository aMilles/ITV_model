# SCRIPT TO ANALYZE TO EXTRACTION OF RESOURCES PER CLUMP AND THE TOTAL NUMBER OF CLUMPS WITH EXTRACTED RESOURCES (Fig 4, B)

library(here)
library(raster)
library(ggplot2)
library(velox)
# veloxURL = "https://cran.r-project.org/src/contrib/Archive/velox/velox_0.2.0.tar.gz"
# install.packages(veloxURL, repos = NULL, type = "source")
#read file names for x y coordinates (csv) and landscape (png) - select sim with 40 or 80 inds
exp.name_ <- "final_trade_off_750inds"
day_ <- "200514"

pngs <- list.files(here("simulations", day_, exp.name_), pattern = ".png", full.names = T)
csvs <- list.files(here("simulations", day_, exp.name_), pattern = ".csv", full.names = T)

#create empty output 
out.list <- vector("list", length = length(csvs))

# check whether pngs and csvs match
pngs <- pngs[match(substr(stringi::stri_reverse(paste0(tools::file_path_sans_ext(basename(csvs)), "_world")), 0, 15),
                   substr(stringi::stri_reverse(paste0(tools::file_path_sans_ext(basename(pngs)), "")), 0, 15))]


# if the csvs and pngs have not been processed before, process and save the output, otherwise just read the output

if(!file.exists(here("simulations", day_, paste0(exp.name_, ".csv")))){

  # iterate over the respective png and csv files
  for(i in seq(length(pngs))){
    
    # NetLogo exports the world with a defined number of pixels per patch, this can be changed with the imager package or by reducing the view in NetLogo..
    
    show <- raster(pngs[i])
    
    
    # classify landscape in resource (1) and no resource (0)
    show <- show < 200
    
    #read png as raster and csv as data.frame
    
    track <- readr::read_csv(csvs[i])
    track$WHO <- as.character(track$WHO)
    
    #read the patchiness of the landscape, only use alpha-level and coordinates
    patchiness <- median(track$PATCHINESS)
    
    ind.ids <- unique(track$WHO)
    
    
    selected.steps <- setNames(vector("list", length = length(ind.ids)), ind.ids)
    for(ind in ind.ids){
      track.ind <- track[track$WHO == ind, ]
      selected.steps[[ind]] <- track.ind[which(diff(track.ind$RESOURCESG) == 1) + 1,]
    }
    
    track <- do.call(rbind, selected.steps)
    
    
    # convert object to spatial object, shift coordinates by 0.5
    subs <- track[, 1:3]
    coordinates(subs) <- ~XCOR+YCOR
    subs@coords <- subs@coords + 0.5
    

    
    
    # calculate the clumps and their size
    rclump <- clump(show)
    frclump <- freq(rclump)
    frclump[,2] <- frclump[,2]
    vrclump <- velox(rclump)
    
    
    # extract the id of the clump where the individual exploited a resource patch, we use a buffer of .9 spatial units around an individuals location to adjust for minor differences between coordinates
    buffered_subs <- rgeos::gBuffer(subs, byid = T, width = .9)
    
    this_clumps <- vrclump$extract(buffered_subs, function(x) median(x, na.rm = T))
    
    clumps <- lapply(this_clumps, function(x) median(x, na.rm = T))
    
    out <- data.frame(subs$ALPHA, clump.size = frclump[unlist(clumps),2], clump.id = unlist(clumps))
    
    out.list[[i]] <- out
    print(i)
  }
}
  
for(i in seq(length(out.list))) out.list[[i]]$ID <- i

all <- do.call(rbind, out.list)
  
data.table::fwrite(all, here("simulations", day_, paste0(exp.name_, "_all.csv")))
  

out.list <- split(all, all$ID)


for(i in seq(length(out.list))){
  out <- out.list[[i]]
  patchiness <- as.numeric(substr(stringi::stri_replace_all_regex(basename(csvs[i]), "generalist_analysis", ""), 26, 27))
  # return the number of resources that have been exploited 
  n.vis.clumps <- aggregate(clump.id ~ clump.size,  out, FUN = function(x) length(x))
  
  # maximum number of resources exploited from only one clump
  n.vis.clumps.max <- aggregate(clump.size ~ subs.ALPHA + clump.id,  out, FUN = function(x) length(x))
  n.vis.clumps.max <- aggregate(clump.size ~ subs.ALPHA, n.vis.clumps.max, median)
  
  # return the number of unique clumps
  n.diff.clumps <- aggregate(clump.id ~ subs.ALPHA,  out, FUN = function(x) length(unique(x))) 
  # return the median size of the utilized clumps
  med.clump.size <- aggregate(clump.size ~ subs.ALPHA,  out, FUN = function(x) median(x)) 
  
  # create the output data.frame
  out.list[[i]] <-  data.frame(ALPHA = n.diff.clumps$subs.ALPHA, n.clumps = n.diff.clumps$clump.id, clump.size = med.clump.size$clump.size, patchiness = patchiness, max.visits = n.vis.clumps.max$clump.size, patchiness = patchiness)
}

# combine and save output

all2 <- do.call(rbind, out.list)

all2$patchiness <- as.factor(all2$patchiness)

all2 <- all2[, c("ALPHA", "n.clumps", "patchiness", "max.visits")]

mall2 <- reshape2::melt(all2, id.vars = c("ALPHA", "patchiness"))


data.table::fwrite(mall2, here("simulations", day_, paste0(exp.name_, "_all.csv")))




