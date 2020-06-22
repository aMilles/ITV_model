# Script to generate Fig. 5 (and supplements)

library(ggplot2)
library(raster)
library(here)
library(doSNOW)
library(parallel)
library(foreach)
library(cowplot)
library(tidyverse)
library(sp)
source(here("code_memory_model", "rscripts", "analysis_tools.R"))

# read alle simulations for the different landscapes (heterogeneous  = patchy, homogeneous = random)
day_ <- "200514"
exp.name_ <- "final_coviability_analysis_heterogeneous"


files <- NULL
for(exp.name_ in c("final_coviability_analysis_heterogeneous",
                   "final_coviability_analysis_intermediate",
                   "final_coviability_analysis_homogeneous")){
  files <- append(files, list.files(here("simulations", day_, exp.name_), full.names = T))
}

list.sims <- lapply(as.list(files), readr::read_csv)

list.sims <- lapply(list.sims, function(x){
  if(NROW(x) > 0){
    x$ID <- runif(1, max = 100)
    return(x)
  }
} )

all.sims <- do.call(rbind, list.sims)

# iterate over different patchiness levels to generate the specific figure for supplement / main text

for(patchy_level in unique(all.sims$PATCHINESS)){
  sim.dfs <- all.sims %>% filter(PATCHINESS == patchy_level)
  
  #calculate the mean beaviour type (trait) for each species in each simulation
  meantrait <- setNames(aggregate(ALPHA ~ ID + SPECIESID + ITV, sim.dfs, mean), c("ID", "SPECIESID", "ITV",  "SPECIES_TRAIT"))
  
  #split focal and competing species
  split.meantrait <- split(meantrait[, -2], meantrait$SPECIESID)
  
  #setup data.frame for each species
  for(i in seq(length(split.meantrait))) assign(paste0("trait", i), setNames(split.meantrait[[i]], c("ID", "ITV", paste0("value", i))))
  
  #merge the data frames with the information about landscape structure, community from sim.dfs
  merged.out <- merge(merge(sim.dfs, trait1[,c("ID", "value1")]), trait2[,c("ID", "value2")])
  
  # convert ITV parameter to a boolean
  merged.out$ITV <- merged.out$ITV > 0
  
  
  #repeat the following for ITV (T) and no ITV (F) scenarios
  
  res <- 0.5 / 15
  
  # loop for simulations with and without ITV
  for(bool in c(T,F)){
    
    # create a raster
    ext <- extent(c(0.25,.75,0.25,.75))
    rast <- raster(ext=ext, resolution = res)
    
    # categorize the output in blocks with the dimensions res x res 
    ras <- merged.out[merged.out$ITV == bool,]
    
    categories <- seq(0.20, 0.8, res)
    
    groups1_res <- cut(ras$value1, categories)
    groups2_res <- cut(ras$value2, categories)
    groups_res <- paste(groups1_res, groups2_res)
    split_res_df <- split(ras, groups_res)
    
    # calculate the resource ratio for the grouped output
    pop_sample <- lapply(split_res_df, function(x){
      df <- data.frame(x, r_1 = 100 * sum(x$RESOURCESG[x$SPECIESID == 2]) / sum(x$RESOURCESG))
      return(df)
    })
    
    # convert to a raster 
    pop_sample <- do.call(rbind, pop_sample)
    
    coordinates(pop_sample)<- ~value1 + value2
    
    # generate rasters of r_ITV, r_noITV ... 
    assign(paste0("r_", ifelse(bool, "ITV", "noITV")), rasterize(pop_sample, rast, pop_sample$r_1, fun = mean))
    
  }
  
  # stack the output
  r_stack <- stack(r_ITV, r_noITV, r_ITV - r_noITV)
  df <- setNames(data.frame(rasterToPoints(r_stack)), c("competing", "focal", "ITV", "no ITV", "ITV-noITV"))
  
  
  df$species_diff_ITV <- abs(2 * (50 - df$ITV))  # total interspecific differences with ITV
  df$species_diff_noITV <- abs(2 * (50 - df$`no ITV`)) # total interspecific differences withou ITV
  df$species_diff_diff <- df$species_diff_ITV - df$species_diff_noITV # difference between total interspecific differences with and without ITV
  
  
  
  # Fig 5., A
  (Focal <- 
      ggplot(df, aes(x = competing, y = focal, fill = `no ITV`))+
      geom_raster()+
      scale_fill_gradientn("Resource ratio\nof focal species\nwithout ITV [%]" , colors = c(RColorBrewer::brewer.pal(7, "RdBu")), values = c(0, 0.15, 0.25, 0.5, 0.75, 0.85, 1))+
      theme_classic()+
      xlab(expression(paste(alpha, "-level of competing species")))+
      ylab(expression(paste(alpha, "-level of focal species")))+
      theme(legend.position = "top", text = element_text(size = 10))+
      scale_x_continuous(breaks = c(0.25, 0.5, 0.75))+
      scale_y_continuous(breaks = c(0.25, 0.5, 0.75)))
  
  # add labels for the plots of the intermediate and the heterogeneous resource distribution experiment
  if(sim.dfs$PATCHINESS[1] >  10){
    Focal <- 
      Focal + 
      geom_text(aes(x = 0.52, y = 0.3, label = "inferior"), size = 2.5, color = "white")+
      geom_text(aes(x = 0.315, y = 0.52, label = "superior"), size = 2.5, color = "white")
  }

  
  # Fig 5., B
  (FD <- 
      ggplot(df, aes(x = competing, y = focal, fill = species_diff_noITV))+
      theme_classic()+
      geom_raster()+
      scale_fill_gradientn("Interspecific differ-\nences in resource\nratio without ITV [%]", colours = rev(c(RColorBrewer::brewer.pal(7, "RdBu"))[1:4]), values = c(0, 0.8, 1))+
      xlab(expression(paste(alpha, "-level of competing species")))+
      ylab(expression(paste(alpha, "-level of focal species")))+
      theme(legend.position = "top", text = element_text(size = 10))+
      scale_x_continuous(breaks = c(0.25, 0.5, 0.75))+
      scale_y_continuous(breaks = c(0.25, 0.5, 0.75)))
  
  
  
  # Fig 5., C
  
  # select color breaks dependent on patchiness level
  if(sim.dfs$PATCHINESS[1] ==  10){
    color.breaks = c(-0.5, 0, 0.5)
  }else{
    if(sim.dfs$PATCHINESS[1] == 90)
    {
      color.breaks = c(-4, -2, -1, 0, 1)
    }else
    {
      color.breaks = c(-1, 0, 1)
    }
  }
  

  (FD_ITV <-
      ggplot(df, aes(x = competing, y = focal, fill = species_diff_diff))+
      geom_raster()+
      scale_fill_gradientn("Change in interspecific\ndifferences in res-\nource ratio with ITV [%]" , colors = c(RColorBrewer::brewer.pal(5, "RdBu")), values = c(0, .4, abs(min(df$species_diff_diff, na.rm = T))/sum(abs(range(df$species_diff_diff, na.rm = T))), .9, 1), breaks = color.breaks, guide = guide_colourbar(title.position = "left"))+
      theme_classic()+
      xlab(expression(paste("Mean ", alpha, "-level of competing species")))+
      ylab(expression(paste("Mean ", alpha, "-level of focal species")))+
      theme(legend.position = "top", text = element_text(size = 10))+
      scale_x_continuous(breaks = c(0.25, 0.5, 0.75))+
      scale_y_continuous(breaks = c(0.25, 0.5, 0.75)))
  
  
  # Fig 5., D
  
  # plot for non-random resource distribution
  if(sim.dfs$PATCHINESS[1] >  10){
    (point_plot <- ggplot(df, aes(x = species_diff_noITV, y = species_diff_diff))+
       geom_point()+ 
       geom_segment(x = 0, y = 0, xend = max(df$species_diff_noITV), yend = 0, color = "gray")+
       geom_label(aes(x = max(df$species_diff_noITV) - 5, y = 0, label = "no equalization"), label.size = NA)+
       geom_segment(x = 0, y = 0, xend = max(df$species_diff_noITV), yend = - max(df$species_diff_noITV), color = "gray")+
       geom_label(aes(x = abs(min(df$species_diff_diff)) - .3, y = min(df$species_diff_diff) + .3, label = "full\nequalization"), label.size = NA)+
       geom_smooth(method = "lm", se = F, color = "coral1")+
       theme_classic()+
       theme(text = element_text(size = 10))+
       xlab("Interspecific differences in\nresource ratio without ITV [%]")+
       ylab("Change in interspecific differ-\nences in resource ratio with ITV"))
  }
  
  # plot for random resource distribution
  if(sim.dfs$PATCHINESS[1] ==  10){
    df$grouping <- factor(df$focal < 0.55 & df$`no ITV`< 0.5 | df$competing < 0.55 & df$`no ITV` > 0.5, labels = c("too slow", "too fast"))
    
    (point_plot <- ggplot(df[df$focal < 0.5, ], aes(x = species_diff_noITV, y = species_diff_diff, group = as.factor(grouping)))+
        geom_point(aes(color = grouping))+
        scale_color_manual("", values = c("orange", "navyblue"))+
        geom_segment(x = 0, y = 0, xend = max(df$species_diff_noITV), yend = 0, color = "gray")+
        geom_label(aes(x = max(df$species_diff_noITV) - 2, y = 0, label = "no equalization"), label.size = NA)+
        geom_segment(x = 0, y = 0, xend = max(df$species_diff_noITV), yend = - max(df$species_diff_noITV), color = "gray")+
        geom_label(aes(x = abs(min(df$species_diff_diff)) - .1, y = min(df$species_diff_diff) + .1, label = "full\nequalization"), label.size = NA)+
        geom_smooth(method = "lm", se = F, color = "black")+
        theme_classic()+
        theme(legend.position = "top", text = element_text(size = 10))+
        xlab("Interspecific differences\nin resource ratio without ITV")+
        ylab("Change in interspecific differ-\nences in resource ratio with ITV"))
  }
  
  # combine all plots and save them to the figs folder
  (graphical_abstract <- point_plot + 
      xlab("Interspecific fitness differences without ITV")+
      ylab("Change in interspecific differences due to ITV"))
  
  
  all.together <- gridExtra::grid.arrange(ggdraw(Focal + coord_fixed()) + draw_plot_label(c("A")),
                                          ggdraw(FD + coord_fixed()) + draw_plot_label("B"),
                                          ggdraw(FD_ITV + coord_fixed()) + draw_plot_label("C"), 
                                          ggdraw(point_plot) + draw_plot_label("D"), 
                                          layout_matrix = matrix(c(1,2,
                                                                   3,4), byrow = T, nc = 2))
  
  
  ggsave(here("figs", paste0(sim.dfs$PATCHINESS[1], "_res-", res,  "_resource_ratio2.jpeg")), all.together, width = 16.6, height = 15, units = "cm")
  
  ggsave(here("figs", paste0(sim.dfs$PATCHINESS[1], "_res-", res,  "_resource_ratio2.pdf")), all.together, width = 16.6, height = 15, units = "cm")
  
  
  ggsave(here("figs", paste0(sim.dfs$PATCHINESS[1], "_res-", res,  "_graphical_abstract.png")), graphical_abstract, width = 5, height = 5)
  
}


