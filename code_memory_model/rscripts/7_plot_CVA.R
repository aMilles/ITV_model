# Rscript to generate Fig. 6 and supplements


library(here)
library(ggplot2)
library(cowplot)
library(distances)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)


source(here("code_memory_model", "rscripts", "analysis_tools.R"))

### Read simulation files ###
L1 <- "simulations"
L2 <- "200520"
L3 <- "Coviability_analysis_homogeneous_noher"

### iterate through different patchiness levels (heterogeneous  = patchy, homogeneous = random)
for(L3 in c("Coviability_analysis_heterogeneous_noher",
            "Coviability_analysis_intermediate_noher",
            "Coviability_analysis_homogeneous_noher")){
  path2sim <- here(L1, L2, L3)
  all.files <- list.files(path2sim, pattern = ".csv", full.names = T)
  alpha.files <- list.files(path2sim, pattern = "_alpha.csv", full.names = T)
  cov.files <- all.files[!all.files %in% alpha.files]
  
  # assign simulations to separater lists
  all_sims <- lapply(as.list(cov.files), function(x) readr::read_csv(x) %>% mutate(ID = runif(1, max = 100)) %>% mutate(SEED = paste0(unique(SPECMEAN), collapse = "_")))
  
  # combine to one data.frame
  sims <- do.call(rbind, all_sims)
  
  # get species-1-mean and species-2-mean
  sims <- 
    sims %>%
    group_by(ID) %>%
    mutate(s2 = unique(SPECMEAN)[length(unique(SPECMEAN))]) %>%
    mutate(s1 = unique(SPECMEAN)[1]) %>% 
    filter(DATETIME > 0)
  
  # for each combination of species-1-mean and species-2-mean and ITV/no ITV, calculate the mean time to extinction of one of the two species
  sims <- sims %>% 
    group_by(s2, s1, ITV) %>% 
    mutate(time2ext = mean(DATETIME)) %>% 
    filter(!duplicated(s1) & !duplicated(s2))
  
  # sims with ITV
  sims.ITV <- sims %>% 
    filter(ITV > 0)
  
  # sims without ITV
  sims.no.ITV <- sims %>% 
    filter(ITV == 0)
  
  # have s1 and s2 in sims without and with ITV in the same order
  sims.ITV <- 
    sims.ITV[match(paste0(sims.no.ITV$s1, sims.no.ITV$s2), paste0(sims.ITV$s1, sims.ITV$s2)),]
  
  # generate one dataset with the times to extinction with and without ITV in separate columns
  sims.all <- sims.ITV %>% 
    mutate(time2ext.ITV = time2ext)
  
  sims.all$time2ext.noITV = sims.no.ITV$time2ext  
  
  # select the times and mean trait levels and change the layout to a ggplot-readable format
  sims.all <-
    sims.all %>% 
    mutate(time2ext.diff =  time2ext.ITV - time2ext.noITV) %>% 
    dplyr::select(s1,s2, time2ext.noITV, time2ext.ITV, time2ext.diff)
  
  
  m.sims.all <- 
    reshape2::melt(sims.all, id.vars = c("ITV", "s1", "s2"))
  
  #################
  ### SUBPLOT A ###
  #################
  # plot mean time to extinction without ITV
  noITV <- 
    ggplot(m.sims.all %>% filter(variable == "time2ext.noITV"), aes(x = s1, y = s2, fill = value))+
    geom_raster()+
    coord_fixed()+
    theme_clean()+
    scale_fill_viridis_c()
  
  
  (no.ITV <- 
      ggplot(m.sims.all %>% filter(variable == "time2ext.noITV"), aes(x = s1, y = s2, fill = value))+
      geom_raster()+
      coord_fixed()+
      theme_classic()+
      scale_fill_gradientn("Mean time\nto extinction\nwithout ITV [%]" , colors = c(RColorBrewer::brewer.pal(7, "Blues")), values = c(0, 0.15, 0.25, 0.5, 0.75, 0.85, 1), breaks = c(20000, 180000))+
      xlab(expression(paste(alpha, "-level of species A")))+
      ylab(expression(paste(alpha, "-level of species B")))+
      theme(text = element_text(size = 10), legend.position = "top"))
  
  
  #################
  ### SUBPLOT B ###
  #################
  # differences due to ITV
  (diff <- 
      ggplot(m.sims.all %>% filter(variable == "time2ext.diff"), aes(x = s1, y = s2, fill = value))+
      geom_raster()+
      coord_fixed()+
      theme_classic()+
      scale_fill_gradientn("Change in mean\n time to extinction\n due to ITV [%]" , colors = c(RColorBrewer::brewer.pal(7, "RdBu")), values = c(0, 0.15, 0.2, 0.35, 0.5, 0.85, 1), breaks = c(-70000, 0, 150000))+
      #scale_fill_gradientn("Change in mean\n time to extinction\n due to ITV [%]" , colors = c(RColorBrewer::brewer.pal(5, "RdBu")), values = c(0, 0.15, 0.47, 0.85, 1), breaks = c(-150000, 0, 150000))+
      xlab(expression(paste(alpha, "-level of species A")))+
      ylab(expression(paste(alpha, "-level of species B")))+
      theme(text = element_text(size = 10), legend.position = "top"))
  
  # adjust visualization for different patchiness levels
  if(grep("intermediate", L3) > 0){
    (diff <- diff +
    scale_fill_gradientn("Change in mean\n time to extinction\n due to ITV [%]" , colors = c(RColorBrewer::brewer.pal(5, "RdBu")), values = c(0, 0.25, 0.45, 0.75, 1), breaks = c(-150000, 0, 150000)))
  } 
  
  
  if(grep("homogeneous", L3) > 0){
    (diff <- diff +
      scale_fill_gradientn("Change in mean\n time to extinction\n due to ITV [%]" , colors = c(RColorBrewer::brewer.pal(5, "RdBu")), values = c(0, 0.001, 0.01, .5, 1), breaks = c(0, 100000)))
  } 
  
  # histogram of differences (not included in MS)
  (gg_hist <- ggplot(m.sims.all %>% filter(variable == "time2ext.diff"), aes(x = value))+
      geom_histogram(fill = "gray50", bins = 20)+
      xlab("Change in mean time to\nextinction due to ITV")+
      ylab("Count")+
      geom_vline(xintercept = 0, linetype = "dashed")+
      theme_classic()+
      theme(text = element_text(size = 10)))
  
  ####################
  ### SAVE FIGURES ###
  ####################
  
  all.together <- 
    gridExtra::grid.arrange(ggdraw(no.ITV)+draw_plot_label("A"), ggdraw(diff)+draw_plot_label("B"), nrow = 2)
  
  
  all.together_with_hist <- 
    gridExtra::grid.arrange(ggdraw(no.ITV)+draw_plot_label("A"), ggdraw(diff)+draw_plot_label("B"), ggdraw(gg_hist)+draw_plot_label("C"), nrow = 3)
  
  ggsave(here("figs", paste0(L2, L3, "_CVA_newplot.png")), all.together, height = 18, width = 8, units = "cm", dpi = 600)
  ggsave(here("figs", paste0(L2, L3, "_CVA_newplot.pdf")), all.together, height = 18, width = 8, units = "cm", dpi = 600)
  
  ggsave(here("figs", paste0(L2, L3, "_CVA_newplot_withhist.png")), all.together_with_hist, height = 18, width = 10, units = "cm", dpi = 600)
  
  ggsave(here("figs", paste0(L2, L3, "_CVA_newhist.png")), gg_hist, height = 8.5, width = 8, units = "cm", dpi = 600)
  
}



