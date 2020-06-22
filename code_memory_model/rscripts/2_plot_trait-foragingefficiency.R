library(here)
library(ggplot2)
library(cowplot)
library(distances)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)

#################
### SUBPLOT B ###
#################


# specify the experiment's name and day of processing

n_inds <- 750
exp.name_ <- paste0("final_fitness_trait_relationship_", n_inds, "inds")
day_ <- "200514"

source(here("code_memory_model", "rscripts", "analysis_tools.R"))

# read simulation run of fitness-trait-relationship or create it if its not there yet 
FTR <- read.sim(day = day_, exp.name = exp.name_, nofcores = 4)


#convert the absolute foraging within 1,000 time steps to foraging efficiency per time step  
FTR$RESOURCESG <- FTR$RESOURCESG / 1000 


# read the trade-off data (number of clumps vs. resources per clump foraged)
gg.data <- readr::read_csv(here("simulations", day_, paste0("final_trade_off_", n_inds, "inds_all.csv")))


# trade-off data
gg.data$variable <- rep(c("Number of different \n clumps utilized [n]", "Mean resources \n per clump"), each = nrow(gg.data)/2)

# TFE relationship data
FTR.data <- cbind(FTR[, c("ALPHA", "PATCHINESS")], "Foraging\nefficiency [n / t]", FTR$RESOURCESG)


names(gg.data) <- names(FTR.data) <- c("ALPHA", "PATCHINESS", "OBSERVATION", "VALUE")


gg.df <- rbind(FTR.data, gg.data)

gg.df$PATCHINESS <- as.factor(gg.df$PATCHINESS)
names(gg.df)


# generate subplot B
gg_FTR_TO <- 
ggplot(gg.df, aes(x = ALPHA, y = VALUE, group = PATCHINESS, linetype = PATCHINESS))+
  theme_classic()+
  geom_smooth(data = gg.df %>% filter(OBSERVATION == "Foraging\nefficiency [n / t]"), method = "lm", formula = y ~ poly(x, 2), se = F, color = "black", size = 0.5)+
geom_smooth(data = gg.df %>% filter(OBSERVATION != "Foraging\nefficiency [n / t]"), method = "loess", formula = y ~ x, se = F, color = "black", size = 0.5)+
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"), "Resource distribution", labels = c("random", "intermediate", "patchy"))+
  ylab("")+
  xlab(expression(paste(alpha, "-level")))+
  scale_x_continuous(breaks = c(.25, .5, .75))+
  theme(legend.position = "bottom", legend.justification = "center", strip.background = element_blank(), strip.placement = "outside", legend.background = element_rect(fill = NA, color = "gray"))+
  facet_wrap(~OBSERVATION, scales = "free_y", nc = 1, strip.position = "left")+
  theme(legend.key.size = ggplot2::unit(4, "mm"))+
  xlim(c(0,1))+
  guides(linetype = guide_legend(title.position = "top"))

gg_FTR_TO
# save plots as .png and .pdf


#################
### SUBPLOT A ###
#################

### select and process some exemplary movement paths and the landscape for ### 

# movement paths
example.move <- read.csv(list.files(here("simulations", day_, "final_trade_off_750inds"), full.names = T, pattern = ".csv$")[1])

example.ls <- raster(list.files(here("simulations", day_, "final_trade_off_750inds"), full.names = T, pattern = ".png$")[1])

img <- imager::load.image(list.files(here("simulations", day_, "final_trade_off_750inds"), full.names = T, pattern = ".png$")[1])

# process image of landscape and convert to raster object
red <- imager::imresize(img, scale = 1/1.5)
imager::save.image(red, here("temp.png"))
example.ls <- raster("temp.png")

# select two different behaviour types 
example.move.ss <- example.move[abs(example.move$ALPHA -   0.65949367) < 0.01 | abs(example.move$ALPHA -  0.34113924) < 0.01, ]

# shift x coordinate for better visualization
example.move.ss$XCOR <- example.move.ss$XCOR + 20 

# convert to spatial object and run process.simulation function to correct for toroidal landscape
coordinates(example.move.ss) <- ~XCOR + YCOR
example.move.proc <- process.simulation(file = example.move.ss, output.dir = here(), save = F, landscape.size = 250, overlap.analysis = F)

# set extent of raster of landscape to extent of movement examples
extent(example.ls) <- extent(example.move.proc)

WHOs <- unique(example.move.proc$WHO)

# distinugish between slower and faster individual
example.ls.gg <- data.frame(rbind(rasterToPoints(example.ls), rasterToPoints(example.ls)), WHO = rep(c("slower", "faster"), each = ncell(example.ls)))

xy <- coordinates(example.move.proc[example.move.proc$WHO == WHOs[1],])

# process spatial object
gg.lines <- 
rbind(fortify(as(as(example.move.proc[example.move.proc$WHO == WHOs[15],], "SpatialLines"), "SpatialLinesDataFrame")) %>% mutate(id = runif(1)),
      fortify(as(as(example.move.proc[example.move.proc$WHO == WHOs[9],], "SpatialLines"), "SpatialLinesDataFrame")) %>% mutate(id = runif(1))) %>%
  mutate(WHO = rep(c("slower", "faster"), each = 1000))

# shift coordinates to improve visualization
gg.lines <- 
gg.lines %>% 
  group_by(id) %>% 
  mutate(long = long - min(long) + 25) %>% 
  mutate(lat = lat - min(lat) + 20)

# visualize movement paths
example.ls.gg <- data.frame(rbind(rasterToPoints(example.ls), rasterToPoints(example.ls)), WHO = rep(c("slower", "faster"), each = ncell(example.ls)))
example.ls.gg <-
example.ls.gg %>% 
  filter(x < (max(gg.lines$long) + 10) & y < (max(gg.lines$lat)+ 10))
(movement_paths <- 
ggplot()+
  geom_raster(data = example.ls.gg, aes(x = x, y = y, fill = as.factor(temp > 200)))+
  geom_path(data = gg.lines, aes(x = long, y = lat, group = group))+
  facet_wrap(~ WHO)+
  scale_fill_manual(values = c("cornsilk4", "white"))+
  theme_map()+
  coord_fixed()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  ggtitle("Behavioural Types")+
  theme(text = element_text(size = 9)))
  

# adjust text size of subplot B
  (gg_FTR_TO <-   
  gg_FTR_TO +
    theme(text = element_text(size = 9)))

#######################################
### combine subplot A and subplot B ###
#######################################

all.together <- 
gridExtra::grid.arrange(ggdraw(movement_paths) + draw_plot_label(c("A"), size = 11),
                        ggdraw(gg_FTR_TO) + draw_plot_label(c("B"), size = 11),
                        layout_matrix = matrix(c(1,1,2, 2,2,2,2), nc = 1))

widths = 7.5
heights = 17

ggsave(here("figs", paste0(n_inds, "inds_FTR_TO.png")), all.together, width = widths, height = heights, units = "cm", dpi = 600) # 8.2, 11.0, 17.3
ggsave(here("figs", paste0(n_inds, "inds_FTR_TO.tiff")), all.together, width = widths, height = heights, units = "cm", dpi = 600)  

ggsave(here("figs", paste0(n_inds, "inds_FTR_TO.pdf")), all.together, width = widths, height = heights, units = "cm", dpi = 600)  

