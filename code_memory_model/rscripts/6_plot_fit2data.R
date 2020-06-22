library(here)
library(ggplot2)
library(cowplot)
library(distances)
library(parallel)
library(foreach)
library(doSNOW)

source(here("code_memory_model",  "rscripts",  "analysis_tools.R"))


######################
##### SOME FUNS ######
######################
# scale values to a range of - 1 to 1
scale.personality <- function(x) (2* (x - min(x))/diff(range(x))) - 1

#scale values to relative proportions of the maximum, the maximum may be the maximum of the vector or the maximum of some vector (e.g. observed / max(simulated))
scale.data <- function(x, maxx = NA) if(is.na(maxx)){x/max(x)}else{x/maxx}

#NRMSE
RMSE <- function(x, y) return(sqrt(mean((x - y)^2))/mean(y))

#Rsq
Rsq <- function(x, y) return(cor(x, y )^2)

##### SOME GENERAL DATASETS ####
predict.seq <- data.frame(boldness.scaled = seq(-1, 1, length.out = 50))
predict.seq.alpha <- data.frame(ALPHA.scaled = seq(-1, 1, length.out = 50))



##### REAL DATA #####
# read data from Schirmer et al. 2019 (can't be provided here)

real.area <- xlsx::read.xlsx(here("data", "SpaceUse_UM16.xlsx"), sheetIndex = 2,endRow = 37, stringsAsFactors = F)
real.distmoved <- xlsx::read.xlsx(here("data", "MovementTracks_Uckermark2016.xlsx"), sheetName = "allTracks",endRow = 37, stringsAsFactors = F)
real.ol <- xlsx::read.xlsx(here("data", "SpaceUse_UM16.xlsx"), sheetIndex = 3,endRow = 37, stringsAsFactors = F)

plot(real.ol$intraspecific_OL ~ real.ol$Boldness)

# merge kernel areas, overlaps and track lengths
real.merged <- Reduce(function(x, y) merge(x = x,y = y, by = c("PIT", "Species", "Grid"), all = T), 
       list(real.area, real.distmoved, real.ol))
plot(intraspecific_OL ~ Boldness, real.merged[real.merged$Species == "SFM",])

real.merged$Boldness.x[is.na(real.merged$Boldness.x)] <- real.merged$Boldness.y[is.na(real.merged$Boldness.x)]

real.merged <- real.merged[real.merged$Species == "BV",]

# choose columns of interest, remove one individual that has an extreme boldness and hence can hardly be used for general pattern analysis. This individual was pregnant and may, hence, have should some other biological signal. 

real.data <- setNames(real.merged[, c("Boldness.x", "Kernel50", "Kernel95", "Track_length", "intraspecific_OL")], c("Boldness", "SIZE50", "SIZE95", "DISTMOVED", "HR50"))

real.data <- real.data[real.data$Boldness < 10, ]

rm(list = c("all.merged", "real.area", "real.ol", "real.distmoved", "real.merged"))

real.data$boldness.scaled <- scale.personality(real.data$Boldness)


# create predictions from simple linear models
real.model.scaled <- real.model <- data.frame("Boldness" = predict.seq,
                         "SIZE50" = predict(lm(SIZE50 ~ boldness.scaled + I(boldness.scaled^2), real.data), newdata = predict.seq),
                         "SIZE95" = predict(lm(SIZE95 ~ boldness.scaled + I(boldness.scaled^2), real.data), newdata = predict.seq),
                         "HR50" =  predict(lm(HR50 ~ boldness.scaled, real.data), newdata = predict.seq),
                         "DISTMOVED" = predict(lm(DISTMOVED ~ boldness.scaled, real.data), newdata = predict.seq))

real.data.scaled <- real.data

for(metric in names(real.model)[-1]){
  real.data.scaled[, metric] <- scale.data(real.data[,metric], maxx = max(real.model[,metric]))
  
  if(metric != "HR50") real.model.scaled[, metric] <- scale.data(real.model[,metric])
  if(metric == "HR50"){
    real.model.scaled[, metric] <- real.model[,metric]
    real.data.scaled[, metric] <- real.data[,metric]
  } 
}



mtext(side = 1, text = "Boldness Scores")
names(real.data)
plot.sc <- real.data.scaled
plot.nsc <- real.data
plot.sc$Boldness <- plot.sc$boldness.scaled
names(plot.sc)[1:5] <- names(plot.nsc)[1:5] <- c("Boldness", "Home Range Size (Kernel Density; 50 %)", "Home Range Size (Kernel Density; 95 %)", "Total Distance Moved", "Mean Home Range Overlap (Kernel Density; 50 %)")

mdf <- reshape2::melt(data.frame(rbind(plot.nsc[,-6], plot.sc[,-6]), scaled = rep(c(F, T), each = nrow(plot.nsc))), id.vars = c("Boldness", "scaled"))

new.levels <- c("Home range area\n(50 % kernel)", "Home range area\n(95 % kernel)", "Total distance moved", "Home range\noverlap (95 % kernel)")
levels(mdf$variable) <- new.levels


scaled <- ggplot(mdf[mdf$scaled,], aes(x = Boldness, y = value))+
  geom_point()+
  facet_wrap( ~ variable, scales = "free", nc = 1)+
 ggtitle("Scaled")+
  theme_classic()+
  ylab("Home range metrics")

not.scaled <- ggplot(mdf[!mdf$scaled,], aes(x = Boldness, y = value))+
  geom_point()+
  facet_wrap( ~ variable, scales = "free", nc = 1)+
  ggtitle("Not scaled")+
  theme_classic()+
  ylab("")

scaled.vs.notscaled<- gridExtra::grid.arrange(not.scaled, scaled, ncol = 2)

ggsave(here("figs", "nonscaled_vs_scaled_real_data.png"), scaled.vs.notscaled, width = 6, height = 8)

png(here("figs", "scaled_data_linear_model.png"))
par(mfrow = c(2,2))
plot(SIZE50 ~ boldness.scaled, real.data.scaled, ylab = "Home Range Size (Kernel Density; 50 %)", xlab = "")
lines(SIZE50 ~ boldness.scaled, real.model.scaled)
plot(SIZE95 ~ boldness.scaled, real.data.scaled, ylab = "Home Range Size (Kernel Density; 95 %)", xlab = "")
lines(SIZE95 ~ boldness.scaled, real.model.scaled)
plot(DISTMOVED ~ boldness.scaled, real.data.scaled, ylab = "Total Distance Moved", xlab = "")
lines(DISTMOVED ~ boldness.scaled, real.model.scaled)
plot(HR50 ~ boldness.scaled, real.data.scaled, ylab = "Mean Home Range Overlap (Kernel Density; 50 %)", xlab = "")
lines(HR50 ~ boldness.scaled, real.model.scaled)
dev.off()

########################
#### SIMULATED DATA ####
########################

#read simulated data from LHS (n = 1,000)
ols <- read.sim(day = "200514", exp.name = "HRanalysis_processed_1000_fit", nofcores = 1)

# split simulation runs into single list elements
split.ols <- split(ols, ols$ID)

# predict the fit of simulated data using the same model as for observed data
real.boldness.scaled <- data.frame(ALPHA.scaled = real.data.scaled$boldness.scaled)

sim.data.scaled <- do.call(rbind, 
                           lapply(split.ols, function(x){
    x$ALPHA.scaled <- -scale.personality(x$ALPHA)
    
    SIZE50.model <- lm(SIZE50 ~ ALPHA.scaled + I(ALPHA.scaled^2), x)
    SIZE50.pred <- scale.data(predict(SIZE50.model, newdata = predict.seq.alpha))
    SIZE95.model <- lm(SIZE95 ~ ALPHA.scaled + I(ALPHA.scaled^2), x)
    SIZE95.pred <- scale.data(predict(SIZE95.model, newdata = predict.seq.alpha))
    DISTMOVED.model <- lm(DISTMOVED ~ ALPHA.scaled, x)
    DISTMOVED.pred <- scale.data(predict(DISTMOVED.model, newdata = predict.seq.alpha))
    HR50.model <- lm(HR95 ~ ALPHA.scaled, x)
    HR50.pred <- predict(HR50.model, newdata = predict.seq.alpha)                   
    
    pred <- data.frame(SIZE50 = SIZE50.pred,
                       SIZE50.RMSE = RMSE(SIZE50.pred, real.model.scaled$SIZE50),
                       SIZE95 = SIZE95.pred,
                       SIZE95.RMSE = RMSE(SIZE95.pred, real.model.scaled$SIZE95),
                       DISTMOVED = DISTMOVED.pred,
                       DISTMOVED.RMSE = RMSE(DISTMOVED.pred, real.model.scaled$DISTMOVED),
                       HR50 = HR50.pred,
                       HR50.RMSE = RMSE(HR50.pred, real.model.scaled$HR50),
                       PATCHINESS = x$patchiness[1],
                       N_INDS = x$n_inds[1],
                       RESOURCE_COVER = x$resource_cover[1],
                       TRAIT_MEAN = x$trait_mean[1],
                       ITV = x$ITV[1],
                       ALPHA.scaled = predict.seq.alpha)
    return(pred)
    } 
  )
)

################################
##### Analyze fit to data ######
################################

names(sim.data.scaled)
summary(sim.data.scaled)
# get the parameter values that are lower than 90 % of RMSE values for each metric
ranks <- setNames(data.frame(apply(sim.data.scaled[,c(2,4,6, 8)], 2, function(x) x)), c("SIZE50.rank", "SIZE95.rank", "DISTMOVED.rank", "HR50.rank"))

sim <- data.frame(sim.data.scaled, data.frame(ranks, rank = apply(ranks, 1, sum)))

############################################################
##### plot the parameter distribution / correlations for the top 25 % #####
############################################################

gg.sim.params <- sim[sim$rank < quantile(sim$rank, .25),]
gg.sim.params <- reshape2::melt(gg.sim.params[,c("N_INDS", "ITV", "TRAIT_MEAN", "PATCHINESS", "RESOURCE_COVER")])


#correlation between parameters
png(here("figs", "parameter_cor.png"))
par(oma = c(0,2,0,0))
psych::cor.plot(gg.sim.params[,9:13], diag = F, stars = F, scale = T, numbers = T, cex = .8)
dev.off()



# parameter distribution of all paramteters
parameter_dist <- ggplot(gg.sim.params, aes(x = value, fill = variable))+
  geom_histogram(bins = 10, fill = "lightgray", color = "gray")+
  facet_wrap(~as.factor(variable), scales = "free")+
  xlab("Paramter value")+
  ylab("Count")+
  theme_classic()

ggsave(here("figs", "parameter_dist.png"), width = 6.9, height = 5.5)
ggsave(here("figs", "parameter_dist.pdf"), width = 6.9, height = 5.5)

# parameter distribution of ITV
ITV_dist <- ggplot(gg.sim.params[gg.sim.params$variable == "ITV",], aes(x = value, fill = variable))+
  geom_density(fill = "lightgray", color = "gray")+
  scale_fill_brewer(palette = "Accent")+
  xlab("ITV")+
  theme(text = element_text(size = 12))+
  geom_vline(xintercept = mean(gg.sim.params[gg.sim.params$variable == "ITV",]$value), linetype = "dashed")

ggplot(gg.sim.params[gg.sim.params$variable == "ITV",], aes(x = value, fill = variable))+
  geom_histogram(fill = "lightgray", color = "gray")+
  scale_fill_brewer(palette = "Accent")+
  xlab("ITV")+
  theme(text = element_text(size = 12))+
  geom_vline(xintercept = median(gg.sim.params[gg.sim.params$variable == "ITV",]$value), linetype = "dashed")

ggsave(here("figs", "ITV_dist.png"), ITV_dist, width = 2, height = 2)
ggsave(here("figs", "ITV_dist.pdf"), ITV_dist, width = 2, height = 2)

###########################################################################
##### format the 4 metrics and the error to a ggplot-readable layout #####
###########################################################################

# convert layout of simulated data
sim$ID <- rep(seq(NROW(sim) / NROW(predict.seq)), each = NROW(predict.seq))
names(sim)
gg.metrics <- reshape2::melt(sim[,c("SIZE50", "SIZE95", "DISTMOVED", "HR50", "ALPHA.scaled", "ID", "ITV", "rank")], id.vars = c("ID", "ALPHA.scaled", "ITV", "rank"), variable.name = "metrics")
gg.error <- reshape2::melt(sim[,c("SIZE50.RMSE", "SIZE95.RMSE", "DISTMOVED.RMSE", "HR50.RMSE", "ID")], id.vars = c("ID"), variable.name = "error.metric", value.name = "RMSE")


# add column which informs about whether parameter set is within the 25-% quantile
gg.merge <- cbind(gg.metrics, gg.error[,-1])
gg.merge$q <- "all"
gg.merge$q[quantile(gg.merge$rank, .25) > gg.merge$rank] <- "25 %"
gg.merge$rank <- NULL


# convert layout of empirical data
gg.real <- reshape2::melt(real.data.scaled[, -1], id.vars = "boldness.scaled", variable.name = "metrics")
new.levels <- c("Home range area (50 % kernel)", "Home range area (95 % kernel)", "Total distance moved", "Home range overlap (95 % kernel)")
levels(gg.merge$metrics) <- levels(gg.real$metrics) <- new.levels

################################################
#### visualize fit of the top 25 % to data ####
################################################
all.fitted <- 
ggplot()+  
  geom_line(data = gg.merge[gg.merge$q != "all",], aes(x = ALPHA.scaled, y = value, group = ID, color = q), alpha = 1)+
  geom_smooth(data = gg.real[gg.real$metrics %in% new.levels[1:2], ], aes(x = boldness.scaled, y = value), alpha = .1, color = "black", se = F, method = "lm", formula = y ~ x + I(x^2))+
  geom_smooth(data = gg.real[gg.real$metrics %in% new.levels[3:4], ], aes(x = boldness.scaled, y = value), alpha = .1, color = "black", se = F, method = "lm", formula = y ~ x)+
  geom_smooth(data = gg.real[gg.real$metrics %in% new.levels[1:2], ], aes(x = boldness.scaled, y = value), alpha = .5, color = "black",  se = T, method = "lm", formula = y ~ x + I(x^2), fill = "cornsilk2")+
  geom_smooth(data = gg.real[gg.real$metrics %in% new.levels[3:4], ], aes(x = boldness.scaled, y = value), alpha = .5, color = "black", se = T, method = "lm", formula = y ~ x, fill = "cornsilk2")+
  geom_point(data = gg.real, aes(x = boldness.scaled, y = value), color = "black")+
  scale_color_manual(values = c("cadetblue", "cadetblue4"))+
  facet_wrap(~metrics)+
  theme_classic()+
  xlab(expression(paste("Scaled ", alpha, " / boldness score")))+
  ylab(expression(paste("Scaled home range metrics")))+
  theme(text = element_text(size = 10), legend.position = "none")


####################
### SAVE FIGURES ###
####################

widths = 16.6
heights = 10

ggsave(here("figs", "all.fitted.png") ,all.fitted, width = widths, height = heights, units = "cm", dpi = 600) # 8.2, 11.0, 17.3
ggsave(here("figs", "all.fitted.pdf") ,all.fitted, width = widths, height = heights, units = "cm", dpi = 600)  
ggsave(here("figs", "all.fitted.tiff") ,all.fitted, width = widths, height = heights, units = "cm", dpi = 600) 


