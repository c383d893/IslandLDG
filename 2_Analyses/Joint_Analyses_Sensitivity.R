############################
###### LOAD PACKAGES #######
############################

library(mgcv); library(gridExtra); library(betareg); library(MASS); library(lme4); library(lmerTest); library(lsmeans); library(ggeffects); library(spdep); library(ggplot2); library(ncf); library(ape); library(sjPlot); library(gridExtra); library(MuMIn); library(tidyverse); library(maps); library(sf); library(tidyverse); library(relaimpo)
options(na.action = "na.fail")

############################
###### LOAD FUNCTIONS ######
############################

#overdispersion function
Check.disp <- function(mod,dat) {
  N <- nrow(dat)
  p <- length(coef(mod))
  E1 <- resid(mod, type = "pearson")
  Dispersion <- sum(E1^2)/ (N-p)
  return(Dispersion)
}

#RAC function
Spat.cor <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

#RAC function when locations repeat (shift latlon)
Spat.cor.rep <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude) + matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

############################
######## READ DATA #########
############################

debt <- readRDS("data/debt_native_sensitivity_myc_latitude_data_2023.RDS")

mycdat <- readRDS("data/GAMexp_native_sensitivity_myc_latitude_data_2023.RDS") %>% select("entity_ID","propAM_exp", "propEM_exp", "propORC_exp") %>% rename(AM.ml = propAM_exp, EM.ml = propEM_exp, ORC.ml = propORC_exp)
nfixdat <- readRDS("data/GAMexp_native_sensitivity_nfix_latitude_data_2023.RDS") %>% select("entity_ID","propnfix_exp") %>% rename(nfix.ml = propnfix_exp)
polldat <- readRDS("data/GAMexp_native_sensitivity_poll_latitude_data_2023.RDS") %>% select("entity_ID","propb_exp") %>% rename(poll.ml = propb_exp)

# join these:
dat <- debt %>%
  full_join(mycdat, by = "entity_ID") %>%
  full_join(nfixdat, by = "entity_ID") %>%
  full_join(polldat, by = "entity_ID") %>%
  #make sum biotic drivers
  #choose one
  #more conservative: average % of species that exhibit one of these syndromes
  mutate(biotic.ml = (nfix.ml + poll.ml + AM.ml)/3) %>%
  #probability that one of 3 independent events occurred: https://www.omnicalculator.com/statistics/probability-three-events
  #mutate(biotic.ml = nfix.ml + poll.ml + AM.ml - nfix.ml * poll.ml - nfix.ml * AM.ml - poll.ml * AM.ml + nfix.ml * poll.ml * AM.ml) %>%
  #area, dist, elev_range already scaled and log transformed
  mutate(abslatitude = as.vector(scale(abslatitude)), nfix.ml = as.vector(scale(log10((nfix.ml)+.01))),poll.ml = as.vector(scale(log10((poll.ml)+.01))),
       AM.ml = as.vector(scale(log10((AM.ml)+.01))), EM.ml = as.vector(scale(log10((EM.ml)+.01))), ORC.ml = as.vector(scale(log10((ORC.ml)+.01))),
       biotic.ml = as.vector(scale(log10((biotic.ml)+.01)))) %>%
  #mutate(abslatitude - as.vector(scale(abslatitude)), nfix.ml = as.vector(scale(logit(nfix.ml))),poll.ml = as.vector(scale(logit(poll.ml))),
  #     AM.ml = as.vector(scale(logit(AM.ml))), EM.ml = as.vector(scale(logit(EM.ml))), ORC.ml = as.vector(scale(logit(ORC.ml)))) 
  drop_na() %>%
  mutate(sprichdiff = ifelse(sprichdiff <0,0, sprichdiff))

dat.br <- dat %>% filter(Tdebt > 0)

############################
######## VAR CORR ##########
############################

corr.vars <- dat %>%
  select(c("abslatitude","dist","area","elev_range","temp","prec","poll.ml", "nfix.ml","AM.ml","EM.ml", "ORC.ml"))                              
corr.vars <- dat %>%
  select(c("abslatitude","dist","area","elev_range","prec","biotic.ml"))                              
corr.mat <- as.matrix(cor(corr.vars))               

############################
##### SPRICH DIFF MODS #####
############################

# no interactions
sprichdiff.mod <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml, data = dat) 
rac <- Spat.cor.rep(sprichdiff.mod, dat, 2000)
sprichdiff.mod.rac <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml + rac, data = dat) 
summary(sprichdiff.mod.rac)
calc.relimp(sprichdiff.mod.rac)
calc.relimp(sprichdiff.mod.rac, rela=TRUE)

# check assumptions
par(mfrow = c(2,2))
plot(sprichdiff.mod.rac)
vif(sprichdiff.mod.rac)

# create models with all possible combinations of predictors
dredged_models <- dredge(sprichdiff.mod.rac)

# perform AIC model averaging, only considering models with AIC values within 4 of the best model
averaged_model_results <- model.avg(dredged_models, delta < 4, fit = TRUE)

# examine results
summary(averaged_model_results)

#############################
###### MODEL SELECTION ######
#############################

# add biotic.ml*area
M1 <- glm(sprichdiff ~ biotic.ml*area + dist +  prec +  abslatitude, data = dat) 
rac <- Spat.cor.rep(M1, dat, 2000)
M1.rac <- glm(sprichdiff ~ biotic.ml*area + dist + prec + abslatitude + rac, data = dat) 
summary(M1.rac)
AIC(sprichdiff.mod.rac, M1.rac) #int improves fit

# add biotic.ml*dist
M2 <- glm(sprichdiff ~  biotic.ml*area + biotic.ml*dist + prec + abslatitude , data = dat) 
rac <- Spat.cor.rep(M2, dat, 2000)
M2.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + abslatitude + rac, data = dat) 
summary(M2.rac)
AIC(sprichdiff.mod.rac, M2.rac) # int improves fit

# add biotic.ml*prec
M3 <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + biotic.ml*prec + abslatitude, data = dat) 
rac <- Spat.cor.rep(M3, dat, 2000)
M3.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + biotic.ml*prec + abslatitude + rac, data = dat) 
summary(M3.rac)
AIC(M2.rac, M3.rac) # int DOES NOT improve fit

# add biotic.ml*abslatitude
M3.2 <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + biotic.ml*abslatitude, data = dat) 
rac <- Spat.cor.rep(M3.2, dat, 2000)
M3.2.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + biotic.ml*abslatitude + rac, data = dat) 
summary(M3.2.rac)
AIC(M2.rac, M3.2.rac)  # int improves fit

# add abslatitude*area
M3.3 <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + abslatitude*area, data = dat) 
rac <- Spat.cor.rep(M3.3, dat, 2000)
M3.3.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + abslatitude*area + rac, data = dat) 
summary(M3.3.rac)
AIC(M2.rac, M3.3.rac) # int DOES NOT improve fit

# add abslatitude*dist
M3.4 <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + abslatitude*dist, data = dat) 
rac <- Spat.cor.rep(M3.4, dat, 2000)
M3.4.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + prec + abslatitude*dist + rac, data = dat) 
summary(M3.4.rac)
AIC(M2.rac, M3.4.rac) # int DOES NOT improve fit

# add abslatitude*prec
M3.5 <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + abslatitude*prec, data = dat) 
rac <- Spat.cor.rep(M3.5, dat, 2000)
M3.5.rac <- glm(sprichdiff ~ biotic.ml*area + biotic.ml*dist + abslatitude*prec + rac, data = dat) 
summary(M3.5.rac)
AIC(M2.rac, M3.5.rac) # int DOES NOT improve fit

calc.relimp(M3.2.rac)
calc.relimp(M3.2.rac, rela =TRUE)

############################
########### PLOT ###########
#### RELATIVE IMPORTANCE ###
############################

mod.dat <- calc.relimp(M3.2.rac)$lmg %>% as.data.frame() %>%
  rownames_to_column(., var = "variable") %>%
  filter(! variable=="rac") %>%
  mutate(variable = case_when(variable=="abslatitude" ~ "Absolute latitude",
                               variable=="area" ~ "Area", 
                               variable=="biotic.ml" ~ "Mutualism filter strength",
                               variable=="dist" ~ "Distance",
                               variable=="prec" ~ "Precipitation",
                               variable=="elev_range" ~ "Elevation range")) %>%
  filter(!is.na(variable))

colnames(mod.dat)<-c("variable", "rel.imp")

mod.dat.ordered <- mod.dat %>% mutate(variable = fct_reorder(variable, rel.imp))

forest.plot <-
  ggplot(data=mod.dat.ordered, aes(x=variable, y=rel.imp, ymin=rel.imp, ymax=rel.imp)) +
  geom_pointrange(alpha = 0.8, size = 2, color="coral3") + 
  #geom_hline(yintercept=0, lty=2, color='darkgrey') +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab(" ") + ylab("Variance explained") +
  theme_classic(base_size = 40) +
  theme(legend.position = 'none') +
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(angle = 45, size = 30)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) 

png("figures/sensitivity_sprichdiff_relimp_forestplot.jpg", width = 11, height = 10, units = 'in', res = 300)
forest.plot
dev.off()

############################
########### PLOT ###########
##### MODEL ESTIMATES ######
############################

mod.dat <- summary(averaged_model_results)$coefmat.full %>% as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  filter(! variable=="(Intercept)") %>%
  filter(! variable=="rac") %>%
  mutate(variable = case_when(variable=="abslatitude" ~ "Absolute latitude",
                              variable=="area" ~ "Area", 
                              variable=="biotic.ml" ~ "Mutualism filter strength",
                              variable=="dist" ~ "Distance",
                              variable=="prec" ~ "Precipitation",
                              variable=="elev_range" ~ "Elevation range")) %>%
  filter(!is.na(variable)) 

colnames(mod.dat)<-c("est", "std.err","adjusted.se","zval","pval", "variable")

mod.dat.ordered <- mod.dat %>% mutate(variable = fct_reorder(variable, est))

forest.plot <-
  ggplot(data=mod.dat.ordered, aes(x=variable, y=est, ymin=est-std.err, ymax=est+std.err), fill = "coral3") +
  geom_pointrange(alpha = 0.8, size = 2, color="coral3") + 
  geom_hline(yintercept=0, lty=2, color='darkgrey') +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab(" ") + ylab("Effect of 1 SD of predictor \non species deficit") +
  theme_classic(base_size = 40) +
  theme(legend.position = 'none') +
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(angle = 45, size = 30)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) 

png("figures/sensitivity_sprichdiff_dredgeest_forestplot.jpg", width = 11, height = 10, units = 'in', res = 300)
forest.plot
dev.off()