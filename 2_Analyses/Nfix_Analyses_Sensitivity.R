############################
###### LOAD PACKAGES #######
############################

library(mgcv); 
library(gridExtra); 
library(betareg); 
library(MASS); 
library(lme4); 
library(lmerTest); 
library(lsmeans); 
library(ggeffects); 
library(spdep); 
library(ggplot2); 
library(ncf); 
library(ape); 
library(sjPlot); 
library(gridExtra); 
library(MuMIn);
library(maps); 
library(sf); 
library(car);
library(viridis);
library(tidyverse)

options(na.action = "na.fail")

packageVersion(c("mgcv")) #‘1.8.41’
packageVersion(c("gridExtra")) # ‘2.3’
packageVersion(c("betareg")) # ‘3.1.4’
packageVersion(c("MASS")) # ‘7.3.58.1’
packageVersion(c("lme4")) # ‘1.1.31’
packageVersion(c("lmerTest")) # ‘3.1.3’
packageVersion(c("lsmeans")) # ‘2.30.0’
packageVersion(c("ggeffects")) # ‘1.1.4’
packageVersion(c("spdep")) # ‘1.2.7’
packageVersion(c("ggplot2")) #'3.4.0'
packageVersion(c("ncf")) # ‘1.3.2’
packageVersion(c("ape")) # ‘5.6.2’
packageVersion(c("sjPlot")) # ‘2.8.12’
packageVersion(c("gridExtra")) # '2.3'
packageVersion(c("MuMIn")) # ‘1.47.1’
packageVersion(c("tidyverse")) #  ‘1.3.2’
packageVersion(c("maps")) # ‘3.4.1’
packageVersion(c("sf")) # ‘1.0.9’
packageVersion(c("car")) # ‘3.1.1’
packageVersion(c("viridis")) # ‘0.6.2’
packageVersion(c("tidyverse")) # ‘1.3.2’
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

dat <-  readRDS("data/native_nfix_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", "elev_range",
           "nfix", "nonfix")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  select(-geology) %>%          
  filter(area > 6) %>% # based on paper Patrick shared
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>% 
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>%
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), 
         temp = as.vector(scale(temp)), 
         prec = as.vector(scale(log10((prec)+.01))),
         elev_range = as.vector(scale(log10((elev_range)+.01)))) %>%
  filter(!entity_class2 == "Non-oceanic") %>%
  drop_na() %>%
  # take 20% of nfix and move to nonfix
  mutate(nfix.switch = nfix*.20) %>%
  mutate(nfix = nfix - nfix.switch) %>%
  mutate(nonfix = nonfix + nfix.switch) %>%
  select(-nfix.switch)

dat2 <-  readRDS("data/native_nfix_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist",   "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec",    
           "nfix", "nonfix")) %>%
  rename(temp = CHELSA_annual_mean_Temp, prec = CHELSA_annual_Prec) %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland")) %>%
  select(-geology) %>%                                                              
  mutate(abslatitude = abs(latitude)) %>%                                             
  mutate(abslatitude = as.vector(abslatitude)) %>%    
  mutate(elev_range = ifelse(elev_range==0,1, elev_range)) %>%                                   
  mutate(elev_range = ifelse(is.na(elev_range),1, elev_range)) %>%
  filter(area > 6) %>% # based on paper Patrick shared
  #for models only; remove for figs:
  mutate(area = as.vector(scale(log10((area)+.01))), dist = as.vector(scale(log10((dist)+.01))), elev_range = as.vector(scale(log10((elev_range)+.01))),
         temp = as.vector(scale(temp)), prec = as.vector(scale(log10((prec)+.01)))) %>%
  #filter(!entity_class2 == "Non-oceanic") %>%
  drop_na() %>%
  # take 20% of nfix and move to nonfix
  mutate(nfix.switch = nfix*.20) %>%
  mutate(nfix = nfix - nfix.switch) %>%
  mutate(nonfix = nonfix + nfix.switch) %>%
  select(-nfix.switch)

############################
######## ML PREDICT ########
############################

############################
######### MAINLAND #########
############################

dat.ml <- dat %>% filter(entity_class2=="Mainland")

gam.mod_sprich <- gam(sprich ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod_sprich)

mod <- gam.mod_sprich
new.dat.sprich <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred_sprich <- predict.gam(mod,newdata = new.dat.sprich, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.sprich$abslatitude) 

# check assumptions
gam.check(gam.mod_sprich)
Check.disp(gam.mod_sprich,dat)

######################################
############ MAINLAND NFIX ###########
######################################

gam.mod.nfix <- gam(nfix ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.nfix)

mod <- gam.mod.nfix
new.dat.nfix <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000)))
pred.nfix <- predict.gam(mod,newdata = new.dat.nfix, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.nfix$abslatitude, nfix = "nfix") 

# check assumptions
gam.check(gam.mod.nfix)
Check.disp(gam.mod.nfix,dat)

###############################
####### MAINLAND NONFIX #######
###############################

gam.mod.nonfix <- gam(nonfix ~ s(abslatitude) , family=nb(link="log"), data = dat.ml) 
summary(gam.mod.nonfix)

mod <- gam.mod.nonfix
new.dat.nonfix <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000)))
pred.nonfix <- predict.gam(mod, newdata = new.dat.nonfix, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.nonfix$abslatitude, nfix ="nonfix") 

# check assumptions
gam.check(gam.mod.nonfix)
Check.disp(gam.mod.nonfix,dat)

############################
########## PLOT ############
##### LAT BY NFIX TYPE #####
############################

pred.mainland <- rbind(pred.nfix,pred.nonfix)

dat.ml.cond <- dat.ml %>%
  select(abslatitude, nfix, nonfix) %>%
  gather(key = "nfix", value = "sprich", nfix, nonfix)

# create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3", "darkgrey"))
fillScale <- scale_fill_manual(values=c("darkseagreen3", "darkgrey"))

pred.mainland$nfix <- ordered(pred.mainland$nfix, levels = c("nfix", "nonfix"))
dat.ml.cond$nfix <- ordered(dat.ml.cond$nfix, levels = c("nfix", "nonfix"))

# rename nfix for legend:

pred.mainland <- pred.mainland %>% mutate(nfix = case_when(nfix=="nfix" ~ "N-fixing",
                                                           nfix=="nonfix" ~ "Non N-fixing"))

dat.ml.cond <- dat.ml.cond %>% mutate(nfix = case_when(nfix=="nfix" ~ "N-fixing",
                                                           nfix=="nonfix" ~ "Non N-fixing"))

lat.nfixtype <-
  ggplot(pred.mainland, aes(x = abslatitude, y = fit, color = nfix, fill = nfix))+
  geom_line(size = 1) +
  geom_point(data = dat.ml.cond, aes(x = abslatitude, y = sprich, color = factor(nfix)), alpha = 0.3, size = 4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,700)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="N-fixing \nType",override.aes=list(fill = NA,size = 3, alpha = 0.7,linetype = c(0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/Sensitivity_Nfix_Latbynfixtype.jpg", width=10, height= 10, units='in', res=300)
lat.nfixtype
dev.off()

############################
###### PREDICT EXP IS ######
####### ALL ISLANDS ########
############################

dat.is.min <- dat %>% filter(entity_class2=="Oceanic") %>% select(c('entity_ID',"abslatitude"))

pred_sprich <- predict.gam(gam.mod_sprich,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_sprich_df <- cbind(dat.is.min,pred_sprich) %>% rename(sprich_exp = fit)

pred_nfix <- predict.gam(gam.mod.nfix,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_nfix_df <- cbind(dat.is.min,pred_nfix) %>% rename(nfix_exp = fit)

pred_nonfix <- predict.gam(gam.mod.nonfix,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_nonfix_df <- cbind(dat.is.min,pred_nonfix) %>% rename(nonfix_exp = fit)

pred_is.dat <- dat %>%
  filter(entity_class2=="Oceanic") %>%
  left_join(pred_sprich_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_nfix_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_nonfix_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('nfix_exp','nonfix_exp','sprich_exp'), as.integer) %>%
  mutate(propnfix_exp = nfix_exp/(nfix_exp + nonfix_exp))

# write out
saveRDS(pred_is.dat,"data/GAMexp_native_sensitivity_nfix_latitude_data_2023.RDS")

############################
####### MAIN MODELS ########
############################

############################
######## CREATE DATA #######
############################

# read data and calc debts
full.dat <- dat2 %>%
  select(c("entity_ID", "dist", "area")) 

# version 1; debt & C_debt neg to 0, >1 to 1:
pred_is.dat.shrink <- readRDS("data/GAMexp_native_sensitivity_nfix_latitude_data_2023.RDS") %>%
  select(-area) %>%
  mutate(sprich = nfix + nonfix, sprich_exp = nfix_exp + nonfix_exp) %>%
  mutate(nfix.diff = nfix_exp - nfix, nonfix.diff = nonfix_exp - nonfix, sprichdiff = sprich_exp - sprich) %>%
  mutate(nfix.debt = (nfix.diff/nfix_exp), nonfix.debt = (nonfix.diff/nonfix_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(nfix.debt = ifelse(nfix.debt < 0, 0, nfix.debt)) %>% mutate(nonfix.debt = ifelse(nonfix.debt <0, 0, nonfix.debt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(nfix.debt = ifelse(nfix.debt > 1, 1, nfix.debt)) %>% mutate(nonfix.debt = ifelse(nonfix.debt >1, 1, nonfix.debt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_nfix.debt = (nfix.diff/sprichdiff), C_nonfix.debt = (nonfix.diff/sprichdiff)) %>% 
  mutate(C_nfix.debt = ifelse(C_nfix.debt < 0, 0, C_nfix.debt)) %>% mutate(C_nonfix.debt = ifelse(C_nonfix.debt <0, 0, C_nonfix.debt)) %>% 
  mutate(C_nfix.debt = ifelse(C_nfix.debt > 1, 1, C_nfix.debt)) %>% mutate(C_nonfix.debt = ifelse(C_nonfix.debt >1, 1, C_nonfix.debt)) %>% 
  left_join(full.dat, by = "entity_ID") 

############################
######### PREP DATA ########
############################

pred_is.dat <- pred_is.dat.shrink

pred_is.dat.alld <- pred_is.dat %>% 
  select(c('entity_ID','nfix.debt','nonfix.debt')) %>%
  gather(key = "nfix", value = "debt", nfix.debt, nonfix.debt) %>%
  mutate(nfix = case_when(nfix == "nfix.debt" ~ "nfix",                                   
                          nfix == "nonfix.debt" ~ "nonfix")) 

pred_is.dat.allcd <- pred_is.dat %>% 
  select(c('entity_ID','C_nfix.debt','C_nonfix.debt')) %>%
  gather(key = "nfix", value = "debt.c", C_nfix.debt, C_nonfix.debt)%>%
  mutate(nfix = case_when(nfix == "C_nfix.debt" ~ "nfix",                                   
                          nfix == "C_nonfix.debt" ~ "nonfix"))

pred_is.dat.all.diff <- pred_is.dat %>% 
  select(c('entity_ID','sprich','latitude','longitude','abslatitude','nfix.diff','nonfix.diff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "nfix", value = "diff", nfix.diff, nonfix.diff) %>%
  mutate(nfix = case_when(nfix == "nfix.diff" ~ "nfix",                                   
                          nfix == "nonfix.diff" ~ "nonfix")) 

pred_is.dat.all.exp <- pred_is.dat %>% 
  select(c('entity_ID','nfix_exp','nonfix_exp')) %>%
  gather(key = "nfix", value = "exp", nfix_exp, nonfix_exp) %>%
  mutate(nfix = case_when(nfix == "nfix_exp" ~ "nfix",                                   
                          nfix == "nonfix_exp" ~ "nonfix")) 

pred_is.dat.all.obs <- pred_is.dat %>% 
  select(c('entity_ID', 'nfix', 'nonfix')) %>%
  gather(key = "nfix", value = "obs", nfix, nonfix) 

pred_is.dat.all.T <- pred_is.dat %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("nfix", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(nfix = as.factor(nfix))

pred_is.dat.all <- within(pred_is.dat.all, nfix <- relevel(nfix, ref = "nonfix"))

############################
######### ONE MODEL ########
############################

############################
########## DEBT C ##########
############################

pred.allcatcd <- glm(debt.c ~ poly(abslatitude,3,raw = TRUE)*nfix + area + dist +elev_range+ prec+ (1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
rac <- Spat.cor.rep(pred.allcatcd,pred_is.dat.all,2000)
pred.allcatcd.rac  <- glm(debt.c ~ poly(abslatitude,3,raw = TRUE)*nfix + area + dist +elev_range+ prec+ rac +(1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
summary(pred.allcatcd.rac)

# glm diagnostics:
par(mfrow=c(3,2))
# homogenetity of variance:
plot(fitted(pred.allcatcd.rac) ~ resid(pred.allcatcd.rac, type = "pearson"))
# independence:
plot(pred_is.dat.all$area ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$dist ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$elev_range ~ resid(pred.allcatcd.rac, type = "pearson"))
boxplot(resid(pred.allcatcd.rac, type = "pearson") ~ pred_is.dat.all$nfix)
# do not look for normality:
# outliers:
plot(cooks.distance(pred.allcatcd.rac), type ="h")
# check dispersion:
Check.disp(pred.allcatcd.rac, pred_is.dat.all)

############################
########### PLOT ###########
#### C DEBT:FULL MODEL #####
############################

new.dat.nfix <- pred_is.dat.all %>% filter(nfix=="nfix") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.nfix <- predict.glm(pred.allcatcd.rac,newdata = new.dat.nfix, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.nfix$abslatitude) 

new.dat.nonfix <- pred_is.dat.all %>% filter(nfix=="nonfix") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.nonfix <- predict.glm(pred.allcatcd.rac,newdata = new.dat.nonfix, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.nonfix$abslatitude)  

allcatcd.plot <- 
  ggplot() +
  #Nfix section
  geom_line(data = pred.nfix, mapping = aes(x = abslatitude, y = fit), color ="darkseagreen3")+
  geom_ribbon(data = pred.nfix, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkseagreen3", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(nfix=="nfix"), aes(x = abslatitude, y=debt.c),color ="darkseagreen3", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(nfix=="nfix"),aes(x = cuts_labeled, y = mean_y),color = "darkseagreen3",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(nfix=="nfix"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkseagreen3", width=3,size=1.5) +
  #no nfix section
  geom_line(data = pred.nonfix, mapping = aes(x = abslatitude, y = fit), color ="darkgrey")+
  geom_ribbon(data = pred.nonfix, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkgrey", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(nfix=="nonfix"), aes(x = abslatitude, y=debt.c),color ="darkgrey", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(nfix=="nonfix"),aes(x = cuts_labeled, y = mean_y),color = "darkgrey",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(nfix=="nonfix"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkgrey", width=3,size=1.5) +
  theme_classic(base_size = 15) +
  #geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylab("Contribution deficit") +
  xlab("Absolute latitude") +
  ylim(0,1.2)

# write out
png("figures/Sensitivity_Nfix_LatPoly_contdebt_fullmodel_shrink.jpg", width = 6, height = 6, units ='in', res = 300)
allcatcd.plot
dev.off()