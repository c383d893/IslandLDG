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

dat <-  readRDS("data/native_myc_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area",  "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec","elev_range",
           "AM","EM","ORC","NM")) %>%
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
  #filter(!entity_class2 == "Non-oceanic") %>% 
  drop_na() %>%
  # take 20% of AM and move to NM
  mutate(myc.switch = AM*.20) %>%
  mutate(AM = AM - myc.switch) %>%
  mutate(NM = NM + myc.switch) %>%
  select(-myc.switch)

dat2 <-  readRDS("data/native_myc_latitude_data_2023.RDS") %>%
  filter(!entity_class == "undetermined") %>%                                        
  select(c("entity_ID","entity_class","sprich","latitude","longitude","geology", "area", "elev_range", "dist", "CHELSA_annual_mean_Temp", "CHELSA_annual_Prec", 
           "AM","EM","ORC","NM")) %>%
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
  # take 20% of AM and move to NM
  mutate(myc.switch = AM*.20) %>%
  mutate(AM = AM - myc.switch) %>%
  mutate(NM = NM + myc.switch) %>%
  select(-myc.switch)

############################
######## ML PREDICT ########
############################

############################
######### MAINLAND #########
############################

dat.ml <- dat %>% filter(entity_class2 == "Mainland")

gam.mod_sprich <- gam(sprich ~ s(abslatitude) , family = nb(link = "log"), data = dat.ml) 
summary(gam.mod_sprich)

mod <- gam.mod_sprich
new.dat.sprich <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred_sprich <- predict.gam(mod,newdata = new.dat.sprich, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.sprich$abslatitude)

# check assumptions
gam.check(gam.mod_sprich)
Check.disp(gam.mod_sprich,dat)

############################
####### MAINLAND AM ########
############################

gam.mod.AM <- gam(AM ~ s(abslatitude),family=nb(link="log"), data = dat.ml) 
summary(gam.mod.AM)

mod <- gam.mod.AM
new.dat.AM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.AM <- predict.gam(mod,newdata = new.dat.AM, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.AM$abslatitude, myctype = "AM") 

# check assumptions
gam.check(gam.mod.AM)
Check.disp(gam.mod.AM,dat)

############################
####### MAINLAND EM ########
############################

gam.mod.EM <- gam(EM ~ s(abslatitude),family=nb(link="log"), data = dat.ml) 
summary(gam.mod.EM)

mod <- gam.mod.EM
new.dat.EM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.EM <- predict.gam(mod,newdata = new.dat.EM, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.EM$abslatitude, myctype = "EM") 

# check assumptions
gam.check(gam.mod.EM)
Check.disp(gam.mod.EM,dat)

############################
####### MAINLAND ORC #######
############################

gam.mod.ORC <- gam(ORC ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.ORC)

mod <- gam.mod.ORC
new.dat.ORC <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.ORC <- predict.gam(mod,newdata = new.dat.ORC, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(abslatitude = new.dat.ORC$abslatitude, myctype = "ORC") 

# check assumptions
gam.check(gam.mod.ORC)
Check.disp(gam.mod.ORC,dat)

############################
####### MAINLAND NM ########
############################

gam.mod.NM <- gam(NM ~ s(abslatitude), family=nb(link="log"), data = dat.ml) 
summary(gam.mod.NM)

mod <- gam.mod.NM
new.dat.NM <- with(mod$model, expand.grid(abslatitude = seq(min(abslatitude), max(abslatitude), length = 1000))) 
pred.NM <- predict.gam(mod,newdata = new.dat.NM, type = "response", se = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude = new.dat.NM$abslatitude, myctype = "NM") 

# check assumptions
gam.check(gam.mod.NM)
Check.disp(gam.mod.NM,dat)

############################
########## PLOT ############
##### LAT BY MYC TYPE ######
############################

#pred.mainland <- rbind(pred.AM,pred.EM,pred.ORC,pred.NM)
pred.mainland <- rbind(pred.AM,pred.EM,pred.NM)

#dat.ml.cond <- dat.ml %>%
#  select(abslatitude, AM,EM,ORC,NM) %>%
#  gather(key="myctype", value="sprich", AM,EM,ORC,NM) 

dat.ml.cond <- dat.ml %>%
  select(abslatitude, AM,EM,NM) %>%
  gather(key="myctype", value="sprich", AM,EM,NM) 

# create a custom color scale
colScale <- scale_colour_manual(values=c("royalblue4", "royalblue3","darkgrey"))
fillScale <- scale_fill_manual(values=c("royalblue4", "royalblue3","darkgrey"))

#pred.mainland$myctype <- ordered(pred.mainland$myctype, levels = c("AM","EM","ORC","NM"))
#dat.ml.cond$myctype <- ordered(dat.ml.cond$myctype, levels = c("AM","EM","ORC","NM"))

pred.mainland$myctype <- ordered(pred.mainland$myctype, levels = c("AM","EM","NM"))
dat.ml.cond$myctype <- ordered(dat.ml.cond$myctype, levels = c("AM","EM","NM"))

lat.myctype <-
ggplot(pred.mainland,aes(x =abslatitude,y=fit,color =myctype, fill =myctype))+
  geom_line(size=1) +
  geom_point(data=dat.ml.cond, aes(x = abslatitude,y=sprich, color =factor(myctype)), alpha=0.3,size=4)+ 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.5) + 
  xlab("Absolute latitude") +
  ylab("Species richness")+
  theme_classic(base_size = 25)+
  ylim(0,4000)+
  colScale+
  fillScale+
  theme(legend.position = 'none')+
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  guides(color = guide_legend(title="Mycorrhizal \nType",override.aes=list(fill =NA,size =3, alpha =0.7,linetype = c(0, 0, 0))))+
  theme(axis.text.x = element_text(size =20),axis.text.y = element_text(angle = 45,size=20))

# write out
png("figures/Sensitivity_Myc_LatbyMycType.jpg", width=10, height= 10, units='in', res=300)
lat.myctype
dev.off()

############################
###### PREDICT EXP IS ######
####### ALL ISLANDS ########
############################

dat.is.min <- dat %>% filter(entity_class2=="Oceanic") %>% select(c('entity_ID',"abslatitude"))

pred_sprich <- predict.gam(gam.mod_sprich, newdata = dat.is.min, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_sprich_df <- cbind(dat.is.min,pred_sprich) %>% rename(sprich_exp = fit)

pred_AM <- predict.gam(gam.mod.AM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_AM_df <- cbind(dat.is.min,pred_AM) %>% rename(AM_exp = fit)

pred_EM <- predict.gam(gam.mod.EM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_EM_df <- cbind(dat.is.min,pred_EM) %>% rename(EM_exp = fit)

pred_ORC <- predict.gam(gam.mod.ORC,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_ORC_df <- cbind(dat.is.min,pred_ORC) %>% rename(ORC_exp = fit)

pred_NM <- predict.gam(gam.mod.NM,newdata = dat.is.min,type = "response", se = TRUE) %>%
  as.data.frame() %>%
  select(fit)
pred_NM_df <- cbind(dat.is.min,pred_NM) %>% rename(NM_exp = fit)

pred_is.dat <- dat %>%
  filter(entity_class2=="Oceanic") %>%
  left_join(pred_sprich_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_AM_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_EM_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_ORC_df, by= c('entity_ID','abslatitude')) %>%
  left_join(pred_NM_df, by= c('entity_ID','abslatitude')) %>%
  mutate_at(c('AM_exp','EM_exp','ORC_exp', 'NM_exp','sprich_exp'), as.integer) %>%
  mutate(propAM_exp = AM_exp/(AM_exp + EM_exp + ORC_exp + NM_exp)) %>%
  mutate(propEM_exp = EM_exp/(AM_exp + EM_exp + ORC_exp + NM_exp)) %>%
  mutate(propORC_exp = ORC_exp/(AM_exp + EM_exp + ORC_exp + NM_exp))

# write out
saveRDS(pred_is.dat,"data/GAMexp_native_sensitivity_myc_latitude_data_2023.RDS")

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
pred_is.dat.shrink <- readRDS("data/GAMexp_native_sensitivity_myc_latitude_data_2023.RDS") %>%
  select(-area) %>%
  mutate(sprich = AM + EM + ORC + NM, sprich_exp = AM_exp + EM_exp + ORC_exp + NM_exp) %>%
  mutate(AMdiff = AM_exp - AM, EMdiff = EM_exp - EM, ORCdiff = ORC_exp - ORC, NMdiff = NM_exp - NM, sprichdiff = sprich_exp - sprich) %>%
  mutate(AMdebt = (AMdiff/AM_exp), EMdebt = (EMdiff/EM_exp), ORCdebt = (ORCdiff/ORC_exp), NMdebt = (NMdiff/NM_exp), Tdebt = (sprichdiff/sprich_exp))  %>%
  mutate(AMdebt = ifelse(AMdebt < 0, 0, AMdebt)) %>% mutate(EMdebt = ifelse(EMdebt <0, 0, EMdebt)) %>% mutate(ORCdebt = ifelse(ORCdebt <0, 0, ORCdebt)) %>% mutate(NMdebt = ifelse(NMdebt <0, 0, NMdebt)) %>% mutate(Tdebt = ifelse(Tdebt <0, 0, Tdebt)) %>%
  mutate(AMdebt = ifelse(AMdebt > 1, 1, AMdebt)) %>% mutate(EMdebt = ifelse(EMdebt >1, 1, EMdebt)) %>% mutate(ORCdebt = ifelse(ORCdebt >1, 1, ORCdebt)) %>% mutate(NMdebt = ifelse(NMdebt >1, 1, NMdebt)) %>% mutate(Tdebt = ifelse(Tdebt >1, 1, Tdebt)) %>%
  mutate(C_AMdebt = (AMdiff/sprichdiff), C_EMdebt = (EMdiff/sprichdiff), C_ORCdebt = (ORCdiff/sprichdiff), C_NMdebt = (NMdiff/sprichdiff)) %>% 
  mutate(C_AMdebt = ifelse(C_AMdebt < 0, 0, C_AMdebt)) %>% mutate(C_EMdebt = ifelse(C_EMdebt <0, 0, C_EMdebt)) %>% mutate(C_ORCdebt = ifelse(C_ORCdebt <0, 0, C_ORCdebt)) %>% mutate(C_NMdebt = ifelse(C_NMdebt <0, 0, C_NMdebt)) %>%
  mutate(C_AMdebt = ifelse(C_AMdebt > 1, 1, C_AMdebt)) %>% mutate(C_EMdebt = ifelse(C_EMdebt >1, 1, C_EMdebt)) %>% mutate(C_ORCdebt = ifelse(C_ORCdebt >1, 1, C_ORCdebt)) %>% mutate(C_NMdebt = ifelse(C_NMdebt >1, 1, C_NMdebt)) %>%
  left_join(full.dat, by = "entity_ID") 

############################
######### PREP DATA ########
############################

pred_is.dat <- pred_is.dat.shrink

# write out
debt <- pred_is.dat.shrink %>% select(c("entity_ID","Tdebt","sprichdiff","dist","area","elev_range","latitude","longitude","abslatitude","sprich_exp", "temp","prec"))
saveRDS(debt,"data/debt_native_sensitivity_myc_latitude_data_2023.RDS")

pred_is.dat.alld <- pred_is.dat %>% 
  select(c('entity_ID','AMdebt','EMdebt','ORCdebt','NMdebt')) %>%
  gather(key = "myctype", value = "debt", AMdebt, EMdebt, ORCdebt, NMdebt) %>%
  mutate(myctype = case_when(myctype=="AMdebt" ~ "AM",                                   
                             myctype=="EMdebt" ~ "EM", 
                             myctype=="ORCdebt" ~ "ORC", 
                             myctype=="NMdebt" ~ "NM")) 

pred_is.dat.allcd <- pred_is.dat %>% 
  select(c('entity_ID','C_AMdebt','C_EMdebt','C_ORCdebt','C_NMdebt')) %>%
  gather(key = "myctype", value = "debt.c", C_AMdebt, C_EMdebt, C_ORCdebt, C_NMdebt)%>%
  mutate(myctype = case_when(myctype == "C_AMdebt" ~ "AM",                                   
                             myctype == "C_EMdebt" ~ "EM", 
                             myctype == "C_ORCdebt" ~ "ORC", 
                             myctype == "C_NMdebt" ~ "NM"))

pred_is.dat.all.diff <- pred_is.dat %>% 
  select(c('entity_ID','sprich','latitude','longitude','abslatitude','AMdiff','EMdiff','ORCdiff','NMdiff','dist','area','elev_range', 'temp','prec')) %>%
  gather(key = "myctype", value = "diff", AMdiff, EMdiff, ORCdiff, NMdiff) %>%
  mutate(myctype = case_when(myctype == "AMdiff" ~ "AM",                                   
                             myctype == "EMdiff" ~ "EM", 
                             myctype == "ORCdiff" ~ "ORC", 
                             myctype == "NMdiff" ~ "NM")) 

pred_is.dat.all.exp <- pred_is.dat %>% 
  select(c('entity_ID','AM_exp','EM_exp','ORC_exp','NM_exp')) %>%
  gather(key = "myctype", value = "exp", AM_exp, EM_exp, ORC_exp,NM_exp) %>%
  mutate(myctype = case_when(myctype == "AM_exp" ~ "AM",                                   
                             myctype == "EM_exp" ~ "EM", 
                             myctype == "ORC_exp" ~ "ORC", 
                             myctype == "NM_exp" ~ "NM")) 


pred_is.dat.all.obs <- pred_is.dat %>% 
  select(c('entity_ID', 'AM', 'EM', 'ORC', 'NM')) %>%
  gather(key = "myctype", value = "obs", AM, EM, ORC, NM) 

pred_is.dat.all.T <- pred_is.dat %>% select(entity_ID, sprichdiff)

pred_is.dat.all <- pred_is.dat.all.diff %>%
  left_join(pred_is.dat.all.exp, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.all.obs, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.alld, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.allcd, by = c("myctype", "entity_ID")) %>%
  left_join(pred_is.dat.all.T, by = c("entity_ID")) %>%
  mutate(debt.weights = exp, debt.c.weights = abs(sprichdiff)) %>%
  mutate(myctype = as.factor(myctype))

pred_is.dat.all <- within(pred_is.dat.all, myctype <- relevel(myctype, ref = "NM"))

############################
######### ONE MODEL ########
############################

############################
########## DEBT C ##########
############################

pred.allcatcd <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*myctype + area + dist + elev_range + prec + (1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
rac <- Spat.cor.rep(pred.allcatcd,pred_is.dat.all,2000)
pred.allcatcd.rac  <- glm(debt.c ~ poly(abslatitude,2,raw = TRUE)*myctype + area + dist + elev_range  + prec + rac +(1|entity_ID),weights = debt.c.weights, data = pred_is.dat.all) 
summary(pred.allcatcd.rac)

# glm diagnostics:
par(mfrow=c(3,2))
# homogenetity of variance:
plot(fitted(pred.allcatcd.rac) ~ resid(pred.allcatcd.rac, type = "pearson"))
# independence:
plot(pred_is.dat.all$area ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$dist ~ resid(pred.allcatcd.rac, type = "pearson"))
plot(pred_is.dat.all$elev_range ~ resid(pred.allcatcd.rac, type = "pearson"))
boxplot(resid(pred.allcatcd.rac, type = "pearson") ~ pred_is.dat.all$myctype)
# do not look for normality:
# outliers:
plot(cooks.distance(pred.allcatcd.rac), type ="h")
# check dispersion:
Check.disp(pred.allcatcd.rac, pred_is.dat.all)

############################
########### PLOT ###########
#### C DEBT:FULL MODEL #####
############################

new.dat.AM <- pred_is.dat.all %>% filter(myctype=="AM") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.AM <- predict.glm(pred.allcatcd.rac,newdata = new.dat.AM, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.AM$abslatitude) 

new.dat.EM <- pred_is.dat.all %>% filter(myctype=="EM") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.EM <- predict.glm(pred.allcatcd.rac,newdata = new.dat.EM, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.EM$abslatitude) 

new.dat.ORC <- pred_is.dat.all %>% filter(myctype=="ORC") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.ORC <- predict.glm(pred.allcatcd.rac,newdata = new.dat.ORC, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.ORC$abslatitude) 

new.dat.NM <- pred_is.dat.all %>% filter(myctype=="NM") %>%
  mutate(area = mean(pred.allcatcd.rac$model$area), dist = mean(pred.allcatcd.rac$model$dist),
         elev_range = mean(pred.allcatcd.rac$model$elev_range), prec = mean(pred.allcatcd.rac$model$prec), rac = mean(pred.allcatcd.rac$model$rac))
pred.NM <- predict.glm(pred.allcatcd.rac,newdata = new.dat.NM, type = "response", se = TRUE, newdata.guaranteed = TRUE) %>%
  as.data.frame() %>% 
  mutate(abslatitude=new.dat.NM$abslatitude) 

allcatcd.plot <- 
ggplot() +
  #AM section
  geom_line(data = pred.AM, mapping = aes(x = abslatitude, y = fit), color ="royalblue4")+
  geom_ribbon(data = pred.AM, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "royalblue4", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(myctype=="AM"), aes(x = abslatitude, y=debt.c),color ="royalblue4", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(myctype=="AM"),aes(x = cuts_labeled, y = mean_y),color = "royalblue4",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(myctype=="AM"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "royalblue4", width=3,size=1.5) +
  #EM section
  geom_line(data = pred.EM, mapping = aes(x = abslatitude, y = fit), color ="royalblue3")+
  geom_ribbon(data = pred.EM, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "royalblue3", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(myctype=="EM"), aes(x = abslatitude, y=debt.c),color ="royalblue3", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(myctype=="EM"),aes(x = cuts_labeled, y = mean_y),color = "royalblue3",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(myctype=="EM"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "royalblue3", width=3,size=1.5) +
  #ORC section
  #geom_line(data = pred.ORC, mapping = aes(x = abslatitude, y = fit), color ="darkorchid3")+
  #geom_ribbon(data = pred.ORC, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkorchid3", alpha= 0.5) +
  #geom_point(data = pred_is.dat.all %>% filter(myctype=="ORC"), aes(x = abslatitude, y=debt.c),color ="darkorchid3", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(myctype=="ORC"),aes(x = cuts_labeled, y = mean_y),color = "darkorchid3",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(myctype=="ORC"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkorchid3", width=3,size=1.5) +
  #NM section
  geom_line(data = pred.NM, mapping = aes(x = abslatitude, y = fit), color ="darkgrey")+
  geom_ribbon(data = pred.NM, aes(x = abslatitude, ymin=fit-se.fit, ymax = fit + se.fit), fill = "darkgrey", alpha= 0.5) +
  geom_point(data = pred_is.dat.all %>% filter(myctype=="NM"), aes(x = abslatitude, y=debt.c),color ="darkgrey", alpha= 0.3, size=2) +
  #geom_point(data = debtc.aggregated_data %>% filter(myctype=="NM"),aes(x = cuts_labeled, y = mean_y),color = "darkgrey",size = 1,shape = 15) +
  #geom_errorbar(data = debtc.aggregated_data %>% filter(myctype=="NM"),aes(x = cuts_labeled, y = mean_y, ymin = low95CI_y,ymax =high95CI_y),color = "darkgrey", width=3,size=1.5) +
  theme_classic(base_size = 15) +
  #geom_abline(intercept = 0, slope = 0, linetype="dashed")+
  ylab("Contribution deficit") +
  xlab("Absolute latitude") +
  ylim(0,1.2)

# write out
png("figures/Sensitivity_Myc_LatPoly_contdebt_fullmodel_shrink.jpg", width = 6, height = 6, units ='in', res = 300)
allcatcd.plot
dev.off()