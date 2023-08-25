# This script outputs counts of each pollination syndrome type (biotic v. abiotic)
# First based on species assignments
# Next based on family level

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)
packageVersion("tidyverse") # ‘1.3.2’

###################
#### LOAD DATA ####
###################

# bring in poll dat species 
pdat.s <- read.csv("data/gift_poll_disp_traits.csv", header = TRUE) %>%                     
  rename(species = work_species, poll.s = trait_value_3.6.1) %>%   
  select(species, poll.s) %>%
  drop_na()

pdat.f <- read.csv("data/gift_poll_disp_traits_tax.csv", header = TRUE) %>%
  filter(trait_ID == "3.6.1") %>%
  rename(poll.f = trait_value, family = taxon_name) %>%   
  select(family, poll.f) %>%
  drop_na() 

# bring in GIFT data:
# RESTRICTED
load("data/GIFT_LatGrad_restricted.RData")
species.dat <- checklists
geo.dat <- env

# EXTENDED (optional)
#load("data/GIFT_LatGrad_extended.RData")
#species.dat <- checklists
#geo.dat <- readRDS("data/geo.dat.extended_foranalyses.rds")

# bring in species dat
species.geo <- species.dat[[1]] %>% as.data.frame() %>%
  select(c('entity_ID', 'entity_class')) %>%
  mutate_at(c('entity_class'), as.character) %>%
  mutate_at(c('entity_ID'), as.integer) 

# bring in geo dat
geo <- geo.dat %>% 
  rename(GDP_satellite = "mean_GDP_satellite", CHELSA_annual_Prec = "mean_wc2.0_bio_30s_12", CHELSA_annual_mean_Temp = "mean_wc2.0_bio_30s_01", 
         Popdensity = "mean_gpw-v4_popdensity_UN_2015", Human_Footprint = "mean_HFP2009_unproj", Human_Influence = "mean_hii_v2geo", 
         elev_range = "range", age_Ma = "age_MA") %>%
  select(c('entity_ID', 'geo_entity', 'area','longitude','latitude', 'biome','dist', 'CHELSA_annual_mean_Temp', 'CHELSA_annual_Prec',
           'GDP_satellite', 'Popdensity', 'Human_Footprint', 'Human_Influence', 'elev_range', 'geology', 'age_Ma')) %>%
  mutate_at(c('area','longitude','latitude', 'dist','CHELSA_annual_mean_Temp', 'CHELSA_annual_Prec',
              'GDP_satellite', 'Popdensity','Human_Footprint', 'Human_Influence', 'elev_range', 'age_Ma'), as.numeric) %>%
  mutate_at(c('geo_entity','biome','geology'), as.character) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  left_join(species.geo, by = "entity_ID") %>%
  distinct(entity_ID, .keep_all = TRUE) %>%
  mutate(entity_class = case_when(entity_class == "Island" ~ "Island",
                                  entity_class == "Island Part" ~ "Island",          
                                  entity_class == "Island Group" ~ "Island",
                                  entity_class == "Mainland" ~ "Mainland",
                                  entity_class == "Island/Mainland" ~ "undetermined")) %>%
  mutate(geology = case_when(geology == "atoll" ~ "dev", 
                             geology == "floor" ~ "dev", 
                             geology == "floor/volcanic" ~ "dev",
                             geology == "volcanic" ~ "dev",
                             geology == "shelf" ~ "nondev",
                             geology == "fragment" ~ "nondev",
                             geology == "atoll/shelf/volcanic" ~ "nondev",
                             geology == "fragment/shelf/volcanic" ~ "nondev")) 

######################################
########### SP RICHNESS ##############
######################################

# with species checklist data, keep native, assign poll preferentially using species, then family
# join with geo data and keep unique sp by loc (entity_ID)
species <- species.dat[[2]] %>% as.data.frame() %>% 
  mutate(species = paste(genus, species_epithet, sep = " ")) %>% 
  select(c("entity_ID", "family", "native", "naturalized", "species", "name_ID", "genus")) %>%
  mutate_at(c('family', 'species', 'name_ID', 'genus'), as.character) %>%
  mutate_at(c('native', 'naturalized'), as.numeric) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  filter(native == 1) %>%                                                                     
  select(-c("native", "naturalized")) %>%                                                  
  left_join(pdat.s, by = c("species")) %>% 
  left_join(pdat.f, by = c("family"))  %>%                                     
  mutate(poll = ifelse(is.na(poll.s), poll.f, poll.s)) %>%    
  mutate(poll_level = ifelse(!is.na(poll.s), "species", "family")) %>%
  drop_na(poll) %>%                                                          
  select(-c("poll.s", "poll.f")) %>%                                            
  left_join(geo, by = "entity_ID") %>%                                          
  group_by(entity_ID) %>% 
  distinct(species, .keep_all = TRUE) %>%             
  ungroup()    

##############################
#### NUMBER SP V FAMILY ######
##############################

species %>% group_by(poll_level) %>% tally()

##############################
#### SP RICHNESS POLL TYPE ###
##############################

sprich.biotic <- species %>%
  filter(poll == "biotic") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(poll.b = n()) 

sprich.abiotic <- species %>%
  filter(poll == "abiotic") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(poll.ab = n()) 

################################################
##### TOTAL SP RICHNESS AND JOIN WITH POLL #####
################################################

# start form species.dat to get even removed data
sprich <- species.dat[[2]] %>% as.data.frame() %>% 
  mutate(species = paste(genus, species_epithet, sep = " ")) %>%                             
  select(c("entity_ID", "family", "native", "naturalized", "species","name_ID", "genus")) %>%        
  mutate_at(c('family', 'species', 'name_ID', 'genus'), as.character) %>%
  mutate_at(c('native', 'naturalized'), as.numeric) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  filter(native == 1) %>%                                                                     
  select(-c("native", "naturalized")) %>%                                                  
  left_join(geo, by = "entity_ID") %>%                                          
  group_by(entity_ID) %>% 
  distinct(species, .keep_all = TRUE) %>%             
  summarise(sprich = n()) %>%
  ungroup()

# join POLL sp richness, replace NAs and join with sprich, and geo
dat <- sprich %>% 
  left_join(sprich.biotic, by = "entity_ID") %>%
  left_join(sprich.abiotic, by = "entity_ID") %>%
  mutate_at(vars(c('poll.b','poll.ab')), ~replace(., is.na(.), 0)) %>%
  left_join(geo, by = "entity_ID") 

########################
####### SAVE DATA ######
########################

# RESTRICTED
saveRDS(dat, "data/native_poll_latitude_data_2023.RDS")

# EXTENDED
#saveRDS(dat, "data/Ext_native_poll_latitude_data_2023.RDS")