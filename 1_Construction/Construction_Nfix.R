# This script outputs counts of each nfixing type (yes or no)
# First based on species assignments
# Next based on family level proportion applied to each flora

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)
packageVersion("tidyverse") # ‘1.3.2’

###################
#### LOAD DATA ####
###################

# bring in nfix dat species, with family
ndat <- read.csv("data/Werner_NFix.csv", header = TRUE) %>%     
  dplyr::select(c("species", "family", "data_fixing")) %>% 
  mutate(species = str_replace(species, "_", " "))      

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

# with species checklist data, keep native
# join with geo data and keep unique sp by loc (entity_ID)
species <- species.dat[[2]] %>% as.data.frame() %>% 
  mutate(species = paste(genus, species_epithet, sep = " ")) %>% 
  select(c("entity_ID", "family", "native", "naturalized", "species", "name_ID", "genus")) %>%
  mutate_at(c('family', 'species', 'name_ID', 'genus'), as.character) %>%
  mutate_at(c('native', 'naturalized'), as.numeric) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  filter(native == 1) %>%                                                                     
  select(-c("native", "naturalized")) %>%                                                  
  left_join(ndat, by = c("species","family")) %>%  
  left_join(geo, by = "entity_ID") %>%                                          
  group_by(entity_ID) %>% 
  distinct(species, .keep_all = TRUE) %>%             
  select("species","family","data_fixing") %>%
  ungroup()   

# subset to known nfix species
# keep only known status
species.known.nfix <- species %>%
  drop_na(data_fixing) 

# subset to the other species
species.unknown.nfix <- species %>%
  anti_join(species.known.nfix)

###############################
#### SP RICHNESS NFIX TYPE ####
###############################

sprich.nfix <- species.known.nfix %>%
  filter(data_fixing == "Yes") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(nfix.y = n())  

sprich.nonfix <- species.known.nfix %>%
  filter(data_fixing == "No") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(nfix.n = n())  

#Join together
sprich.sp <- sprich.nfix %>% 
  full_join(sprich.nonfix, by = "entity_ID") %>%                   
  mutate(nfix.y  = ifelse(is.na(nfix.y), 0, nfix.y)) %>%  
  mutate(nfix.n = ifelse(is.na(nfix.n), 0, nfix.n))  

######################################
######## FAMILY PROPORTIONS ##########
######## (NUMBER SP ASSIGNED) ########
######################################

# take species and collapse to single rep of each species
# keep only known status

unique.species <- species %>%
  distinct(species, .keep_all = TRUE) %>% 
  drop_na(data_fixing) %>%
  select("species","family","data_fixing") %>%
  ungroup()    

# get family prop nfix
fam.nfix <- unique.species %>%
  group_by(family, data_fixing) %>%                 
  summarise(n = n()) %>%                            
  mutate(p = n / sum(n)) %>%                        
  filter(data_fixing == "Yes") %>%                  
  dplyr::select(-data_fixing) %>%                   
  rename(fam.nfix.sum = n, fam.nfix.prop = p) %>%       
  ungroup() 

# get family prop nonfix
fam.nonfix <- unique.species %>%
  group_by(family, data_fixing) %>%                 
  summarise(n = n()) %>%                            
  mutate(p = n / sum(n)) %>%                        
  filter(data_fixing == "No") %>%                  
  dplyr::select(-data_fixing) %>%                   
  rename(fam.nonfix.sum = n, fam.nonfix.prop = p) %>%       
  ungroup() 

# join together
merge.fam.fix <- fam.nfix %>%
  full_join(fam.nonfix, by = "family") %>%                                        
  dplyr::select(-c(fam.nfix.sum, fam.nonfix.sum)) %>%                                 
  mutate(fam.nfix.prop = ifelse(is.na(fam.nfix.prop), 0, fam.nfix.prop)) %>%            
  mutate(fam.nonfix.prop = ifelse(is.na(fam.nonfix.prop), 0, fam.nonfix.prop)) %>%      
  mutate(fam.nfix.prop = ifelse(fam.nfix.prop == 1, 0, fam.nfix.prop)) %>%            
  mutate(fam.nonfix.prop = ifelse(fam.nonfix.prop == 1, 0, fam.nonfix.prop))             

###############################
#### SP RICHNESS NFIX TYPE ####
####### FOR SP UNKOWN #########
#### BASED ON FAMILY PROPS ####
###############################

# find counts based on prop for unknown nfix status at the species level
sprich.nfix.fam <- species.unknown.nfix %>% 
  left_join(merge.fam.fix, by = "family") %>%             
  filter(is.na(data_fixing)) %>%                         
  group_by(entity_ID) %>%                                
  summarise(sum = sum(fam.nfix.prop, na.rm = T)) %>%       
  mutate(sum = ceiling(sum)) %>%
  rename(nfix.y.f = sum) %>%
  ungroup()

sprich.nonfix.fam <- species.unknown.nfix %>% 
  left_join(merge.fam.fix, by = "family") %>%             
  filter(is.na(data_fixing)) %>%                          
  group_by(entity_ID) %>%                                
  summarise(sum = sum(fam.nonfix.prop, na.rm = T)) %>%      
  mutate(sum = ceiling(sum)) %>%
  rename(nfix.n.f = sum) %>%
  ungroup()

#Merge numbers of N-fixing and non-Nfixing families per entity
sprich.fam <- sprich.nfix.fam %>%
  full_join(sprich.nonfix.fam, by = "entity_ID")

################################
#### NUMBER ASSIGNED FAMILY ####
################################

sprich.nfix.fam.species <- species.unknown.nfix %>% 
  left_join(merge.fam.fix, by = "family") %>%             
  filter(is.na(data_fixing)) %>%
  distinct(species, .keep_all = TRUE) %>%
  filter(!is.na(fam.nfix.prop))

################################################
##### TOTAL SP RICHNESS AND JOIN WITH NFIX #####
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

dat <- sprich %>%
  full_join(sprich.sp, by = "entity_ID") %>%
  full_join(sprich.fam, by = "entity_ID") %>%
  filter(!is.na(entity_ID)) %>%                                       
  replace(is.na(.), 0) %>%                                           
  mutate(nfix.y.final = nfix.y + nfix.y.f) %>%                       
  mutate(nfix.n.final= nfix.n + nfix.n.f) %>%                  
  dplyr::select(c(entity_ID, sprich, nfix.y.final, nfix.n.final)) %>%       
  rename(nfix = nfix.y.final, nonfix = nfix.n.final) %>%
  left_join(geo, by = "entity_ID")          

########################
####### SAVE DATA ######
########################

# RESTRICTED
saveRDS(dat, "data/native_nfix_latitude_data_2023.RDS")

# EXTENDED
#saveRDS(dat, "data/Ext_native_nfix_latitude_data_2023.RDS")
