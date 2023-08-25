# This script outputs counts of each myc type
# First based on species assignments
# Next based on genus level

##########################
###### LOAD PACKAGES #####
##########################

library(tidyverse)
packageVersion("tidyverse") # ‘1.3.2’

###################
#### LOAD DATA ####
###################

# bring in myc dat species (vascular plants only)
mdat.sp <- read.csv("data/Funroot_species_08092021.csv", header = TRUE) %>%                     
  rename(species = latin, myc.s = myc) %>%                                                     
  mutate(species = str_replace(species, "_", " "))                                         

# bring in myc dat genus (vascular plants only)
mdat.g <- read.csv("data/Funroot_genus_08092021.csv", header = TRUE) %>%                       
  rename(myc.g = myc)  

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

# with species checklist data, keep native, assign myc preferentially using species, then genus, condense myc types
# join with geo data and keep unique sp by loc (entity_ID)
species <- species.dat[[2]] %>% as.data.frame() %>% 
  mutate(species = paste(genus, species_epithet, sep = " ")) %>% 
  select(c("entity_ID", "family", "native", "naturalized", "species", "name_ID", "genus")) %>%
  mutate_at(c('family', 'species', 'name_ID', 'genus'), as.character) %>%
  mutate_at(c('native', 'naturalized'), as.numeric) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  filter(native == 1) %>%                                                                     
  select(-c("native", "naturalized")) %>%                                                  
  left_join(mdat.sp, by = c("species")) %>%                                    
  left_join(mdat.g, by = c("genus"))  %>%                                     
  mutate(myc = ifelse(is.na(myc.s), myc.g, myc.s)) %>%  
  mutate(myc_level = ifelse(!is.na(myc.s), "species", "genus")) %>%
  drop_na(myc) %>%                                                          
  mutate(myc = case_when(myc=="AM" ~ "AM",                                   
                         myc=="AMNM" ~ "AM", # if plant is AMNM: some dependency on AM
                         myc=="AMEM" ~ "AM", # same above.
                         myc=="EM" ~ "EM",
                         myc=="NM" ~ "NM",
                         myc=="ORC" ~ "ORC")) %>%
  select(-c("myc.s", "myc.g")) %>%                                            
  left_join(geo, by = "entity_ID") %>%                                          
  group_by(entity_ID) %>% 
  distinct(species, .keep_all = TRUE) %>%             
  ungroup()    

##############################
#### NUMBER SP V GENUS #######
##############################

species %>% group_by(myc_level) %>% tally()

##############################
#### SP RICHNESS MYC TYPE ####
##############################

sprich.AM <- species %>%
  filter(myc == "AM") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(AM = n()) 

sprich.EM <- species %>%
  filter(myc == "EM") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(EM = n())   

sprich.ORC <- species %>%
  filter(myc == "ORC") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(ORC = n())

sprich.NM <- species %>%
  filter(myc == "NM") %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(NM = n()) 

###############################################
##### TOTAL SP RICHNESS AND JOIN WITH MYC #####
###############################################

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

# join myc sp richness, replace NAs and join with sprich, and geo
dat <- sprich %>% 
  left_join(sprich.AM, by = "entity_ID") %>%
  left_join(sprich.EM, by = "entity_ID") %>%
  left_join(sprich.ORC, by = "entity_ID") %>%
  left_join(sprich.NM, by = "entity_ID") %>%
  mutate_at(vars(c('AM','EM','ORC','NM')), ~replace(., is.na(.), 0)) %>%
  left_join(geo, by = "entity_ID") 

########################
####### SAVE DATA ######
########################

# RESTRICTED
saveRDS(dat, "data/native_myc_latitude_data_2023.RDS")

# EXTENDED
#saveRDS(dat, "data/Ext_native_myc_latitude_data_2023.RDS")
