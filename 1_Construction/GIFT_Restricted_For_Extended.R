# Extract geology and age from full dataset for extended dataset
library("tidyverse")

# bring in GIFT data:
load("data/GIFT_LatGrad_restricted.RData")
species.dat <- checklists
geo.dat <- env

# bring in GIFT data:
load("data/GIFT_LatGrad_extended.RData") 
species.dat.ext <- checklists 
geo.dat.ext <- env

colnames(species.dat[[1]])
colnames(species.dat.ext[[1]])
colnames(geo.dat)
colnames(geo.dat.ext)

# extract age_MA and geology from restricted for locations in extended
geo.dat.forext <- geo.dat %>% filter(entity_ID %in% geo.dat.ext$entity_ID)

# save
saveRDS(geo.dat.forext,"data/geo.dat.extended_foranalyses.rds")

length(unique(geo.dat$entity_ID)) # 1780
length(unique(geo.dat.ext$entity_ID)) #1726
length(unique(geo.dat.forext$entity_ID)) #1709