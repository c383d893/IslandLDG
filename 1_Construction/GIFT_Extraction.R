library("GIFT")
library("tidyverse")

api <- "https://gift.uni-goettingen.de/api/extended/"
# Using the "extended" api of GIFT offers access to all publicly available data 
# in GIFT. The analyses of this paper have been based on these data plus a few 
# restricted resources. Access to these requires contacting the authors of this 
# paper to get permission from the primary data contributors and access to the 
# password protected restricted api:

# api <- "https://gift.uni-goettingen.de/api/restricted/"

# References
ref <- GIFT_references(GIFT_version = "3.0", api = api)

head(ref)
ref <- ref[which(ref$restricted == 1), ] # We exclude most restricted resources
#  But we keep: WCVP and World Ferns
ref <- ref[which(!ref$ref_ID %in% c(10647, 10601)), ]

# Download all Checklists of covered for native Tracheophyta
checklists <- GIFT_checklists(
  taxon_name = "Tracheophyta",
  complete_taxon = TRUE,
  floristic_group = "native",
  complete_floristic = TRUE,
  geo_type = "All",
  ref_excluded = unique(ref$ref_ID),
  suit_geo = TRUE,
  shp = NULL, coordinates = NULL, overlap = "centroid_inside",
  remove_overlap = TRUE,
  area_threshold_island = 0, 
  area_threshold_mainland = 100, 
  overlap_th = 0.1,
  namesmatched = TRUE,
  list_set_only = FALSE,
  GIFT_version = "3.0",
  api = api)

# Environmental data
env <- GIFT_env(
  entity_ID = unique(checklists[[1]]$entity_ID),
  miscellaneous = c("area", "longitude", "latitude", "botanical_continent",
                    "no_polygons", "biome", "dist"),
  rasterlayer = c("wc2.0_bio_30s_01", "wc2.0_bio_30s_12", "GDP_satellite",
                  "gpw-v4_popdensity_UN_2015", "HFP2009_unproj",
                  "hii_v2geo"),
  sumstat = "mean",
  GIFT_version = "3.0",
  api = api)

env2 <- GIFT_env(
  entity_ID = unique(checklists[[1]]$entity_ID),
  rasterlayer = c("mx30_grd", "mi30_grd"),
  sumstat = list("max", "min"),
  GIFT_version = "3.0",
  api = api)
env2$range <- env2$max_mx30_grd - env2$min_mi30_grd

env <- full_join(env, env2[, c("entity_ID", "range")], by = "entity_ID")

# native and naturalized species richness
tra_meta <- GIFT_traits_meta(GIFT_version = "3.0",
                             api = api)
rich <- GIFT_richness(taxon_name = "Tracheophyta",
                      GIFT_version = "3.0", api = api)

env <- left_join(env, rich[, c("entity_ID", "native", "naturalized")],
                 by = "entity_ID")

# Geology (restricted)
geol <- GIFT_regions(GIFT_version = "3.0", api = api)
geol <- geol[,-c(2,4,5,6,7,9,11,12,13)]

env <- left_join(env, geol,
                 by = "entity_ID")



# Pollination and dispersal syndrome
# species level
gift_traits <- GIFT_traits(trait_IDs = c("3.6.1"), GIFT_version = "3.0", api = api)

library(GIFT)
gift_traits_tax <- GIFT_traits_tax(trait_IDs = c("3.6.1"),
                 GIFT_version = "3.0", api = api)



# Get references for traits
ref_all <- GIFT_references(GIFT_version = "3.0", api = api)

refs_traits_3.6.1 <- ref_all[which(ref_all$ref_ID %in% 
                                     unique(unlist(strsplit(paste(gift_traits$references_3.6.1[!is.na(gift_traits$references_3.6.1)], collapse = ","),",")))
),]


save(checklists, env, gift_traits, gift_traits_tax, refs_traits_3.6.1, 
    file = "GIFT_LatGrad_extended.RData")