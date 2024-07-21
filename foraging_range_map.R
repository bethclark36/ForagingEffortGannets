#PLot each colony foraging range (mean max dist)
#Beth Clark
library(sf)     #for handling spatial points and polygons
library(rnaturalearth)# rnaturalearth package for geographic basemaps
library(tidyverse) #for plots and data wrangling

sessionInfo() 
#R version 4.3.3 (2024-02-29 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 11 x64 (build 22000)

#Matrix products: default

#locale:
#[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8
#[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.utf8    

#time zone: Europe/London
#tzcode source: internal

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] sf_1.0-14           lubridate_1.9.3     forcats_1.0.0       stringr_1.5.0       dplyr_1.1.2         purrr_1.0.1        
#[7] readr_2.1.4         tidyr_1.3.0         tibble_3.2.1        tidyverse_2.0.0     rnaturalearth_0.3.4 ggrepel_0.9.3      
#[13] ggplot2_3.4.4      

#loaded via a namespace (and not attached):
#[1] gtable_0.3.3             lattice_0.22-5           tzdb_0.4.0               vctrs_0.6.3              tools_4.3.3             
#[6] generics_0.1.3           proxy_0.4-27             fansi_1.0.4              pkgconfig_2.0.3          Matrix_1.6-1            
#[11] KernSmooth_2.23-22       DHARMa_0.4.6             lifecycle_1.0.3          compiler_4.3.3           farver_2.1.1            
#[16] munsell_0.5.0            mitools_2.4              survey_4.2-1             class_7.3-22             pillar_1.9.0            
#[21] nloptr_2.0.3             crayon_1.5.2             MASS_7.3-60.0.1          classInt_0.4-9           wk_0.7.3                
#[26] boot_1.3-29              nlme_3.1-164             tidyselect_1.2.0         digest_0.6.33            stringi_1.7.12          
#[31] pander_0.6.5             splines_4.3.3            grid_4.3.3               colorspace_2.1-0         cli_3.6.1               
#[36] magrittr_2.0.3           survival_3.5-8           utf8_1.2.3               e1071_1.7-13             withr_2.5.0             
#[41] scales_1.2.1             jtools_2.2.2             sp_2.0-0                 timechange_0.2.0         httr_1.4.6              
#[46] lme4_1.1-34              rnaturalearthhires_0.2.1 hms_1.1.3                s2_1.1.4                 rlang_1.1.1             
#[51] Rcpp_1.0.11              glue_1.6.2               DBI_1.1.3                rstudioapi_0.15.0        minqa_1.2.5             
#[56] jsonlite_1.8.7           R6_2.5.1                 units_0.8-2             

rm(list=ls())

#Read in the data
dat <- read.csv("data_gannet_foraging_by_colony.csv")
head(dat)

#Set coordinate reference system (CRS) for colony locations (example for lat/lon, WGS84)
col_locs <- st_as_sf(dat,coords = c("Lon","Lat"), crs = 4326) 

#create a buffer, dist is in metres. 
max_dist_buffers <- st_buffer(col_locs, dist = dat$Maxdist*1000)

## Land polygon
#Include the name of the country you wish to import a map for
#To include multiple countries use for example "country = c("Antarctica","Chile")"
land <- ne_countries(scale = "large", returnclass = "sf",
                     continent = c("Europe", "North America")) %>% st_make_valid()

#Remove part due to lines of crossing date line
land_cropped <- st_crop(land, xmin = -70, xmax = 30,
                          ymin = 26, ymax = 72)

map1 <- ggplot(max_dist_buffers)+
  geom_sf()+
  geom_sf(colour = "blue", fill = "blue", alpha = 0.1)+
  geom_sf(data = land_cropped, colour = "grey75", fill = "grey75")+
  geom_sf(data = col_locs, size = 0.5)+
  coord_sf(xlim = c(-64, 25), ylim = c(46.5, 71))+
  theme_bw(); map1

png("col_range_size.png",
    width=2900,height=1600,res=300)
map1
dev.off()
