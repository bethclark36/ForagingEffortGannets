#PLot each colony foraging range (mean max dist)
library(ggrepel)
library(sf)     #for handling spatial points and polygons
library(rnaturalearth)# rnaturalearth package for geographic basemaps
library(tidyverse) #for plots and data wrangling
library(stringr) #for text edits

rm(list=ls())

dat <- read.csv("data_gannet_foraging_by_colony.csv")
head(dat)

#Set coordinate reference system (CRS) for colony locations (example for lat/lon, WGS84)
col_locs <- st_as_sf(dat,coords = c("Lon","Lat"), crs = 4326) 

max_dist_buffers <- st_buffer(col_locs, dist = dat$Maxdist*1000)

plot(max_dist_buffers)

ggplot(max_dist_buffers)+
  geom_sf()+
  geom_sf(data = col_locs)+
  theme_bw()

## Land polygon
#Include the name of the country you wish to import a map for
#To include multiple countries use for example "country = c("Antarctica","Chile")"
land <- ne_countries(scale = "large", returnclass = "sf",
                     continent = c("Europe", "North America")) %>% st_make_valid()

ggplot(land)+
  geom_sf()

#Remove part due to lines of crossing date line
land_cropped <- st_crop(land, xmin = -70, xmax = 30,
                          ymin = 26, ymax = 72)

map1 <- ggplot(max_dist_buffers)+
  geom_sf()+
  geom_sf(colour = "blue", fill = "blue", alpha = 0.1)+
  geom_sf(data = land_cropped, colour = "grey50", fill = "grey75")+
  #geom_sf_label(data = max_dist_buffers,aes(label = label),
   #             nudge_x = -1,nudge_y = -1) +
  geom_sf(data = col_locs, size = 0.5)+
  coord_sf(xlim = c(-64, 25), ylim = c(46, 71))+
  theme_bw()


png("col_range_size.png",
    width=2900,height=1700,res=300)
map1
dev.off()


