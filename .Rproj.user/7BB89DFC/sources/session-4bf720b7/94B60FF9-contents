#Sep 2023 Beth Clark & David Pascall

#Load packages, data and set up ####
install.packages("tidyverse") 
install.packages("cowplot") 
install.packages("DHARMa")
install.packages("spaMM")
library(tidyverse) 
library(cowplot) 
library(DHARMa)
library(spaMM)

sessionInfo() #The output is below showing the versions of R and rpackages used when this script was originally run:

#R version 4.3.3 (2024-02-29 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 11 x64 (build 22000)

#Matrix products: default


#locale:
#[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
#[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
#[5] LC_TIME=English_United Kingdom.utf8    

#time zone: Europe/London
#tzcode source: internal

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] spaMM_4.4.0     DHARMa_0.4.6    cowplot_1.1.3   lubridate_1.9.3 forcats_1.0.0   stringr_1.5.0  
#[7] dplyr_1.1.2     purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4  
#[13] tidyverse_2.0.0

#loaded via a namespace (and not attached):
#[1] utf8_1.2.3          generics_0.1.3      slam_0.1-50         stringi_1.7.12      lattice_0.22-5     
#[6] lme4_1.1-34         hms_1.1.3           magrittr_2.0.3      grid_4.3.3          timechange_0.2.0   
#[11] Matrix_1.6-1        backports_1.4.1     fansi_1.0.4         scales_1.2.1        pbapply_1.7-2      
#[16] numDeriv_2016.8-1.1 registry_0.5-1      cli_3.6.1           crayon_1.5.2        rlang_1.1.1        
#[21] ROI_1.0-1           munsell_0.5.0       splines_4.3.3       withr_2.5.0         parallel_4.3.3     
#[26] tools_4.3.3         tzdb_0.4.0          checkmate_2.2.0     nloptr_2.0.3        minqa_1.2.5        
#[31] colorspace_2.1-0    boot_1.3-29         vctrs_0.6.3         R6_2.5.1            proxy_0.4-27       
#[36] lifecycle_1.0.3     MASS_7.3-60.0.1     pkgconfig_2.0.3     pillar_1.9.0        gtable_0.3.3       
#[41] glue_1.6.2          Rcpp_1.0.11         tidyselect_1.2.0    rstudioapi_0.15.0   nlme_3.1-164       
#[46] compiler_4.3.3 

#Clear R environment
rm(list=ls())

#read in data 
dat <- read.csv("data_gannet_foraging_by_year.csv")
head(dat)

#Calculate s=the square root of colony size
dat$sqrt.size <- sqrt(dat$Size)
head(dat)

#check response variable distributions (annual means per colony)
hist(dat$Dur)     #Foraging trip duration
hist(dat$Maxdist) #Maximum distance reached on foraging trip

#Check difference between colony count & gps year####
dat$GPS_year
dat$Year

range(dat$GPS_year)

dat$GPS_year_num <- dat$GPS_year - min(dat$GPS_year) + 1

dat$count_year <- as.numeric(ifelse(dat$Year == "2013/14",2013.5,as.character(dat$Year)))

year_diff <- abs(dat$count_year - dat$GPS_year)

range(year_diff)
mean(year_diff)

plot(dat$Lon, dat$Lat)
cor.test(dat$Lon, dat$Lat)
#correlated, so consider not including Lat & Lon in same models

#Duration models ####
#Using the spaMM package to fit mixed models with spatial autocorrelation
citation("spaMM")

#Colony is a random intercept
#Lat and Lon with a matern covariance function to account for spatial autocorrelation
#"Earth" distance is used because locations are latlon and so on a sphere

size_lat_yr <- fitme(Dur ~ sqrt.size + Lat + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                     data=dat, control.dist = list(dist.method = "Earth"))

size_lon_yr <- fitme(Dur ~ sqrt.size + Lon + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                     data=dat, control.dist = list(dist.method = "Earth"))

size_lat <- fitme(Dur ~ sqrt.size + Lat + (1|colony) + Matern(1|Lon + Lat),
                  data=dat, control.dist = list(dist.method = "Earth"))

size_lon <- fitme(Dur ~ sqrt.size + Lon + (1|colony) + Matern(1|Lon + Lat), 
                  data=dat, control.dist = list(dist.method = "Earth"))

size_yr <- fitme(Dur ~ sqrt.size + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                 data=dat, control.dist = list(dist.method = "Earth"))

size <- fitme(Dur ~ sqrt.size + (1|colony) + Matern(1|Lon + Lat), 
              data=dat, control.dist = list(dist.method = "Earth"))

lat_yr <- fitme(Dur ~ Lat + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                data=dat, control.dist = list(dist.method = "Earth"))

lon_yr <- fitme(Dur ~ Lon + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                data=dat, control.dist = list(dist.method = "Earth"))

lat <- fitme(Dur ~ Lat + (1|colony) + Matern(1|Lon + Lat),  
             data=dat, control.dist = list(dist.method = "Earth"))

lon <- fitme(Dur ~ Lon + (1|colony) + Matern(1|Lon + Lat), 
             data=dat, control.dist = list(dist.method = "Earth"))

yr <- fitme(Dur ~ GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
            data=dat, control.dist = list(dist.method = "Earth"))

null <- fitme(Dur ~ (1|colony) + Matern(1|Lon + Lat),
              data=dat, control.dist = list(dist.method = "Earth"))


#store results
fixed_effects <- c("size_lat_yr","size_lon_yr","size_lat","size_lon","size_yr",
                 "size","lat_yr","lon_yr","lat","lon","yr","none")
model_results <- as.data.frame(fixed_effects)
model_results$response <- "trip_duration"

#calculate AIC - we are using marginal
AIC(size_lat_yr)

model_results$mAIC <- c(AIC(size_lat_yr)[[1]],AIC(size_lon_yr)[[1]],AIC(size_lat)[[1]],
                        AIC(size_lon)[[1]],AIC(size_yr)[[1]],AIC(size)[[1]],
                        AIC(lat_yr)[[1]],AIC(lon_yr)[[1]],AIC(lat)[[1]],AIC(lon)[[1]],
                        AIC(yr)[[1]],AIC(null)[[1]])

#calculate pseudo R2 
model_results$pseudoR2 <- c(pseudoR2(size_lat_yr)[[1]],pseudoR2(size_lon_yr)[[1]],
                            pseudoR2(size_lat)[[1]],pseudoR2(size_lon)[[1]],
                            pseudoR2(size_yr)[[1]],pseudoR2(size)[[1]],
                            pseudoR2(lat_yr)[[1]],pseudoR2(lon_yr)[[1]],
                            pseudoR2(lat)[[1]],pseudoR2(lon)[[1]],
                            pseudoR2(yr)[[1]],pseudoR2(null)[[1]])
warnings()
#Default null model formula may not be appropriate for mixed-effect models.
#Note this in the results table

#Round results
model_results$mAIC <- round(model_results$mAIC, digits = 2)
model_results$pseudoR2 <- round(model_results$pseudoR2, digits = 3)

#Sort by AIC
model_results <- model_results %>% arrange(mAIC)
model_results$delta_mAIC <- model_results$mAIC - min(model_results$mAIC)

model_results

#size_lat is the top model, so check the fit using DHARMa
plot(DHARMa::simulateResiduals(size_lat))

#no significant problems, so use to predict for plots

#Duration plots predicting latitude ####

dat$Dur_original <- dat$Dur
dat$Dur <- dat$Dur_original*10 #rescale for plotting
size_lat <- fitme(Dur ~ sqrt.size + Lat + (1|colony) + Matern(1|Lon + Lat),
                  data=dat, control.dist = list(dist.method = "Earth"))

border <- 1
newdat <- tidyr::expand(data = dat,
                        Lat = seq((min(Lat)-border), 
                                  (max(Lat)+border),length=100),
                        sqrt.size = mean(sqrt.size)) %>%
  data.frame()

pred <- predict(size_lat, newdat, 
                re.form=NA, intervals = "fixefVar")
pred <- as.data.frame(cbind(pred, attr(pred, "intervals")))
newdat$pred <- pred$V1 
newdat$lower <- pred$fixefVar_0.025
newdat$upper <- pred$fixefVar_0.975

unique(dat$colony)

#correct character coding error
dat$colony <- ifelse(dat$colony == "Skr\xfa\xf0ur" ,"Skr??ur",dat$colony)
dat$colony <- ifelse(dat$colony == "Store Ulv\xf8yhomen" ,"Store Ulv?yhomen",dat$colony)

#set colours
dat$col_col <- ifelse(as.character(dat$colony) %in% 
                        c("Ailsa Craig","Baccalieu Island","Bonaventure Island",
                          "Bull Rock","Heligoland","Lambay","Little Skellig",
                          "St Kilda","Sule Skerry"),"grey","black")

#set shapes
dat$col_shape <- ifelse(as.character(dat$colony) %in% 
                          c("Ailsa Craig","Baccalieu Island","Bonaventure Island",
                            "Bull Rock","Heligoland","Lambay","Little Skellig",
                            "St Kilda","Sule Skerry"),16,0)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Alderney",18, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Bass Rock",2, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Bempton",3, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Funk Island",4, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Grassholm",1, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Great Saltee",17, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Rouzic",6, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Skr??ur",0, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Store Ulv?yhomen",25, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Storstappen",8, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Vestmannaeyjar",15, dat$col_shape)
dat$col_shape <- ifelse(as.character(dat$colony) == 
                          "Cape St Mary's",5, dat$col_shape)

#add theme elements for all plots
ms_theme <- theme_bw()+
  theme(text = element_text(size=20)) +
  theme(axis.text=element_text(colour="black")) +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())



head(newdat)
dl3 <-ggplot(newdat, aes(x=Lat, y=pred)) + 
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              colour = "grey", alpha = 0.2) +  
  geom_point(data = dat, aes(x=Lat, y=Dur),
             shape = dat$col_shape,
             fill = dat$col_col,
             color = dat$col_col,
             size=3.5,stroke=1.1)+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 425)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(min(dat$Lat)- border, 
                              max(dat$Lat)+ border)) +
  ms_theme;dl3
#note duration is hours * 10, to match the number of digits to max dist
#axis labelled are added later

#Duration plots predict colony size ####
border <- 11.74734
newdat <- tidyr::expand(data = dat,
                        sqrt.size = seq(0, 
                                        (max(sqrt.size) + border),length=100),
                        Lat = mean(Lat))
pred <- predict(size_lat, newdat, 
                re.form=NA, intervals = "fixefVar")
pred <- as.data.frame(cbind(pred, attr(pred, "intervals")))
newdat$pred <- pred$V1 
newdat$lower <- pred$fixefVar_0.025
newdat$upper <- pred$fixefVar_0.975

ds3 <-ggplot(newdat, aes(x=sqrt.size, y=pred)) + 
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              colour = "grey", alpha = 0.2) +  
  geom_point(data = dat, aes(x=sqrt.size, y=Dur),
             shape = dat$col_shape,
             fill = dat$col_col,
             color = dat$col_col,
             size=3.5,stroke=1.1)+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 425)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(0, 
                              max(dat$sqrt.size)+ border)) +
  ms_theme;ds3

plot_grid(ds3,dl3,ncol = 2)

#Maximum distance models####
size_lat_yr <- fitme(Maxdist ~ sqrt.size + Lat + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                     data=dat, control.dist = list(dist.method = "Earth"))

size_lon_yr <- fitme(Maxdist ~ sqrt.size + Lon + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                     data=dat, control.dist = list(dist.method = "Earth"))

size_lat <- fitme(Maxdist ~ sqrt.size + Lat + (1|colony) + Matern(1|Lon + Lat),
                  data=dat, control.dist = list(dist.method = "Earth"))

size_lon <- fitme(Maxdist ~ sqrt.size + Lon + (1|colony) + Matern(1|Lon + Lat), 
                  data=dat, control.dist = list(dist.method = "Earth"))

size_yr <- fitme(Maxdist ~ sqrt.size + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                 data=dat, control.dist = list(dist.method = "Earth"))

size <- fitme(Maxdist ~ sqrt.size + (1|colony) + Matern(1|Lon + Lat), 
              data=dat, control.dist = list(dist.method = "Earth"))

lat_yr <- fitme(Maxdist ~ Lat + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                data=dat, control.dist = list(dist.method = "Earth"))

lon_yr <- fitme(Maxdist ~ Lon + GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
                data=dat, control.dist = list(dist.method = "Earth"))

lat <- fitme(Maxdist ~ Lat + (1|colony) + Matern(1|Lon + Lat),  
             data=dat, control.dist = list(dist.method = "Earth"))

lon <- fitme(Maxdist ~ Lon + (1|colony) + Matern(1|Lon + Lat), 
             data=dat, control.dist = list(dist.method = "Earth"))

yr <- fitme(Maxdist ~ GPS_year_num + (1|colony) + Matern(1|Lon + Lat), 
            data=dat, control.dist = list(dist.method = "Earth"))

null <- fitme(Maxdist ~ (1|colony) + Matern(1|Lon + Lat),
              data=dat, control.dist = list(dist.method = "Earth"))

#store results
model_results_dur <- model_results

model_results <- as.data.frame(fixed_effects)
model_results$response <- "maximum_distance"

#calculate AIC - we are using marginal
model_results$mAIC <- c(AIC(size_lat_yr)[[1]],AIC(size_lon_yr)[[1]],AIC(size_lat)[[1]],
                        AIC(size_lon)[[1]],AIC(size_yr)[[1]],AIC(size)[[1]],
                        AIC(lat_yr)[[1]],AIC(lon_yr)[[1]],AIC(lat)[[1]],AIC(lon)[[1]],
                        AIC(yr)[[1]],AIC(null)[[1]])

#calculate pseudo R2 
model_results$pseudoR2 <- c(pseudoR2(size_lat_yr)[[1]],pseudoR2(size_lon_yr)[[1]],
                            pseudoR2(size_lat)[[1]],pseudoR2(size_lon)[[1]],
                            pseudoR2(size_yr)[[1]],pseudoR2(size)[[1]],
                            pseudoR2(lat_yr)[[1]],pseudoR2(lon_yr)[[1]],
                            pseudoR2(lat)[[1]],pseudoR2(lon)[[1]],
                            pseudoR2(yr)[[1]],pseudoR2(null)[[1]])
warnings()
#Default null model formula may not be appropriate for mixed-effect models.
#Note this in the results table

#Round results
model_results$mAIC <- round(model_results$mAIC, digits = 2)
model_results$pseudoR2 <- round(model_results$pseudoR2, digits = 3)

#Sort by AIC
model_results <- model_results %>% arrange(mAIC)
model_results$delta_mAIC <- model_results$mAIC - min(model_results$mAIC)

model_results

#size_lat is not the top model, but is within 0.18 AIC and has a bigger pseudo r2,
#so plot to match trip duration
#check the fit using DHARMa
plot(DHARMa::simulateResiduals(size_lat))

#good, so plot
#Predict max dist v latitude ####
border <- 1
newdat <- tidyr::expand(data = dat,
                        Lat = seq((min(Lat)-border), 
                                  (max(Lat)+border),length=100),
                        sqrt.size = mean(sqrt.size)) %>%
  data.frame()
pred <- predict(size_lat, newdat, 
                re.form=NA, intervals = "fixefVar")
pred <- as.data.frame(cbind(pred, attr(pred, "intervals")))
newdat$pred <- pred$V1 
newdat$lower <- pred$fixefVar_0.025
newdat$upper <- pred$fixefVar_0.975


head(newdat)
ml3 <-ggplot(newdat, aes(x=Lat, y=pred)) + 
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              colour = "grey", alpha = 0.2) +  
  geom_point(data = dat, aes(x=Lat, y=Maxdist),
             shape = dat$col_shape,
             fill = dat$col_col,
             color = dat$col_col,
             size=3.5,stroke=1.1)+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 330)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(min(dat$Lat)- border, 
                              max(dat$Lat)+ border)) +
  ms_theme;ml3

#Predict size max ####
border <- 11.74734
newdat <- tidyr::expand(data = dat,
                        sqrt.size = seq(0, 
                                        (max(sqrt.size) + border),length=100),
                        Lat = mean(Lat))
pred <- predict(size_lat, newdat, 
                re.form=NA, intervals = "fixefVar")
pred <- as.data.frame(cbind(pred, attr(pred, "intervals")))
newdat$pred <- pred$V1 
newdat$lower <- pred$fixefVar_0.025
newdat$upper <- pred$fixefVar_0.975

ms3 <-ggplot(newdat, aes(x=sqrt.size, y=pred)) + 
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              colour = "grey", alpha = 0.2) +  
  geom_point(data = dat, aes(x=sqrt.size, y=Maxdist),
             shape = dat$col_shape,
             fill = dat$col_col,
             color = dat$col_col,
             size=3.5,stroke=1.1)+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 330)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(0, 
                              max(dat$sqrt.size)+ border)) +
  ms_theme;ms3

#plot grids ####
plot_grid(ds3, dl3, ms3, ml3, ncol=2)

#save plots
png("plots_by_year.png",
    width=5000,height=5000,res=400)
plot_grid(ds3, dl3, ms3, ml3, ncol=2)
dev.off()

#save model results
model_results_all <- rbind(model_results_dur,model_results)
write.csv(model_results_all,"result_by_year.csv",
          row.names = F)
