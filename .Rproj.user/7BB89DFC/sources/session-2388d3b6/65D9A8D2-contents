#Sep 2023 Beth Clark & David Pascall

#Load packages, data and set up ####
#install.packages("jtools") 
#install.packages("survey") 
#install.packages("tidyverse") 
#install.packages("DHARMa")
#install.packages("ggrepel")
#install.packages("cowplot") 

library(jtools)
library(survey)
library(tidyverse)
library(DHARMa)
library(ggrepel)
library(cowplot)

sessionInfo() #The output is below showing the versions of R and rpackages 
              #used when this script was originally run:

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
#[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggrepel_0.9.3   survey_4.2-1    survival_3.5-8  Matrix_1.6-1    jtools_2.2.2    spaMM_4.4.0    
#[7] DHARMa_0.4.6    cowplot_1.1.3   lubridate_1.9.3 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2    
#[13] purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0

#loaded via a namespace (and not attached):
#[1] gtable_0.3.3        lattice_0.22-5      tzdb_0.4.0          numDeriv_2016.8-1.1 vctrs_0.6.3        
#[6] tools_4.3.3         generics_0.1.3      parallel_4.3.3      proxy_0.4-27        fansi_1.0.4        
#[11] pkgconfig_2.0.3     checkmate_2.2.0     lifecycle_1.0.3     farver_2.1.1        compiler_4.3.3     
#[16] munsell_0.5.0       mitools_2.4         gap.datasets_0.0.5  codetools_0.2-19    httpuv_1.6.11      
#[21] htmltools_0.5.5     later_1.3.1         pillar_1.9.0        nloptr_2.0.3        crayon_1.5.2       
#[26] ellipsis_0.3.2      MASS_7.3-60.0.1     iterators_1.0.14    boot_1.3-29         foreach_1.5.2      
#[31] mime_0.12           nlme_3.1-164        digest_0.6.33       tidyselect_1.2.0    stringi_1.7.12     
#[36] slam_0.1-50         pander_0.6.5        labeling_0.4.2      splines_4.3.3       fastmap_1.1.1      
#[41] colorspace_2.1-0    cli_3.6.1           magrittr_2.0.3      utf8_1.2.3          withr_2.5.0        
#[46] promises_1.2.0.1    scales_1.2.1        backports_1.4.1     registry_0.5-1      timechange_0.2.0   
#[51] lme4_1.1-34         hms_1.1.3           shiny_1.7.4.1       ROI_1.0-1           pbapply_1.7-2      
#[56] gap_1.5-1           doParallel_1.0.17   mgcv_1.9-1          rlang_1.1.1         Rcpp_1.0.11        
#[61] DBI_1.1.3           xtable_1.8-4        glue_1.6.2          qgam_1.3.4          rstudioapi_0.15.0  
#[66] minqa_1.2.5         R6_2.5.1            plyr_1.8.8   

rm(list=ls())
df <- read.csv("data_gannet_foraging_by_colony.csv")
head(df)

df$sqrt.size <- sqrt(df$Size)

plot(df$Lon,df$Lat)
cor.test(df$Lon,df$Lat)
#correlated, so consider not including Lat & Lon in same models

#Survey design set up to use the survey package, which allows for 
#the finite population correction to be applied
dat <- df[,c("Lat","Lon","sqrt.size","Dur","Maxdist","GPS_n_birds")]
design <- survey::svydesign(ids=~0, fpc=rep(60,nrow(dat)), data = dat)
#60 refers to a maximum possible number of northern gannet colonies worldwide

#duration models ####
fit1<-svyglm(Dur ~ sqrt.size + Lat, design = design)
fit2<-svyglm(Dur ~ sqrt.size, design = design)
fit3<-svyglm(Dur ~ Lat, design = design)

anova(fit1, fit2) #lat effect
anova(fit1, fit3) #size effect

par(mfrow=c(2,2));plot(fit1);par(mfrow=c(1,1))

#Check residuals for the best model (fit1)
plot(simulateResiduals(fit1))
#DHARMs results are good, but give a warning:
#"fittedModel not in class of supported models. Absolutely no guarantee that this will work!"

#check with a normal lm of the same structure 
fitlm<-lm(Dur ~ sqrt.size + Lat, data = dat)
plot(simulateResiduals(fitlm)) #same result

#record model statistics
summary(fit1) 

jtools::summ(fit1) #size+lat
summ(fit2) #size
summ(fit3) #lat

summ(fitlm)

#all
0.62  #adj r2
0.62-0.21 #delta size
0.62-0.54 #delta lat

confint(fit1)

#max dist ####
fit1<-svyglm(Maxdist ~ sqrt.size + Lat, design = design)
fit2<-svyglm(Maxdist ~ sqrt.size, design = design)
fit3<-svyglm(Maxdist ~ Lat, design = design)

anova(fit1, fit2) #lat effect
anova(fit1, fit3) #size effect

par(mfrow=c(2,2));plot(fit1);par(mfrow=c(1,1))

#Check residuals for the best model (fit1)
plot(simulateResiduals(fit1))
#DHARMs results are good, but give a warning:
#"fittedModel not in class of supported models. Absolutely no guarantee that this will work!"

summary(fit1) 
summ(fit1)
summ(fit2)
summ(fit3)

#all
0.72  #adj r2
0.72-0.70 #delta lat
0.72-0.10 #delta size

confint(fit1)

#dist v dur ####
fit1<-svyglm(Maxdist ~ Dur, design = design)
fit2<-svyglm(Maxdist ~ 1, design = design)
anova(fit1,fit2)
summary(fit1)
summ(fit1)
plot(simulateResiduals(fit1))


#Plots ####
df$label <- c("Cape St Mary's",
               "Baccalieu Is",
               "Bonaventure Is",
               "Rouzic",
               "Alderney",
               "Funk Is",
               "Bull Rock",
               "Grassholm",
               "Little Skellig",
               "Great Saltee",
               "Lambay",
               "Bempton",
               "Heligoland",
               "Ailsa Craig",
               "Bass Rock",
               "St Kilda",
               "Sule Skerry",
               "Vestmannaeyjar",
               "Skrúður",
               "St Ulvøyhomen",
               "Storstappen") 


dat <- df[,c("Lat","Lon","sqrt.size","Dur","Maxdist","GPS_n_birds","label")]
dat$Dur_original <- dat$Dur
dat$Dur <- dat$Dur_original*10 #rescale for plotting
design <- survey::svydesign(ids=~0, fpc=rep(60,nrow(dat)), data = dat)

Max_max <- 212.5
Dur_max <- 300

#Predict latitude
border <- 1
newdat <- tidyr::expand(data = dat,
                        Lat = seq((min(Lat)-border), 
                                  (max(Lat)+border),length=100),
                        sqrt.size = mean(sqrt.size))

fit_dur <- svyglm(Dur ~ sqrt.size + Lat, design = design)

summary(fit_dur)
summ(fit_dur)
pred <- predict(fit_dur, newdat, re.form=NA, se = TRUE)
pred <- as.data.frame(pred)

newdat$newDur <- pred$link
newdat$lwrDur <- pred$link - (1.96*pred$SE)
newdat$uprDur <- pred$link + (1.96*pred$SE)

#subtract the effect of colony size
dat$dur_corrected_for_colsize <- dat$Dur - 
  (fit_dur$coefficients[[2]]*dat$sqrt.size) + 
  (fit_dur$coefficients[[2]]*mean(dat$sqrt.size))

dat2 <- gather(dat,value = trip_dur,
               key = type,
               Dur,
               dur_corrected_for_colsize)

ms_theme <- theme_bw()+
  theme(text = element_text(size=20)) +
  theme(axis.text=element_text(colour="black")) +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())

d3 <- ggplot(dat2) + 
  geom_ribbon(data=newdat, aes(x=Lat, ymin = lwrDur, ymax = uprDur), 
              fill = "grey", alpha = 0.2) +  
  geom_line(data = dat2, aes(x=Lat, y = trip_dur, group=label),colour="lightgrey",linetype = 2) +
  geom_line(data=newdat,aes(x=Lat, y=newDur)) +  
  geom_point(data=dat2,aes(x=Lat, y = trip_dur, group=label, colour=type),size=3) +
  #geom_text_repel(data = dat, aes(x=Lat, y = dur_corrected_for_colsize, label = label)) +
  scale_color_manual(values=c("lightgrey","black"))+
  ms_theme+
  theme(axis.title = element_blank()) +
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0), limits=c(0, Dur_max)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c((min(dat$Lat)-border), 
                              (max(dat$Lat)+border))); d3 


fit_max <- svyglm(Maxdist ~ sqrt.size + Lat, design = design)
summary(fit_max)
pred <- predict(fit_max, newdat, re.form=NA, se = TRUE)
pred <- as.data.frame(pred)
newdat$newMax <- pred$link
newdat$lwrMax <- pred$link - (1.96*pred$SE)
newdat$uprMax <- pred$link + (1.96*pred$SE)

dat$Maxdist_corrected_for_colsize <- dat$Maxdist - 
  (fit_max$coefficients[[2]]*dat$sqrt.size) + 
  (fit_max$coefficients[[2]]*mean(dat$sqrt.size))

dat3 <- gather(dat,value = trip_Max,
               key = type,
               Maxdist,
               Maxdist_corrected_for_colsize)

m3 <- ggplot(dat3) + 
  geom_ribbon(data=newdat, aes(x=Lat, ymin = lwrMax, ymax = uprMax), 
              fill = "grey", alpha = 0.2) +  
  geom_line(data = dat3, aes(x=Lat, y = trip_Max, group=label),colour="lightgrey",linetype = 2) +
  geom_line(data=newdat,aes(x=Lat, y=newMax)) +  
  geom_point(data=dat3,aes(x=Lat, y = trip_Max, group=label, colour=type),size=3) +
  #geom_text_repel(data = dat, aes(x=Lat, y = Maxdist_corrected_for_colsize, label = label)) +
  scale_color_manual(values=c("lightgrey","black"))+
  ms_theme+
  theme(axis.title = element_blank()) +
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0), limits=c(0, Max_max)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c((min(dat$Lat)-border), 
                              (max(dat$Lat)+border))); m3 
plot_grid(d3, m3, ncol=2)

#Predict sqrt.size
border <- 11.74734
newdat <- tidyr::expand(data = dat,
                        sqrt.size = seq((min(sqrt.size)-border), (max(sqrt.size)+border),length=100),
                        Lat = mean(Lat))

summary(fit_dur)
pred <- predict(fit_dur, newdat, re.form=NA, se = TRUE)
pred <- as.data.frame(pred)
newdat$newDur <- pred$link
newdat$lwrDur <- pred$link - (1.96*pred$SE)
newdat$uprDur <- pred$link + (1.96*pred$SE)

pred <- predict(fit_max, newdat, re.form=NA, se = TRUE)
pred <- as.data.frame(pred)
newdat$newMax <- pred$link
newdat$lwrMax <- pred$link - (1.96*pred$SE)
newdat$uprMax <- pred$link + (1.96*pred$SE)

d2 <- ggplot(newdat, aes(x=sqrt.size, y=newDur)) + 
  geom_ribbon(data=newdat, aes(ymin = lwrDur, ymax = uprDur), 
              fill = "grey", alpha = 0.3) +  
  geom_point(data = dat, aes(x=sqrt.size, y = Dur),size=3) +
  #geom_text_repel(data = dat, aes(x=sqrt.size, y = Dur, label = label)) +
  geom_line()+
  ms_theme+
  theme(axis.title = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, Dur_max)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(0,(max(dat$sqrt.size)+border))); d2

m2 <- ggplot(newdat, aes(x=sqrt.size, y=newMax)) + 
  geom_ribbon(data=newdat, aes(ymin = lwrMax, ymax = uprMax), 
              fill = "grey", alpha = 0.4) +
  geom_point(data = dat, aes(x=sqrt.size, y = Maxdist),size=3) +
  geom_line() +
  ms_theme+
  theme(axis.title = element_blank()) +
  #geom_text_repel(data = dat, aes(x=sqrt.size, y = Maxdist, label = label)) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, Max_max)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits=c(0,(max(dat$sqrt.size)+border))); m2

plot_grid(d2, m2, ncol=2)



#plot grids
plot_grid(d2, d3, m2, m3, ncol=2)

png(paste0(folder,
           "newdata/model_plot_size_cor_nolabels_capesm.png"),
    width=3200,height=3200,res=300)
plot_grid(d2, d3, m2, m3, ncol=2)
dev.off()

# dur v maxdist ####

moddat<-dat[,c("Lat","sqrt.size","Dur_original","Maxdist")]
design<-svydesign(ids=~0, fpc=rep(60,nrow(moddat)), data = moddat)


#Predict dist
border <- 1
newdat <- tidyr::expand(data = dat,
                        Dur_original = seq((min(Dur_original)-border), 
                                           (max(Dur_original)+border),length=100))

fit_dur <- svyglm(Maxdist ~ Dur_original, design = design)
fit2 <- svyglm(Maxdist ~ 1, design = design)
anova(fit_dur,fit2)
summary(fit_dur)
summ(fit_dur)
pred <- predict(fit_dur, newdat, re.form=NA, se = TRUE)
pred <- as.data.frame(pred)

newdat$newMax <- pred$link
newdat$lwrMax <- pred$link - (1.96*pred$SE)
newdat$uprMax <- pred$link + (1.96*pred$SE)


p1 <- ggplot(dat2) + 
  geom_ribbon(data=newdat, aes(x=Dur_original, ymin = lwrMax, ymax = uprMax), 
              fill = "grey", alpha = 0.2) +  
  geom_line(data=newdat,aes(x=Dur_original, y=newMax)) +  
  geom_point(data=dat,aes(x=Dur_original, y = Maxdist),size=3) +
  geom_text_repel(data = dat, aes(x=Dur_original, y = Maxdist, label = label)) +
  ms_theme+
  ylab("Maximum distance from the colony (km)") +
  xlab("Trip duration (h)")  +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 205)) +
  scale_x_continuous(expand = c(0, 0)) ; p1

png("dur_maxdist_capesm.png",
    width=2400,height=2400,res=300)
p1
dev.off()

#dat <- subset(dat,GPS_n_birds > 9) ####
dat <- subset(dat,GPS_n_birds > 9) 

dat <- dat[,c("Lat","Lon","sqrt.size","Dur","Maxdist","GPS_n_birds","Lon")]
design <- svydesign(ids=~0, fpc=rep(60,nrow(dat)), data = dat)

#dur 
fit1<-svyglm(Dur ~ sqrt.size + Lat, design = design)
fit2<-svyglm(Dur ~ sqrt.size , design = design)
fit3<-svyglm(Dur ~ Lat , design = design)

anova(fit1, fit2) #lat effect
anova(fit1, fit3) #size effect

fitll<-lm(Dur ~ sqrt.size + Lat, data = dat)
#plot(fitll)
plot(simulateResiduals(fitll))

summary(fit1) 

summ(fit1) #size+lat
summ(fit2) #size
summ(fit3) #lat

#all
0.67  #adj r2
0.67-0.34 #delta size
0.67-0.57 #delta lat

confint(fit1)

par(mfrow=c(2,2));plot(fit1);par(mfrow=c(1,1))

#max dist 
fit1<-svyglm(Maxdist ~ sqrt.size + Lat, design = design)
fit2<-svyglm(Maxdist ~ sqrt.size , design = design)
fit3<-svyglm(Maxdist ~ Lat , design = design)

anova(fit1, fit2)
anova(fit1, fit3)

fitll<-lm(Maxdist ~ sqrt.size + Lat, data = dat)
#plot(fitll)
plot(simulateResiduals(fitll))
par(mfrow=c(2,2));plot(fit1);par(mfrow=c(1,1))

summary(fit1) 
summ(fit1)
summ(fit2)
summ(fit3)

#>9
0.73  #adj r2
0.73-0.72 #delta lat
0.73-0.16 #delta size

confint(fit1)

