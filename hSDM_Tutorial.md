Introduction to hSDM
====================

Objectives
----------

-   Use opportunistic species occurrence data for occupancy modelling
-   Use `hSDM` R package to fit hierarchical distribution model
-   Compare output from models built with interpolated and satellite-derived environmental data

> One could easily teach an full semester course on the use of this package and modelling framework. Today we will provide only a quick introduction.

This script is available:

-   [MD format with images/plots](https://github.com/adammwilson/hSDM_Tutorial/blob/master/hSDM_Tutorial.md)
-   [Plain text (.R) with commented text](https://raw.githubusercontent.com/adammwilson/hSDM_Tutorial/master/hSDM_Tutorial.R)

You are welcome to follow along on the screen or download the R script and run it on your computer.

The data used below are available in a [public dropbox folder](https://www.dropbox.com/sh/2q0k7qn5rxz0bis/AAB42fVn-s4Teqynrs6rgzR3a?dl=0), though they will be downloaded using code below. You may want to download them now by visiting

Species Distribution Modeling
=============================

Two important problems which can bias model results:

1.  imperfect detections ("false absences")
2.  spatial correlation of the observations.

hSDM R Package
--------------

Developed by [Ghislain Vieilledent](mailto:ghislain.vieilledent@cirad.fr) with Cory Merow, Jérôme Guélat, Andrew M. Latimer, Marc Kéry, Alan E. Gelfand, Adam M. Wilson, Frédéric Mortier & John A. Silander Jr

-   User-friendly statistical functions to overcome limitations above.
-   Developed in a hierarchical Bayesian framework.
-   Call a Metropolis-within-Gibbs algorithm (coded in C) to estimate model parameters and drastically the computation time compared to other methods (e.g. \~2-10x faster than OpenBUGS).

### Software for modeling species distribution including imperfect detection.

![hSDM](assets/DetectionMods.png)

The problem of imperfect detection
----------------------------------

Site-occupancy models (MacKenzie et al., 2002, *aka* zero inflated binomial (ZIB) models) for presence-absence data and Nmixture models (Royle, 2004) or zero inflated Poisson (ZIP) models for abundance data (Flores et al., 2009), were developed to solve the problems created by imperfect detection.

* * * * *

Example application
===================

Load libraries
--------------

``` {.r}
library(hSDM)
library(ggplot2)
library(rasterVis)
library(raster)
library(maptools)
library(dplyr)

library(RCurl)  # Downloading from DropBox

library(coda)  # Summarizing model output

library(doParallel)  #running models in parallel
ncores=2  # number of processor cores you would like to use
registerDoParallel(ncores)
```

If you don't have the packages above, install them in the package manager or by running `install.packages("doParallel")`.

* * * * *

Example Species: *Montane Woodcreeper* (*Lepidocolaptes lacrymiger*)
--------------------------------------------------------------------

<img src="assets/Lepidocolaptes_lacrymiger.jpg" alt="Lepidocolaptes_lacrymiger Photo" width="250px" />

<br><span style="color:grey; font-size:1em;">Figure from [hbw.com](http://www.hbw.com/species/montane-woodcreeper-lepidocolaptes-lacrymiger) </span>

> This species has a large range, occurring from the coastal cordillera of Venezuela along the Andes south to south-east Peru and central Bolivia. [birdlife.org](http://www.birdlife.org/datazone/speciesfactsheet.php?id=31946)

<img src="assets/Lepidocolaptes_lacrymiger_range.png" alt="Lepidocolaptes_lacrymiger Photo" width="200px" />

<br><span style="color:grey; font-size:1em;">Data via [MOL.org](http://map.mol.org/maps/Lepidocolaptes%20lacrymiger) </span>

Set species name:

``` {.r}
sp="Lepidocolaptes_lacrymiger"
```

Data Directory
--------------

The script below will download data to the directory specified below. Feel free to change it as desired.

``` {.r}
## set working directory and path where you want to download data
## you may want to change this to your home directory or elsewhere
datadir=paste0("hSDM_Tutorialdata/")
if(!file.exists(datadir)) dir.create(datadir,recursive = T)

setwd(datadir)
```

Extract species 'expert range' via MOL.
---------------------------------------

``` {.r}
fExpertRange=paste0(sp,".shp")
if(!file.exists(fExpertRange)){
  download.file(paste0("http://mol.cartodb.com/api/v2/sql?",
                     "q=SELECT%20ST_TRANSFORM(the_geom_webmercator,4326)%20as%20the_geom,",
                     "%20seasonality%20FROM%20get_tile('jetz','range','",
                     paste(strsplit(sp,"_")[[1]],collapse="%20"),
                     "','jetz_maps')&format=shp&filename=",sp),
              destfile=sub("shp","zip",fExpertRange))
  unzip(sub("shp","zip",fExpertRange))
}
```

> Full documentation and release of the MOL API in the works.

Load the expert range from the downloaded shapefile.

``` {.r}
reg=readShapePoly(fExpertRange)
## extract bounding box of Expert Range
ereg=extent(reg)
## adjust bbox if desired
ereg@xmin=-81.4
```

Query eBird data contained in MOL
---------------------------------

-   Find all observations of our species
-   Find all unique observation locations for any species limited to bounding box of expert range
-   Filter to where observer indicated recording all observed species (`all_species_reported='t'`)
-   Filter to lists that do not correspond to an observation of our species

> The best method for selecting data to use for *non-detections* is very case and dataset specific.

Metadata for eBird[1] is [available here](http://ebird.org/ebird/data/download)

For this example we'll use data that has been precompiled using the criteria above.

Download species occurrence data
--------------------------------

We've made this exampled dataset available *via* the DropBox links below. If you have the `RCurl` package installed, the following commands should run. If these do not work, [you can also download these datasets from here](https://www.dropbox.com/sh/2q0k7qn5rxz0bis/AAB42fVn-s4Teqynrs6rgzR3a?dl=0).

``` {.r}
fspData="Lepidocolaptes_lacrymiger_points_env.csv"

if(!file.exists(fspData)) {
    URL <- paste0("https://www.dropbox.com/s/1jsyl0xkxnzk20k/",
                  "Lepidocolaptes_lacrymiger_points_env.csv?dl=1")

    download.file(URL,
              destfile=fspData,
              method='curl',extra='-L')
}

spd_all=read.csv(fspData)
```

Check out the data structure:

``` {.r}
kable(head(spd_all[,-1]))
```

|latitude|longitude|observation\_date|presence|effort\_distance\_km|duration\_minutes|effort\_area\_ha|lon|lat|cells|PPTJAN|PPTJUL|PPTSEAS|MAT|ALT|CLDJAN|CLDJUL|CLDSEAS|cell\_lon|cell\_lat|
|-------:|--------:|:----------------|-------:|-------------------:|----------------:|---------------:|--:|--:|----:|-----:|-----:|------:|--:|--:|-----:|-----:|------:|--------:|--------:|
|3.518712|-76.71831|1989-07-06|1|NA|NA|NA|-76.71831|3.518712|2044450|-0.5854745|-0.6289205|-0.6954783|-0.4205597|0.6166451|0.9091247|1.0589650|-1.4256417|-76.72083|3.520833|
|-10.553145|-75.31591|2004-05-05|1|NA|NA|NA|-75.31591|-10.553145|5780687|-0.1790140|-0.9208561|0.1404697|-1.2675758|1.5233957|1.3211621|1.0634702|-1.2863199|-75.31250|-10.554167|
|4.647434|-75.64572|2004-04-06|1|2.000|150|NA|-75.64572|4.647434|1745959|-0.0011875|-0.3512257|-0.6954783|-0.8191555|0.8737833|1.1372553|1.0013975|-1.1757711|-75.64583|4.645833|
|-13.053392|-71.53610|1994-09-03|1|NA|NA|NA|-71.53610|-13.053392|6444740|1.6839302|-0.4722722|-0.1502948|-0.3541270|0.5413645|1.2025126|0.9788711|-0.7154036|-71.53750|-13.054167|
|-16.193000|-67.88500|1999-09-23|1|1.931|390|NA|-67.88500|-16.193000|7279102|0.2189787|-1.0276618|0.7219987|-1.0184534|1.3753721|1.3373415|-0.0508354|0.2674206|-67.88750|-16.195833|
|8.603094|-71.24751|1985-04-18|1|16.093|420|NA|-71.24751|8.603094|695787|-1.1612936|-0.3369850|-0.1139493|-1.5665226|1.4303523|1.1361767|1.3207717|-1.2469464|-71.24583|8.604167|

Explore observer effort: sampling duration, distance travelled, and area surveyed.

``` {.r}
ggplot(spd_all,aes(
  y=duration_minutes/60,
  x=effort_distance_km,
  colour=presence==1,
  order=as.factor(presence)))+
  scale_colour_manual(values=c("grey","red"))+
  ylab("Sampling Duration\n(hours)")+
  xlab("Sampling Distance (km)")+
  labs(colour = "Observed\nPresence")+
  geom_point()+scale_x_log10()
```

![](hSDM_Tutorial_files/figure-markdown_github/spdDurationPlot-1.png)

Also note that there are many records with missing duration and distance values.

``` {.r}
table("Duration"=!is.na(spd_all$duration_minutes),
      "Distance/Area"=!is.na(spd_all$effort_distance_km)|
        !is.na(spd_all$effort_area_ha))
```

    ##         Distance/Area
    ## Duration FALSE  TRUE
    ##    FALSE 13280    93
    ##    TRUE   6713 26842

> For this exercise, we'll simply remove points with large or unknown spatial uncertainty. Incorporating this spatial uncertainty into distribution models is an active area of research.

Filter the data below thresholds specified above.

``` {.r}
cdur=4*60  # Duration in minutes
cdis=5     # Distance in km
care=500   # Area in Hectares
 
spd=filter(spd_all,
           duration_minutes<=cdur&
          (effort_distance_km<=cdis|effort_area_ha<=care))
```

Convert to a spatialDataFrame to faciliate linking with georeferenced environmental data.

``` {.r}
coordinates(spd)=c("longitude","latitude")
projection(spd)="+proj=longlat +datum=WGS84 +ellps=WGS84"
spd@data[,c("lon","lat")]=coordinates(spd)   
```

### Load coastline from maptools package for plotting.

``` {.r}
coast <- map_data("world",
                  xlim=c(ereg@xmin-1,ereg@xmax+1),
                  ylim=c(ereg@ymin-1,ereg@ymax+1))

ggcoast=geom_path(data=coast,
                  aes(x=long,y=lat,group = group),lwd=.1)
```

Available Species Data
----------------------

``` {.r}
ggplot(spd@data,aes(y=lat,x=lon))+
  ggcoast+ 
  geom_path(data=fortify(reg),
            aes(y=lat,x=long,group=piece),
            col="green")+
  geom_point(data=filter(spd@data,presence==0),
             aes(x=lon,y=lat),pch=1,
             col="black",cex=.8,lwd=2,alpha=.3)+
  geom_point(data=filter(spd@data,presence==1),
             aes(x=lon,y=lat),pch=3,
             col="red",cex=2,lwd=3,alpha=1)+
  ylab("Latitude")+xlab("Longitude")+
  coord_equal()+
  xlim(c(min(spd$lon),max(spd$lon)))+
  ylim(c(min(spd$lat),max(spd$lat)))
```

![](hSDM_Tutorial_files/figure-markdown_github/spdPlot-1.png)

* * * * *

Environmental Data
------------------

We've also pre-compiled environmental data for the region and made it available in the shared DropBox folder.

-   **PPTJAN**: Mean January Precipitation (mm, WorldClim)
-   **PPTJUL**: Mean January Precipitation (mm, WorldClim)
-   **PPTSEAS**: Precipitation Seasonality (WorldClim)
-   **MAT**: Mean Annual Temperature (C, WorldClim)
-   **ALT**: Elevation (m, WorldClim)
-   **CLDJAN**: Mean January Cloud Frequency (1000s %, Wilson&Jetz)
-   **CLDJUL**: Mean July Cloud Frequency (1000s %, Wilson&Jetz)
-   **CLDSEAS**: Cloud Seasonality (1000s %, Wilson&Jetz)

Download a single geotif with 8 bands corresponding to the data above for a subset of the species' range.

``` {.r}
fenvdata="Lepidocolaptes_lacrymiger_env_scaled.tif"

if(!file.exists(fenvdata)) {

URL <- paste0("https://www.dropbox.com/s/7i5hl3gv53l8m4v/",
           "Lepidocolaptes_lacrymiger_env_scaled_small.tif?dl=1")

download.file(URL,
              destfile=fenvdata,
              method='curl',extra='-L')
}
 
env=stack(fenvdata)
gain(env)=0.01
names(env)=c("PPTJAN","PPTJUL","PPTSEAS","MAT","ALT","CLDJAN","CLDJUL","CLDSEAS")
```

### Visualize the environmental data

``` {.r}
## set plotting limits using expert range above

gx=xlim(extent(env)@xmin,extent(env)@xmax)
gy=ylim(extent(env)@ymin,extent(env)@ymax)


## Environmental data
gplot(env)+
  geom_raster(aes(fill=value)) + 
  facet_wrap(~variable,nrow=2) +
  scale_fill_gradientn(
    colours=c('blue','white','red','darkred'),
    breaks=c(-3,0,3,6),na.value="transparent")+
  ylab("")+xlab("")+labs(fill = "Standardized\nValue")+
    geom_path(data=fortify(reg),aes(y=lat,x=long,group=piece),fill="green",col="black")+
  ggcoast+gx+gy+coord_equal()
```

![](hSDM_Tutorial_files/figure-markdown_github/plotEnv-1.png) \#\#\# Covariate correlation Scatterplot matrix of the available environmental data.

``` {.r}
splom(env,varname.cex=2)
```

![](hSDM_Tutorial_files/figure-markdown_github/envCor-1.png)

### Generate `data.frame` for model fitting

First we need to 'grid' the point-level species observations to match the environmental data.

``` {.r}
fdata=spd@data%>%
  group_by(cell_lon,cell_lat,cells)%>%
  summarize(
            presences=sum(presence),  # Number of observations
            trials=sum(presence==0),  # Number of trials in each cell
            PPTJAN=mean(PPTJAN),
            PPTJUL=mean(PPTJUL),
            PPTSEAS=mean(PPTSEAS),
            MAT=mean(MAT),
            ALT=mean(ALT),
            CLDJAN=mean(CLDJAN),
            CLDJUL=mean(CLDJUL),
            CLDSEAS=mean(CLDSEAS))

## Due to oppotunistic observations of the species, a few grid cells
## have more observations than trials.  Could remove them or adjust 
## trials up in these cells
fdata$trials=ifelse(fdata$presences>fdata$trials,
                     fdata$presences,fdata$trials)

## omit rows with missing data (primarily ocean pixels)
fdata=na.omit(fdata)
```

Then transform the full gridded dataset into a `data.frame` with associated environmental data for predicting across space.

``` {.r}
pdata=cbind.data.frame(
  coordinates(env),
  values(env))


## omit rows with missing data (primarily ocean pixels)
pdata=na.omit(pdata)

kable(head(pdata))
```

||x|y|PPTJAN|PPTJUL|PPTSEAS|MAT|ALT|CLDJAN|CLDJUL|CLDSEAS|
|---|--:|--:|-----:|-----:|------:|--:|--:|-----:|-----:|------:|
|859|-74.24583|11.22083|-1.42|-0.81|0.97|0.84|-0.62|-3.09|-0.55|1.62|
|860|-74.23750|11.22083|-1.42|-0.80|0.94|0.82|-0.61|-3.09|-0.53|1.72|
|861|-74.22917|11.22083|-1.40|-0.70|0.90|0.72|-0.53|-3.11|-0.41|1.83|
|862|-74.22083|11.22083|-1.38|-0.60|0.94|0.64|-0.47|-3.09|-0.18|2.00|
|863|-74.21250|11.22083|-1.41|-0.78|0.90|0.82|-0.61|-3.12|0.11|2.30|
|864|-74.20417|11.22083|-1.41|-0.77|0.90|0.82|-0.61|-3.10|0.27|2.47|

This table is similar to the data available from the "Annotate" function in MOL, with the exception that it contains the point data aggregated to the resolution of the environmental data.

Model Comparison
----------------

Let's compare two models, one using interpolated precipitation and the other using satellite-derived cloud data.

``` {.r}
# Set number of chains to fit.
mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c( "Precipitation",
          "Cloud"))

kable(mods)
```

|model|formula|name|
|:----|:------|:---|
|m1|\~PPTJAN+PPTJUL+PPTSEAS+MAT|Precipitation|
|m2|\~CLDJAN+CLDJUL+CLDSEAS+MAT|Cloud|

Specify model run-lengths.

``` {.r}
burnin=1000
mcmc=1000
thin=1
```

Fit the models
--------------

Both site-occupancy or ZIB models (with `hSDM.siteocc()` or `hSDM.ZIB()` functions respectively) can be used to model the presence-absence of a species taking into account imperfect detection.

The site-occupancy model can be used in all cases but can be less convenient and slower to fit when the repeated visits at each site are made under the same observation conditions. While this is likely not true in this situation (the observations occurred in different years, etc.), we'll use the simpler model today. For more information about the differences, see the hSDM Vignette Section 4.3.

### Example: `hSDM.ZIB`

The model integrates two processes, an ecological process associated to the presence or absence of the species due to habitat suitability and an observation process that takes into account the fact that the probability of detection of the species is less than one.

If the species has been observed at least once during multiple visits, we can assert that the habitat at this site is suitable. And the fact that the species can be unobserved at this site is only due to imperfect detection.

**Ecological process:**

![](assets/M1.png)

**Observation process:**

![](assets/M2.png)

In this example we'll assume a spatially constant p(observation|presence), but it's also possible to put in covariates for this parameter.

Run the model
-------------

``` {.r}
results=foreach(m=1:nrow(mods)) %dopar% { 
  ## if foreach/doParallel are not installed, you can use this line instead
  # for(m in 1:nrow(mods)) { 
  tres=hSDM.ZIB(
    suitability=as.character(mods$formula[m]),  #Formula for suitability
    presences=fdata$presences,    # Number of Observed Presences
    observability=~1,             # Formula for p(observation|presence)
    trials=fdata$trials,          # Number of Trials
    data=fdata,                   # Covariates for fitting model
    suitability.pred=pdata,        # Covariates for prediction 
    mugamma=0, Vgamma=1.0E6,      # Priors on Gamma
    gamma.start=0,                # Gamma initial Value
    burnin=burnin, mcmc=mcmc, thin=thin,  # MCMC parameters
    beta.start=0,                 # Initial values for betas
    mubeta=0, Vbeta=1.0E6,        # Priors on Beta
    save.p=0,                     # Don't save full posterior on p
    verbose=1,                    # Report progress
    seed=round(runif(1,0,1e6)))   # Random seed
  tres$model=mods$formula[m]      # Add model formula to result
  tres$modelname=mods$name[m]     # Add model name to result
  return(tres)
  }
```

Summarize posterior parameters
------------------------------

The model returns full posterior distributions for all model parameters. To summarize them you need to choose your summary metric (e.g. mean/median/quantiles).

``` {.r}
params=foreach(r1=results,.combine=rbind.data.frame)%do% {
  data.frame(model=r1$model,
             modelname=r1$modelname,
             parameter=colnames(r1$mcmc),
             mean=summary(r1$mcmc)$statistics[,"Mean"],
             sd=summary(r1$mcmc)$statistics[,"SD"],
             median=summary(r1$mcmc)$quantiles[,"50%"],
             HPDinterval(mcmc(as.matrix(r1$mcmc))),
             RejectionRate=rejectionRate(r1$mcmc))}
```

``` {.r}
## plot it
ggplot(params[!grepl("Deviance*",rownames(params)),],
       aes(x=mean,y=parameter,colour=modelname))+
  geom_errorbarh(aes(xmin=lower,xmax=upper,height=.1),lwd=1.5)+
  geom_point(size = 3)+xlab("Standardized Coefficient")+
  ylab("Parameter")+
  theme(legend.position="bottom")
```

![](hSDM_Tutorial_files/figure-markdown_github/PlotPosteriors-1.png)

### Detection probability

The model uses repeat obserations within cells to estimate the probabiliy observation given that the species was present.

``` {.r}
pDetect <-   params[params$parameter=="gamma.(Intercept)",
                    c("modelname","mean")]
pDetect$delta.est <- inv.logit(pDetect$mean)
colnames(pDetect)[2]="gamma.hat"
kable(pDetect,row.names=F)
```

|modelname|gamma.hat|delta.est|
|:--------|--------:|--------:|
|Precipitation|-0.9347880|0.2819543|
|Cloud|-0.8972593|0.2896140|

> How does this change if you add environmental covariates to the observability regression?

Predictions for each cell
-------------------------

``` {.r}
pred=foreach(r1=results,.combine=stack)%dopar% {
  tr=rasterFromXYZ(cbind(x=pdata$x,
                         y=pdata$y,
                         pred=r1$prob.p.pred))
  names(tr)=r1$modelname    
  return(tr)
  }
```

Compare the model predictions

``` {.r}
predscale=scale_fill_gradientn(values=c(0,.5,1),colours=c('white','darkgreen','green'),na.value="transparent")

gplot(pred)+geom_raster(aes(fill=value)) +
  facet_wrap(~ variable) +
  predscale+
  coord_equal()+
  geom_path(data=fortify(reg),
            aes(y=lat,x=long,group=piece),
            fill="red",col="red")+
  geom_point(data=filter(spd@data,presence==0),
             aes(x=lon,y=lat,fill=1),pch=1,col="black",cex=.8,lwd=2,alpha=.3)+
  geom_point(data=filter(spd@data,presence==1),
             aes(x=lon,y=lat,fill=1),pch=3,col="red",cex=2,lwd=3,alpha=1)+
  ggcoast+gx+gy+ylab("Latitude")+xlab("Longitude")+
  labs(col = "p(presence)")+
  coord_equal()
```

![](hSDM_Tutorial_files/figure-markdown_github/plotmodel-1.png)

Additional Models in hSDM
-------------------------

`*.icar` The `*.icar` functions in `hSDM` add *spatial effects* to the model as well, accounting for spatial autocorrelation of species occurrence.

`hSDM.binomial` & `hSDM.binomial.iCAR` Simple and spatial binomial model (perfect detection).

`hSDM.ZIB` & `hSDM.ZIB.iCAR` & `hSDM.ZIB.iCAR.alteration` Zero-inflated Binomial (example we used today).

`hSDM.ZIP` & `hSDM.ZIP.iCAR` & `hSDM.ZIP.iCAR.alteration` Zero-inflated Poisson (Abundance data with imperfect detection).

`hSDM.siteocc` & `hSDM.siteocc.iCAR` Incorporates temporally varying environment to account for changing observation conditions.

`hSDM.poisson` & `hSDM.poisson.iCAR` Simple and spatial poisson model for species abundance (perfect detection).

`hSDM.Nmixture` & `hSDM.Nmixture.iCAR` Poisson model for abundance with imperfect detection.

Looking forward
---------------

-   Incorporate multiple scales of observations (e.g. points & polygons)
-   Account directly for spatial uncertainties in point observations
-   Time-varying covariates with `hSDM.siteocc` or similar

[1] M. Arthur Munson, Kevin Webb, Daniel Sheldon, Daniel Fink, Wesley M. Hochachka, Marshall Iliff, Mirek Riedewald, Daria Sorokina, Brian Sullivan, Christopher Wood, and Steve Kelling. The eBird Reference Dataset, Version 5.0. Cornell Lab of Ornithology and National Audubon Society, Ithaca, NY, January 2013.
