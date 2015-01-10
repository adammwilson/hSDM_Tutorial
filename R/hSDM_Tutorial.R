#' ---
#' title: "Overview of hSDM"
#' author: "Adam M. Wilson"
#' date: "January 6, 2015"
#' output:
#'   md_document:
#'     variant: markdown_github
#' ---
#' 
## ----, echo=FALSE, warning=FALSE-----------------------------------------
library(knitr)
library(rmarkdown)

opts_chunk$set(cache=TRUE,
               root.dir="/Users/adamw/repos/hSDM_Tutorial",
               warning=FALSE,
               message=F,
               fig.width=12,
               fig.height=15)

# purl("R/hSDM_Tutorial.Rmd","R/hSDM_Tutorial.R",documentation=2)
# rmarkdown::render("R/hSDM_Tutorial.Rmd", "all")


#' # Introduction to hSDM
#' 
#' ## Objectives
#' 
#' * Use opportunistic species occurrence data for occupancy modelling
#' * Use `hSDM` R package to fit hierarchical distribution model
#' * Compare output from models built with interpolated and satellite-derived environmental data
#' 
#' > One could easily teach an full semester course on the use of this package and modelling framework.  Today we will provide only a quick introduction.
#' 
#' This script is available:
#' 
#'   * [hSDM Tutorial on GitHub](https://github.com/adammwilson/hSDM_Tutorial)
#'   * [HTML format with images/plots](https://rawgit.com/adammwilson/hSDM_Tutorial/master/R/hSDM_Tutorial.html)
#'   * [Plain text (.R) with commented text](https://raw.githubusercontent.com/adammwilson/hSDM_Tutorial/master/R/hSDM_Tutorial.R)
#' 
#' The data used below are available in a [public dropbox folder](https://www.dropbox.com/sh/2q0k7qn5rxz0bis/AAB42fVn-s4Teqynrs6rgzR3a?dl=0), though they will be downloaded using code below.  
#' 
#' # Species Distribution Modeling
#' 
#' Two major problems which can bias model results: 
#' 
#' 1. imperfect detections ("false absences")
#' 2. spatial correlation of the observations.
#' 
#' ## hSDM R Package
#' 
#' Developed by [Ghislain Vieilledent](mailto:ghislain.vieilledent@cirad.fr) with Cory Merow, Jérôme Guélat, Andrew M. Latimer, Marc Kéry, Alan E. Gelfand, Adam M. Wilson, Frédéric Mortier & John A. Silander Jr
#' 
#' * User-friendly statistical functions to overcome limitations above.
#' * Developed in a hierarchical Bayesian framework. 
#' * Call a Metropolis-within-Gibbs algorithm (coded in C) to estimate model parameters and drastically the computation time compared to other methods (e.g. ~2-10x faster than OpenBUGS).
#' 
#' ### Software for modeling species distribution including imperfect detection.
#' ![hSDM](../assets/DetectionMods.png)
#' 
#' ## The problem of imperfect detection
#' 
#' Site-occupancy models (MacKenzie et al., 2002, _aka_ zero inflated binomial (ZIB) models) for presence-absence data and Nmixture models (Royle, 2004) or zero inflated Poisson (ZIP) models for abundance data (Flores et al., 2009), were developed to solve the problems created by imperfect detection.
#' 
#' _________________
#' 
#' # Example application
#' 
#' ## Load libraries
## ----loadLibraries,message=FALSE,warning=FALSE---------------------------
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

## set light theme for ggplot
theme_set(theme_light()+
            theme(text=element_text(size = 38)))


#' If you don't have the packages above, install them in the package manager or by running `install.packages("doParallel")`. 
#' 
#' 
#' __________
#' 
#' ## Example Species: *Montane Woodcreeper* (_Lepidocolaptes lacrymiger_)
#' 
#' ![Lepidocolaptes_lacrymiger Photo](../assets/Lepidocolaptes_lacrymiger.jpg)
#' <br><span style="color:grey; font-size:1em;">Figure from [hbw.com](http://www.hbw.com/species/montane-woodcreeper-lepidocolaptes-lacrymiger) </span>
#' 
#' > This species has a large range, occurring from the coastal cordillera of Venezuela along the Andes south to south-east Peru and central Bolivia. [birdlife.org](http://www.birdlife.org/datazone/speciesfactsheet.php?id=31946)
#' 
#' ![Lepidocolaptes_lacrymiger Data](../assets/Lepidocolaptes_lacrymiger_range.png)
#' <br><span style="color:grey; font-size:1em;">Data via [MOL.org](http://map.mol.org/maps/Lepidocolaptes%20lacrymiger) </span>
#' 
#' 
#' Set species name:
## ------------------------------------------------------------------------
sp="Lepidocolaptes_lacrymiger"

#' 
#' ## Data Directory
#' The script below will download data to the directory specified below.  Feel free to change it as desired.
## ------------------------------------------------------------------------
## set path to data folder
datadir=paste0("../data/",sp,"/")
if(!file.exists(datadir)) dir.create(datadir,recursive = T)

#' 
#' 
#' ## Extract species 'expert range' via MOL.
#' 
## ----getRange------------------------------------------------------------
fExpertRange=paste0(datadir,"/",sp,".shp")
if(!file.exists(fExpertRange)){
  download.file(paste0("http://mol.cartodb.com/api/v2/sql?",
                     "q=SELECT%20ST_TRANSFORM(the_geom_webmercator,4326)%20as%20the_geom,",
                     "%20seasonality%20FROM%20get_tile('jetz','range','",
                     paste(strsplit(sp,"_")[[1]],collapse="%20"),
                     "','jetz_maps')&format=shp&filename=",sp),
              destfile=sub("shp","zip",fExpertRange))
  unzip(sub("shp","zip",fExpertRange),exdir=datadir)
}

#' 
#' > Full documentation and release of the MOL API in the works.
#' 
#' Load the expert range from the downloaded shapefile.
## ----readRange-----------------------------------------------------------
reg=readShapePoly(fExpertRange)
## extract bounding box of Expert Range
ereg=extent(reg)
## adjust bbox if desired
ereg@xmin=-81.4

#' 
#' 
#' ## Query eBird data contained in MOL
#' 
#' * Find all observations of our species
#' * Find all unique observation locations for any species limited to bounding box of expert range
#' * Filter to where observer indicated recording all observed species (`all_species_reported='t'`)
#' * Filter to lists that do not correspond to an observation of our species
#' 
#' > The best method for selecting data to use for _non-detections_ is very case and dataset specific.
#' 
#' Metadata for eBird^[M. Arthur Munson, Kevin Webb, Daniel Sheldon, Daniel Fink, Wesley M. Hochachka, Marshall Iliff, Mirek Riedewald, Daria Sorokina, Brian Sullivan, Christopher Wood, and Steve Kelling. The eBird Reference Dataset, Version 5.0. Cornell Lab of Ornithology and National Audubon Society, Ithaca, NY, January 2013.] is [available here](http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz)
#' 
#' For this example we'll use data that has been precompiled using the criteria above. If you'd like to see how we compiled these data, [see here](https://github.com/adammwilson/hSDM_Tutorial/blob/master/R/hSDM_DataPrep.md)
#' 
#' ## Download species occurrence data
#' 
#' We've made this exampled dataset available _via_ the DropBox links below.  If you have the `RCurl` package installed, the following commands should run. If these do not work, [you can also download these datasets from here](https://www.dropbox.com/sh/2q0k7qn5rxz0bis/AAB42fVn-s4Teqynrs6rgzR3a?dl=0).
#' 
## ----loadSpd-------------------------------------------------------------
fspData=paste0(datadir,"Lepidocolaptes_lacrymiger_points.csv")

if(!file.exists(fspData)) {
    URL <- paste0("https://www.dropbox.com/s/y9np2fdpw5lg5jw/",
                  "Lepidocolaptes_lacrymiger_points.csv?dl=1")

    download.file(URL,
              destfile=fspData,
              method='curl',extra='-L')
}

spd_all=read.csv(fspData)

#' 
#' Check out the data structure:
## ----headSpd, results='asis'---------------------------------------------
kable(head(spd_all[,-1]))

#' 
#' Explore  observer effort: sampling duration, distance travelled, and area surveyed.
## ----spdDurationPlot-----------------------------------------------------
ggplot(spd_all,aes(
  y=duration_minutes/60,
  x=effort_distance_km,
  colour=presence==1,
  order=as.factor(presence)))+
  ylab("Sampling Duration (hours)")+
  xlab("Sampling Distance (km)")+
  labs(col = "Observed\nPresence")+
  geom_point()+scale_x_log10()

#' 
#' 
#' Also note that there are many records with missing duration and distance values.
## ----spNA----------------------------------------------------------------
table("Duration"=!is.na(spd_all$duration_minutes),
      "Distance/Area"=!is.na(spd_all$effort_distance_km)|
        !is.na(spd_all$effort_area_ha))

#' 
#' > For this exercise, we'll simply remove points with large or unknown spatial uncertainty.  Incorporating this spatial uncertainty into distribution models is an active area of research.
#' 
#' 
#' Filter the data below thresholds specified above.
## ----filterSpd-----------------------------------------------------------

cdur=4*60  # Duration in minutes
cdis=5     # Distance in km
care=500   # Area in Hectares

spd=filter(spd_all,
           duration_minutes<=cdur&
          (effort_distance_km<=cdis|effort_area_ha<=care))

#' 
#' 
#' Convert to a spatialDataFrame to faciliate linking with georeferenced environmental data.
#' 
## ----spSpd---------------------------------------------------------------
coordinates(spd)=c("longitude","latitude")
projection(spd)="+proj=longlat +datum=WGS84 +ellps=WGS84"
spd@data[,c("lon","lat")]=coordinates(spd)  

#' 
#' ### Load coastline from maptools package for plotting.
## ----loadCoast-----------------------------------------------------------
coast <- map_data("world",
                  xlim=c(ereg@xmin-1,ereg@xmax+1),
                  ylim=c(ereg@ymin-1,ereg@ymax+1))
ggcoast=geom_path(data=coast,aes(x=long,y=lat,group = group),lwd=.1)

## set plotting limits using expert range above
gx=xlim(ereg@xmin-1,ereg@xmax+1)
gy=ylim(ereg@ymin-1,ereg@ymax+1)


#' 
#' ## Available Species Data
## ----spdPlot, message=FALSE----------------------------------------------
ggplot(spd@data,aes(y=lat,x=lon))+
  geom_path(data=fortify(reg),aes(y=lat,x=long,group=piece),fill="green",col="green")+
  geom_point(aes(colour=presence>=1,order=as.factor(presence)))+
  ggcoast+gx+gy+ylab("Latitude")+xlab("Longitude")+
  labs(col = "Species\nObserved")+
  coord_equal()


#' 
#' ______________________
#' 
#' ## Environmental Data
#' 
#' We've also pre-compiled environmental data for the region and made it available in the shared DropBox folder.
#' 
#' * **PPTJAN**: Mean January Precipitation (mm, WorldClim)
#' * **PPTJUL**: Mean January Precipitation (mm, WorldClim)
#' * **PPTSEAS**: Precipitation Seasonality (WorldClim)
#' * **MAT**: Mean Annual Temperature (C, WorldClim)
#' * **ALT**: Elevation (m, WorldClim)
#' * **CLDJAN**:  Mean January Cloud Frequency (1000s %, Wilson&Jetz)
#' * **CLDJUL**:  Mean July Cloud Frequency (1000s %, Wilson&Jetz)
#' * **CLDSEAS**: Cloud Seasonality (1000s %, Wilson&Jetz)
#' 
#' Download a single geotif with 8 bands corresponding to the data above.
## ----loadEnv-------------------------------------------------------------
fenvdata="data/Lepidocolaptes_lacrymiger.tif"

if(!file.exists(fenvdata)) {
    URL <- paste0("https://www.dropbox.com/s/1vjw848ehyzybl1/",
              "Lepidocolaptes_lacrymiger_env.tif?dl=1")

    download.file(URL,
              destfile=fenvdata,
              method='curl',extra='-L')
}

env=stack(paste0(datadir,sp,"env.tif"))
names(env)=c("PPTJAN","PPTJUL","PPTSEAS","MAT","ALT","CLDJAN","CLDJUL","CLDSEAS")

#' 
#' ### Scale environmental data
#' Scaling covariate data^[subtract the mean and divide by the standard deviation] results in standardized parameter values and also can speed up modeling convergence.  It's possible to _unscale_ the results later if desired.
#' 
## ----scaleEnv------------------------------------------------------------
## Create a 'scaled' version for modeling
senv=raster::scale(env)

#' 
#' ### Visualize the environmental data
#' 
## ----plotEnv-------------------------------------------------------------
## Environmental data
gplot(senv)+
  geom_raster(aes(fill=value)) + 
  facet_wrap(~variable,nrow=2) +
  scale_fill_gradientn(colours=c('blue','white','red','darkred'),breaks=c(-3,0,3,6),na.value="transparent")+
  ylab("")+xlab("")+labs(fill = "Standardized\nValue")+
  ggcoast+gx+gy

#' ### Covariate correlation
#' Scatterplot matrix of the available environmental data.
#' 
## ----envCor--------------------------------------------------------------
splom(senv,varname.cex=2)

#' 
#' ### Generate `data.frame` for model fitting
#' 
#' First we need to 'grid' the point-level species observations to match the environmental data.
#' 
## ----gridSpd-------------------------------------------------------------
## add cell id to facilitate linking points to raster
cell=env[[1]]
raster::values(cell)=1:ncell(cell)
names(cell)="cell"

## rasterize points
presences=rasterize(spd,env,fun="sum",field="presence",background=0)
trials=rasterize(spd,env,fun="count",field="presence",background=0)

#' 
#' Then transform the full gridded dataset into a point-level dataset with the number of `trials`, `presences` and associated environmental data.
#' 
## ----spdReformat---------------------------------------------------------
data=cbind.data.frame(
  coordinates(senv),
  trials=values(trials),
  presences=values(presences),
  cell=values(cell),
  values(senv))

## omit rows with missing data (primarily ocean pixels)
data=na.omit(data)

kable(head(data))

#' 
#' ### Select data for fitting
#' Create 'fitting' dataset limited to locations with at least one `trial`.
#' 
## ----MakeFittingdata-----------------------------------------------------
fdata=data[data$trials>0,]

#' 
#' ## Model Comparison
#' Let's compare two models, one using interpolated precipitation and the other using satellite-derived cloud data.
#' 
## ----modelNames----------------------------------------------------------
# Set number of chains to fit.
mods=data.frame(
  model=c("m1","m2"),
  formula=c("~PPTJAN+PPTJUL+PPTSEAS+MAT",
            "~CLDJAN+CLDJUL+CLDSEAS+MAT"),
  name=c( "Precipitation",
          "Cloud"))

kable(mods)

#' 
#' Specify model run-lengths.  
## ----modelSetup----------------------------------------------------------
burnin=1000
mcmc=1000
thin=1

#' ## Fit the models
#' 
#' Both site-occupancy or ZIB models (with `hSDM.siteocc()` or `hSDM.ZIB()` functions respectively) can be used to model the presence-absence of a species taking into account imperfect detection. 
#' 
#' The site-occupancy model can be used in all cases but can be less convenient and slower to fit when the repeated visits at each site are made under the same observation conditions. While this is likely not true in this situation (the observations occurred in different years, etc.), we'll use the simpler model today.  For more information about the differences, see the hSDM Vignette Section 4.3.  
#' 
#' ### Example: `hSDM.ZIB`
#' The model integrates two processes, an ecological process associated to the presence or absence of the species due to habitat suitability and an observation process that takes into account the fact that
#' the probability of detection of the species is less than one.
#' 
#' If the species has been observed at least once during multiple visits, we can assert that the habitat at this site is suitable. And the fact that the species can be unobserved at this site is only due to imperfect detection.
#' 
#' **Ecological process:**
#' 
#'  ![](../assets/M1.png)
#' 
#' **Observation process:**
#' 
#'  ![](../assets/M2.png)
#' 
#' 
## ----runmodel------------------------------------------------------------
results=foreach(m=1:nrow(mods)) %dopar% { 
  ## if foreach/doParallel are not installed, you can use this line instead
  # for(m in 1:nrow(mods)) { 
  tres=hSDM.ZIB(
    suitability=as.character(mods$formula[m]),  #Formula for suitability
    presences=fdata$presences,    # Number of Observed Presences
    observability=~1,             # Formula for p(observation|presence)
    trials=fdata$trials,          # Number of Trials
    data=fdata,                   # Covariates for fitting model
    suitability.pred=data,        # Covariates for prediction 
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


#' 
#' ## Summarize posterior parameters
#' The model returns full posterior distributions for all model parameters.  To summarize them you need to choose your summary metric (e.g. mean/median/quantiles). 
#' 
## ----SummarizePosteriors-------------------------------------------------

params=foreach(r1=results,.combine=rbind.data.frame)%do% {
  data.frame(model=r1$model,
             modelname=r1$modelname,
             parameter=colnames(r1$mcmc),
             mean=summary(r1$mcmc)$statistics[,"Mean"],
             sd=summary(r1$mcmc)$statistics[,"SD"],
             median=summary(r1$mcmc)$quantiles[,"50%"],
             HPDinterval(mcmc(as.matrix(r1$mcmc))),
             RejectionRate=rejectionRate(r1$mcmc))}

## plot it
ggplot(params[!grepl("Deviance*",rownames(params)),],
       aes(x=mean,y=parameter,colour=modelname))+
  geom_point()+
  geom_errorbarh(aes(xmin=lower,xmax=upper,height=.1))+
  theme(legend.position="bottom")

#' 
#' 
#' ### Detection probability
#' The model uses repeat obserations within cells to estimate the probabiliy observation given that the species was present.  
#' 
## ----pDetect-------------------------------------------------------------
pDetect <-   params[params$parameter=="gamma.(Intercept)",
                    c("modelname","mean")]
pDetect$delta.est <- inv.logit(pDetect$mean)
colnames(pDetect)[2]="gamma.hat"
kable(pDetect,row.names=F)

#' 
#' >  How does this change if you add environmental covariates to the observability regression?
#' 
#' ## Predictions for each cell
## ----pPred---------------------------------------------------------------

pred=foreach(r1=results,.combine=stack)%dopar% {
  tr=rasterFromXYZ(cbind(x=data$x,
                         y=data$y,
                         pred=r1$prob.p.pred))
  names(tr)=r1$modelname    
  return(tr)
  }

#' 
#' Compare the model predictions
## ----plotmodel-----------------------------------------------------------

predscale=scale_fill_gradientn(values=c(0,.5,1),colours=c('white','darkgreen','green'),na.value="transparent")

gplot(pred)+geom_raster(aes(fill=value)) +
  facet_wrap(~ variable) +
  predscale+
  coord_equal()+
  geom_path(data=fortify(reg),
            aes(y=lat,x=long,group=piece),fill="red",col="red")+
  ggcoast+gx+gy+ylab("Latitude")+xlab("Longitude")+
  labs(col = "p(presence)")+
  coord_equal()


#' 
#' ## Additional Models in hSDM
#' 
#' ### `*.icar`
#' The `*.icar` functions in `hSDM` add _spatial effects_ to the model as well, accounting for spatial autocorrelation of species occurrence.  
#' 
#' 
#' ### hSDM.binomial & hSDM.binomial.iCAR
#' Simple and spatial binomial model (perfect detection).
#' 
#' ### hSDM.ZIB & hSDM.ZIB.iCAR & hSDM.ZIB.iCAR.alteration
#' Zero-inflated Binomial (example we used today).
#' 
#' ### hSDM.ZIP & hSDM.ZIP.iCAR & hSDM.ZIP.iCAR.alteration
#' Zero-inflated Poisson (Abundance data with imperfect detection).
#' 
#' ### hSDM.siteocc & hSDM.siteocc.iCAR
#' Incorporates temporally varying environment to account for changing observation conditions.  
#' 
#' ### hSDM.poisson & hSDM.poisson.iCAR
#' Simple and spatial poisson model for species abundance (perfect detection).
#' 
#' ### hSDM.Nmixture & hSDM.Nmixture.iCAR
#' Poisson model for abundance with imperfect detection.
#' 
#' 
#' ## Looking forward
#' 
#' * Incorporate multiple scales of observations (e.g. points & polygons)
#' * Account directly for spatial uncertainties in point observations
#' * Time-varying covariates with `hSDM.siteocc` or similar
