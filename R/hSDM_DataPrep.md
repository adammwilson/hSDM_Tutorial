Data Preparation for hSDM tutorial
==================================

Load libraries
--------------

``` {.r}
library(hSDM)
library(ggplot2)
library(rasterVis)
library(raster)
library(maptools)
library(dplyr)
```

* * * * *

Example Species: *Montane Woodcreeper* (*Lepidocolaptes lacrymiger*)
--------------------------------------------------------------------

![Lepidocolaptes\_lacrymiger Photo](../assets/Lepidocolaptes_lacrymiger.jpg) <br><span style="color:grey; font-size:1em;">Figure from [hbw.com](http://www.hbw.com/species/montane-woodcreeper-lepidocolaptes-lacrymiger) </span>

> This species has a large range, occurring from the coastal cordillera of Venezuela along the Andes south to south-east Peru and central Bolivia. [birdlife.org](http://www.birdlife.org/datazone/speciesfactsheet.php?id=31946)

``` {.r}
sp="Lepidocolaptes_lacrymiger"

## set path to data folder
datadir=paste0("../data/",sp,"/")
```

Query eBird data contained in MOL
---------------------------------

-   Find all observations of our species
-   Find all unique observation locations for any species limited to bounding box of expert range
-   Filter to where observer indicated recording all observed species (`all_species_reported='t'`)
-   Filter to lists that do not correspond to an observation of our species

> The best method for selecting data to use for *non-detections* is very case and dataset specific.

Metadata for eBird[1] is [available here](http://ebirddata.ornith.cornell.edu/downloads/erd/ebird_all_species/erd_western_hemisphere_data_grouped_by_year_v5.0.tar.gz)

Below is an R function that queries the eBird data and summarized as bulleted above. This database is not currently publically accessible, so we're providing the summarized data below.

``` {.r}
# install_github("pingles/redshift-r")
getebird=function(con, sptaxon, region){
  print(paste("Extracting data, this can take a few minutes..."))
    dbGetQuery(conn,paste(   
        "WITH ebird_subset as (SELECT 
            all_species_reported,
            taxonomic_order,
            latitude,
            longitude,
            observation_date,
            sampling_event_identifier,
            group_identifier,
            effort_distance_km,
            effort_area_ha,
            duration_minutes
          FROM ebird
          WHERE latitude BETWEEN ",paste(bbox(region)["y",],collapse=" AND "),"
          AND longitude BETWEEN ",paste(bbox(region)["x",],collapse=" AND "),")
        presence as (SELECT DISTINCT 
            latitude,
            longitude,
            observation_date,
            sampling_event_identifier,
            group_identifier,
            effort_distance_km,
            effort_area_ha,
            duration_minutes,
            1 AS presence
          FROM ebird_subset
          WHERE floor(taxonomic_order) IN (",paste(sptaxon,collapse=","),"))
          absence as (SELECT DISTINCT
            latitude,
            longitude,
            observation_date,
            sampling_event_identifier,
            group_identifier,
            effort_distance_km,
            effort_area_ha,
            duration_minutes,
            0 AS presence
          FROM ebird_subset
          WHERE all_species_reported='t'
          AND sampling_event_identifier NOT IN (SELECT
            sampling_event_identifier FROM presence))
          SELECT 
            latitude,
            longitude,
            observation_date,
            presence,
            effort_distance_km,
            duration_minutes,
            effort_area_ha 
          FROM presence
         UNION
          SELECT
              latitude,
              longitude,
              observation_date,
              presence,
              effort_distance_km,
              duration_minutes,
              effort_area_ha 
          FROM absence"))
    }
```

Use the `getebird()` function to query the database and return the summarized data frame.

``` {.r}
## get species data
require(redshift)
rs_url="jdbc:postgresql://***redshift.amazonaws.com:5439/mol?tcpKeepAlive=true"

conn <- redshift.connect(rs_url)

spd_all=getebird(
  con=conn,
  sptaxon=sptaxon,
  nulltaxon=NULL,
  region=reg
  )
```

[1] M. Arthur Munson, Kevin Webb, Daniel Sheldon, Daniel Fink, Wesley M. Hochachka, Marshall Iliff, Mirek Riedewald, Daria Sorokina, Brian Sullivan, Christopher Wood, and Steve Kelling. The eBird Reference Dataset, Version 5.0. Cornell Lab of Ornithology and National Audubon Society, Ithaca, NY, January 2013.
