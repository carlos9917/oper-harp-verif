# Repository for the developments of the ACCORD VS of Carlos Peralta in GeoSphere (July 2025)

**Topics of VS**

Spatial verifications for weather extremes

  - Review of spatial verification methods, with a focus on weather extremes.

  - Implementation of new scores in `harpSpatial`
  
  - Example scripts to run the verification with the new scores.


**Data used for spatial verification**

Observations:
  
  - DMI's radar precipitation product: Surface Quantitative Precipitation Estimation (SQPE) using both rain gauge and radar data.
  
  - INCA data

NWP:

  - Grib files output of the DINI model
  - Grib files output of the CLAEF model

### Installation instructions

Development was done on ATOS (using group "accord" for common access)

Refer to the [installation instructions](INSTALLATION.md) for details of how to install the libraries.

### Overview of the repository

* **docs** 
   
- `general_rev`: General review on spatial methods.
- `review_sass`: General review and validation of the structure of local extrems (SLX) score
- `review_dey`: General review and validation of the agreement scales score
- `review_keil`: General review of the DES score 
- `case_studies`: Some case studies
- `report`: Report documenting the work 

  
* **scripts** 
  - ``example_call_slx_geosphere.R``: Example to call `spatial_verif` with the `SLX` score GeoSphere data (INCA and CLAEF)
  - ``example_call_slx_dmi.R``: Example to call `spatial_verif` using DMI data (SQPE and DINI)
  - ``example_call_.R``: Example to call `spatial_verif` with `agreement_scales` score GeoSphere data (INCA and CLAEF)
  - ``example_call_.R``: Example to call `spatial_verif` with `agreement_scales` score GeoSphere data (INCA and CLAEF)
    
  
## Libraries installation
Currently installing these.
Maybe makes more sense using my fork of harpSpatial dev version?

```
Start R session

install.packages("remotes")
library(remotes)
install_github("harphub/harp")
install_github("harphub/harpIO", "develop")
install_github("harphub/harpCore", "develop")
install_github("harphub/harpSpatial", "develop")
```

