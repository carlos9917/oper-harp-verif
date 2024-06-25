# Developments done during Polly Schmederer's ACCORD VS at DMI (June 2024)

Polly Schmederer (Geosphere), Carlos Peralta (DMI) and Fabrizio Baordo (DMI)

**Topics of VS**

Generalizing spatial verifications: improving R scripting; use of reticulate package to interface R with Python; testing 'panelification' tool.

**Data used for spatial verification**

Observations:

DMI's radar precipitation product: Surface Quantitative Precipitation Estimation (SQPE) using both rain guage and radar data

EUMETSAT SEVIRI data (https://api.eumetsat.int/data/browse/collections): High Rate SEVIRI Level 1.5 Image Data - MSG - 0 degree (native), e.g. MSG3-SEVI-MSG15-0100-NA-20240102235743.693000000Z-NA.nat

NWP: 

Grib files output of the DEODE workflow running HARMONIE cy46h1 (total presipitation and FULL POS simulated radiances channels WV_062 & IR_108)

**Content of the repository**

### Installation instructions
Refer to the [installation instructions](INSTALLATION.md) for details of how to install different libraries.

