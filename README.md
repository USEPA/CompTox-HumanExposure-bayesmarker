# CompTox-HumanExposure-bayesmarker

## R Package "bayesmarker"

This package provides a modeling approach that uses biomonitoring data and chemical metabolism 
information to estimate chemical exposure intake rates while carefully characterizing the uncertainties. 
Bayesian methods are provided to infer ranges of exposure for parent chemicals consistent with 
biomarkers identified in urine samples from the U.S population by the National Health and 
Nutrition Examination Survey (NHANES). Metabolites are linked to potential parent chemicals 
using the NHANES reports and text mining of PubMed abstracts. The inferred parent chemical
exposures have been incorporated into the US EPA's CompTox Chemicals Dashboard for chemicals
with available data. The exposure predictions can be found for an individual chemical by clicking 
Exposure -> Monitoring Data (for example: comptox.epa.gov/dashboard/chemical/monitoring-data/DTXSID7020182). 


## Background

Biomonitoring data provides insight into the chemicals commonly present in humans, but it is 
important to know which environmental chemicals contribute to the observed metabolites for 
estimating exposure and risk. These inferred exposures, when coupled with predicted toxic doses and high throughput TK, 
allow the identification of public health priority chemicals via risk-based bioactivity:exposure ratios (BER).


## Getting Started

### Dependencies

* Users will need the freely available R statistical computing language: <https://www.r-project.org/>
* Users will likely want a development environment like RStudio: <https://www.rstudio.com/products/rstudio/download/>

### Installing

* Getting Started with R Package bayesmarker from the R command line
``` 
library(devtools)
install_github("USEPA/CompTox-HumanExposure-bayesmarker")
```
* RStudio provides a menu ‘Install Packages’ under ‘Tools’ tab
* Load the bayesmarker data and functions
``` 
library(bayesmarker)
```
* Check what version you are using
```
packageVersion(bayesmarker)
```

## Authors

Lead package developer Zach Stanfield 
[@stanfield.zachary@epa.gov]

Woodrow Setzer 
[@woodrow.setzer@epa.gov]

Victoria Hull 
[@hull.victoria@epa.gov]

Risa Sayre 
[@sayre.risa@epa.gov]

Kristin Isaacs 
[@isaacs.kristin@epa.gov]

John Wambaugh 
[@wambaugh.john@epa.gov]



## License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>

