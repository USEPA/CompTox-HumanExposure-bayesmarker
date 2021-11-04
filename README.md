# R Package "bayesmarker"

This package provides a modeling approach that uses biomonitoring data and chemical metabolism 
information to estimate chemical exposure intake rates while carefully characterizing the uncertainties.
Bayesian methods are provided to infer ranges of exposure for parent chemicals consistent with 
biomarkers identified in urine samples from the U.S population by the National Health and 
Nutrition Examination Survey (NHANES). Metabolites arelinked to potential parent chemicals 
using the NHANES reports and text mining of PubMed abstracts. 

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

Adapted from <a href="https://doi.org/10.1080/17425255.2021.1935867">Breen et al. (2021)</a>
* Getting Started with R Package httk from the R command line
```
install.packages(httk)
```
* RStudio provides a menu ‘Install Packages’ under ‘Tools’ tab
* Load the HTTK data, models, and functions
```
library(httk)
```
* Check what version you are using 
```
packageVersion(httk)
```

## Authors

lead package developer Zach Stanfield
[@stanfield.zach@epa.gov]

Woodrow Setzer
[@setzer.woodrow@epa.gov]

Victoria Hull
[@Hull.Victoria@epa.gov]

Risa Sayre
[@sayre.risa@epa.gov]

Kristin Isaacs
[@isaacs.kristin@epa.gov]

John Wambaugh
[@wambaugh.john@epa.gov]



## License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>