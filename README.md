# R Package "bayesmarker"

This package provides a modeling approach that uses biomonitoring data and chemical metabolism 
<<<<<<< HEAD
information to estimate chemical exposure intake rates while carefully characterizing the uncertainties. 
Bayesian methods are provided to infer ranges of exposure for parent chemicals consistent with 
biomarkers identified in urine samples from the U.S population by the National Health and 
Nutrition Examination Survey (NHANES). Metabolites arelinked to potential parent chemicals 
using the NHANES reports and text mining of PubMed abstracts.
=======
information to estimate chemical exposure intake rates while carefully characterizing the uncertainties.
Bayesian methods are provided to infer ranges of exposure for parent chemicals consistent with 
biomarkers identified in urine samples from the U.S population by the National Health and 
Nutrition Examination Survey (NHANES). Metabolites arelinked to potential parent chemicals 
using the NHANES reports and text mining of PubMed abstracts. 
>>>>>>> 5d1cb5e57807c898dab17daf0c05bf1bacb4eadf

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
<<<<<<< HEAD
``` 
=======
```
>>>>>>> 5d1cb5e57807c898dab17daf0c05bf1bacb4eadf
library(devtools)
install_github("HumanExposure/bayesmarker")
```
* RStudio provides a menu ‘Install Packages’ under ‘Tools’ tab
* Load the bayesmarker data and functions
<<<<<<< HEAD
``` 
library(bayesmarker)
```
* Check what version you are using
``` 
=======
```
library(bayesmarker)
```
* Check what version you are using 
```
>>>>>>> 5d1cb5e57807c898dab17daf0c05bf1bacb4eadf
packageVersion(bayesmarker)
```

## Authors

<<<<<<< HEAD
lead package developer Zach Stanfield 
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
=======
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
>>>>>>> 5d1cb5e57807c898dab17daf0c05bf1bacb4eadf
[@wambaugh.john@epa.gov]



## License

<<<<<<< HEAD
License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>
=======
License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>
>>>>>>> 5d1cb5e57807c898dab17daf0c05bf1bacb4eadf
