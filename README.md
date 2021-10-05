---
  output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "inst/images/README-",
  out.width = "100%"
)
library(badger)
```

# IAM : Impact Assessment bio-economic Model for fisheries management
<!--<img src="https://raw.githubusercontent.com/usr/repo/master/inst/images/pic.png" alt="logo" align="right" height="200px/"/>
-->

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) [![Dev Version](https://img.shields.io/badge/devel%20version-2.0.0-blue.svg)](https://gitlab.ifremer.fr/iam/iam) [![License: CeCILL-2](https://img.shields.io/badge/license-CeCILL--2-blue.svg)](https://cran.r-project.org/web/licenses/CeCILL-2)
<!-- badges: end -->

By Mathieu Merzéréaud, Claire Macher, Michel Bertignac, 
Marjolaine Frésard, Olivier Guyader, Christelle Le Grand

IAM (Impact Assessment Model) is a bio-economic model for the simulation of fisheries dynamics, integrating specific decision support tools within the framework of the theoretical implementation of management measures

It is a discrete time, multi-fleet, multi-trade, multi-species annual bio-economic model with "age" components for the biological part, and "commercial category" components for the economic part.


You are free to copy, modify, and distribute IAM with attribution under the terms of the CECILL-2.1 license. See the LICENSE-CECILL-2.1 file for details.

## Usage

Before using IAM, you need few software :
- R >= 3.6

<!--
### Dependencies

This package relies on very few packages listed below, that you can install with the following code.

methods, lattice, stats, grDevices, abind, MASS,
        openxlsx, reshape2

```{r, eval = FALSE}
for (i in c('graphics', 'stats', 'viridisLite') ){
  if(!require(i,character.only = TRUE))
    install.packages(i)
}
```

### Development version

You can install the development version of `{DiveR}` from [github](https://github.com/gowachin/DiveR) with:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github('https://github.com/gowachin/DiveR')
# or 
# install.packages("remotes")
remotes::install_github("gowachin/DiveR")
```

<!--## Usage-->

This is a simple example where we simulate a dive. This show also the desaturation stops due in the table model.

```{r example_dive, dev='png', out.width="100%"}
# Simulation of a dive
library(DiveR)
dive <- dive(depth = 39, time = 22, secu = TRUE, 
             ascent_speed = 10, desat_model = "table")
summary(dive)
```


-->
