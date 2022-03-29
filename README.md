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
- Rtools (can be installed with `{devtools}` R package)

<!--
## Installation

### From Gitlab (recommended)

To install the released version of IAM from https://gitlab.ifremer.fr/iam/iam:

1. install "remotes" R package by executing in a R console:

`install.packages("remotes")` 

2. install IAM package by executing in a R console:

`remotes::install_gitlab(repo="iam/iam",host="https://gitlab.ifremer.fr")`

### From a package archive file

1. Install IAM dependencies by running in a R console:
```{r, eval = FALSE}
for (i in c('utils', 'stats', 'methods', 'grDevices', 'abind',
         'reshape2', 'openxlsx', 'lattice', 'tcltk', 'tcltk2',
         'Rcpp') ){
  if(!require(i,character.only = TRUE))
    install.packages(i)
}
```

2. Go to the [Project 'Releases' page](https://gitlab.ifremer.fr/iam/iam/releases) and download an IAM binary package (.zip format under Window, tar.gz format under Linux)
3. Install IAM from the archive file by either entering the following command in the R prompt: 
`install.packages(path_to_file, repos = NULL, type="source")`
or use e.g. the RStudio graphical interface
-->


<!-- Not working because private repo

### Development version

You can install the development version of `{IAM}` from [github](https://gitlab.ifremer.fr/iam/iam) with:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_gitlab('gitlab.ifremer.fr/iam/iam')
# or 
# install.packages("remotes")
remotes::install_gitlab("gitlab.ifremer.fr/iam/iam")
```

-->


## See the Wiki section for more support

If you want to contribute code, please email Mathieu.Merzereaud[at]ifremer.fr
