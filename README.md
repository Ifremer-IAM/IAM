
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # IAM : Impact Assessment bio-economic Model for fisheries management <img src="https://gitlab.ifremer.fr/iam/iam/-/raw/dev/inst/fig/IAM_hex.png?inline=false" alt="IAM logo" align="right" height="200px/"/> -->

# IAM : Impact Assessment bio-economic Model for fisheries management <img src="man/figures/logo.png" align="right" height="200px/"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Dev
Version](https://img.shields.io/badge/devel%20version-2.0.0-blue.svg)](https://gitlab.ifremer.fr/iam/iam)
[![License:
CeCILL-2](https://img.shields.io/badge/license-CeCILL--2-blue.svg)](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html)

<!-- [![](https://img.shields.io/badge/devel%20version-2.0.0-blue.svg)](https://github.com/iam/iam) -->
<!-- badges: end -->

IAM (Impact Assessment Model) is a bio-economic model for the simulation
of fisheries dynamics, integrating specific decision support tools
within the framework of the theoretical implementation of management
measures.

It is a discrete time, multi-fleet, multi-trade, multi-species annual
bio-economic model with “age” components for the biological part, and
“commercial category” components for the economic part.

You are free to copy, modify, and distribute IAM with attribution under
the terms of the [CECILL-2
Licence](https://gitlab.ifremer.fr/iam/iam/-/blob/main/LICENCE-CECILL-2.1.txt).
See the LICENSE-CECILL-2.1 file for details.

## Usage

Before using IAM, you need few software :

-   R &gt;= 3.6

-   Rtools (can be installed with `{devtools}` R package)

## Installation

### From Gitlab (recommended)

To install the released version of IAM from
<https://gitlab.ifremer.fr/iam/iam>:

1.  install “remotes” R package by executing in a R console:

``` r
install.packages("remotes")
```

2.  install IAM package by executing in a R console:

``` r
remotes::install_gitlab(
  repo = "iam/iam", 
  host = "https://gitlab.ifremer.fr"
)
```

<!--
### From a package archive file

1. Install IAM dependencies by running in a R console:

```r
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

## Documentation

For further information, multiple vignettes are available along with
function help.

To access a vignette you can run this code
locally<!-- or check them with provided links-->. Please note that to
this day, vignettes are written in french.

``` r
vignette("use_IAM", package = "IAM")
```

A more in depth notice in french is available if you need to check the
model implementation and mathematical formula. The pdf file is available
for [download
here](https://gitlab.ifremer.fr/iam/iam/-/tree/dev/inst/notice/Notice.pdf)

## Credits

IAM is written by Mathieu Merzéréaud, Claire Macher, Michel Bertignac,
Marjolaine Frésard, Olivier Guyader, Christelle Le Grand, Sophie
Gourgue, Florence Briton and Maxime Jaunatre.

## Support

For further information or help, don’t hesitate to fill an [issue on
gitlab](https://gitlab.ifremer.fr/iam/iam/-/issues)
