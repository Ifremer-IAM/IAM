
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # IAM : Impact Assessment bio-economic Model for fisheries management <img src="https://gitlab.ifremer.fr/iam/iam/-/raw/dev/inst/fig/IAM_hex.png?inline=false" alt="IAM logo" align="right" height="200px/"/> -->

# IAM : Impact Assessment bio-economic Model for fisheries management <img src="man/figures/logo.png" align="right" height="200px/"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![IAM status
badge](https://ifremer-iam.r-universe.dev/badges/IAM)](https://ifremer-iam.r-universe.dev)
[![License:
CeCILL-2](https://img.shields.io/badge/license-CeCILL--2-blue.svg)](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html)

<!-- badges: end -->
<!-- old badge -->
<!-- [![Dev Version](https://img.shields.io/badge/dev%20version-2.0.0-blue.svg)](https://github.com/https://gitlab.ifremer.fr/iam/iam)  -->

IAM (Impact Assessment Model for fisheries management) is a bio-economic
model developed as part of a partnership with stakeholders to support
fisheries management. It is a tool for academic and non academic
knowledge integration which models dynamics and interactions between
fish stocks, vessels or fleets, fisheries governance and fish market. It
is dedicated to scenario simulations and optimization, impact assessment
of management strategies (transition to MSY, fisheries Management Plans,
socio-economic consequences of alternative TAC and quotas allocation
options) and exploration of conditions for fisheries viability and
sustainability. It enables stochastic simulations of biological and
socio-economic consequences of scenarios to compare trade-offs of
aletrnative options from a multi-criteria perspective.

It is a discrete time (annual), multi-fleet or multi-vessel,
multi-métier, multi-species bio-economic model with “age” components for
the biological part, and “commercial category” components for the
economic part.

You are free to copy, modify, and distribute `{IAM}` with attribution
under the terms of the [CECILL-2
Licence](https://gitlab.ifremer.fr/iam/iam/-/blob/main/LICENCE-CECILL-2.1.txt).
See the LICENSE-CECILL-2.1 file for details.

Note : the github repository is only a mirror from
[gitlab](https://gitlab.ifremer.fr/iam/iam)

## Installation

`{IAM}` is a `R` package and thus require you to install it with a
version superior to 3.6. Updating `R` (and Rstudio) is recommended
before running a lot of package updates.

### From package repository (recommended)

``` r
install.packages("IAM", repos = "https://ifremer-iam.r-universe.dev")
```

<!-- This is CRAN-like -->

### From a package archive file

1.  Go to the [Project ‘Releases’
    page](https://gitlab.ifremer.fr/iam/iam/-/releases/) and download an
    `{IAM}` binary package (.zip format under Window, tar.gz format
    under Linux) coressponding with your `R` version. These are **not**
    the files listed in assets !

2.  Install IAM dependencies if needed:

``` r
# install.packages("remotes")
remotes::install_deps("--path to zipfile--/IAM_2.0.0.zip")
```

3.  Install IAM from the archive file by entering the following command
    in the `R` prompt:

``` r
install.packages("--path to zipfile--/IAM_2.0.0.zip",
                 repos = NULL, type = "win.binary")
```

<!-- ```{r, installlink, eval = FALSE} -->
<!-- # Also work with direct link  -->
<!-- install.packages("https://gitlab.ifremer.fr/iam/iam/uploads/1252e068a81c5c282bb1599686ef67df/IAM_2.0.0.zip", -->
<!--                  repos = NULL, type = "win.binary") -->
<!-- ``` -->

### From package source code

In order to build packages from sources on Windows, you will need
`Rtools`, which is an independent software from `R`. You can install
`Rtools` from the [CRAN project
website](https://cran.r-project.org/bin/windows/Rtools/) after selecting
the right version. Note that `Rtools 4.2` is still in development
version. You can check your `R` version with `R.version`.

Install IAM from the sources files by entering the following command in
the `R` prompt:

``` r
remotes::install_git(url = "https://gitlab.ifremer.fr/iam/iam")

# You can install development version with the ref option
remotes::install_git(
  url = "https://gitlab.ifremer.fr/iam/iam", 
  ref = "dev", 
  build_vignettes = require(rmarkdown)
)
```

## Documentation

For further information, all documentation is centralized in a
[{pkgdown} website](https://ifremer-iam.github.io/IAM/index.html). IAM
package has multiple vignettes available along with functions help.
Please note that to this day, vignettes are written in french.

To get started, read [Use IAM
vignette](https://ifremer-iam.github.io/IAM/articles/use_IAM.html). You
can also access a vignette by running this code locally :

``` r
vignette("use_IAM", package = "IAM")
```

A more in depth notice in french is available if you need to check the
model implementation and mathematical formula. The pdf file is available
for [download
here](https://gitlab.ifremer.fr/iam/iam/-/raw/main/inst/notice/Notice_IAM.pdf?inline=false)

## Credits

IAM is written by Mathieu Merzéréaud, Claire Macher, Michel Bertignac,
Marjolaine Frésard, Olivier Guyader, Christelle Le Grand, Sophie
Gourguet, Florence Briton and Maxime Jaunatre.

## Support

For further information or help, don’t hesitate to fill an [issue on
gitlab](https://gitlab.ifremer.fr/iam/iam/-/issues)
