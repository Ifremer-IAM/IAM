% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_args_integration.r
\docType{methods}
\name{IAM.args-methods}
\alias{IAM.args-methods}
\alias{IAM.args}
\alias{IAM.args,iamInput,missing-method}
\alias{IAM.args,iamArgs,missing-method}
\alias{IAM.args,character,character-method}
\title{'iamArgs' objects creator by GUI}
\usage{
IAM.args(object, specific, ...)

\S4method{IAM.args}{iamInput,missing}(object, desc = as.character(NA), ...)

\S4method{IAM.args}{iamArgs,missing}(object, desc = as.character(NA), ...)

\S4method{IAM.args}{character,character}(object, specific, desc = as.character(NA), ...)
}
\arguments{
\item{object}{\code{\link[IAM]{iamInput-class}} or
\code{\link[IAM]{iamArgs-class}} object.}

\item{specific}{this is for the specific file, a deprecated use of
this function.}

\item{...}{desc parameter described below.

 @return \code{\link[IAM]{iamArgs-class}} object after modification by a GUI.}

\item{desc}{Object descriptor (default value : \code{as.character(NA)}).
If not provided, copied the description slot of object.}
}
\description{
'iamArgs' objects creator by GUI

Etape d'initialisation

Etape de modification

Etape d'initialisation from txt file.
}
\examples{
\dontrun{
data(IAM_input_2009)
# example not run because depend on shiny app.
# Initiate argument from input.
IAM_argum_2009 <- IAM.args(IAM_input_2009)
summary(IAM_argum_2009)
# Modify argument from argument.
IAM_argum_2009 <- IAM.args(IAM_argum_2009)
summary(IAM_argum_2009)
}

}
