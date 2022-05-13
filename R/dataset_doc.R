#' Complete IAM dataset
#'
#' Dataset anonymised, with 7 fleets and 3 dynamic species,
#' one being with a SS3 dynamic. Time begin in 2009 and end in 2020
#'
#' @format Two different part of parameter for a IAM.model simulation :\cr
#' \itemize{
#'  \item{IAM_input_2009 - a \code{\link[IAM]{iamInput-class}} object.
#'  Created from a \code{.xlsx file. See Details.}}
#'  \item{IAM_argum_2009 - a \code{\link[IAM]{iamArgs-class}} object.
#'  Created from the input, to modify Recruitment, Gestion and
#'  Economic parameters }
#' }
#'
#' @details
#' This dataset is created with \code{\link[IAM]{IAM.input-methods}} function.
#' This function load a \code{.xlsx} file, require a \code{t_init} and \code{nbStep}.
#' If a SS3 species is included in the \code{.xlsx} file, then additionnal
#' data must be provided.
#'
#' @author Jaunatre Maxime <maxime.jaunatre@yahoo.fr>
#'
#' @examples
#' \dontrun{
#'
#' # Code used to create this object from raw .xlsx file. More information in vignette.
#' load(IAM_example("IAM_SS3_2009.RData"))
#' IAM_input_2009 <- IAM::IAM.input(
#'   fileIN = IAM_example("inputFile.xlsx"),
#'   t_init = 2009, nbStep = 12, folderFleet = IAM_example("fleets"),
#'   Fq_i = list(DAR = iniFq_i), iniFq_i = list(DAR = iniFq_i),
#'   Fq_fmi = list(DAR = iniFq_fmi), iniFq_fmi = list(DAR = iniFq_fmi),
#'   FqLwt_i = list(DAR = iniFqLwt_i),iniFqLwt_i = list(DAR = iniFqLwt_i),
#'   FqLwt_fmi = list(DAR = iniFqLwt_fmi), iniFqLwt_fmi = list(DAR = iniFqLwt_fmi),
#'   FqDwt_i = list(DAR = iniFqDwt_i), iniFqDwt_i = list(DAR = iniFqDwt_i),
#'   FqDwt_fmi = list(DAR = iniFqDwt_fmi), iniFqDwt_fmi = list(DAR = iniFqDwt_fmi),
#'   Nt0s1q = list(DAR = Nt0s1q), Ni0q = list(DAR = Ni0q),
#'   iniNt0q = list(DAR = iniNt0q), matwt = list(DAR = mat_morphage)
#' )
#'
#' IAM_argum_2009 <- IAM.input2args(IAM_input_2009)
#' IAM_argum_2009@arguments$Eco$dr <- 0.04
#' IAM_argum_2009@arguments$Eco$perscCalc <- as.integer(1)
#' }
#'
#' data(IAM_input_2009)
#' data(IAM_argum_2009)
#'
#' @name IAM_dataset
#' @rdname IAM_dataset
#' @source
#' Based on Florence Briton dataset.
"IAM_input_2009"


#' @rdname IAM_dataset
"IAM_argum_2009"
