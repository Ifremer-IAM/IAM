#' Complete IAM dataset
#'
#' Dataset anonymised, with 7 fleets and 3 dynamic species,
#' one being with a SS3 dynamic. Time begin in 1984 and end in 1995
#'
#' @format Two different part of parameter for a IAM.model simulation :\cr
#' \itemize{
#'  \item{IAM_input_1984 - a \code{\link[IAM]{iamInput-class}} object.
#'  Created from a \code{.xlsx file. See Details.}}
#'  \item{IAM_argum_1984 - a \code{\link[IAM]{iamArgs-class}} object.
#'  Created from the input, to modify Recruitment, Gestion and
#'  Economic parameters }
#' }
#'
#' @details
#' This dataset is created with \code{\link[IAM]{IAM.input}} function.
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
#' load(IAM_example("IAM_SS3_1984.RData"))
#' IAM_input_1984 <- IAM::IAM.input(
#'   fileIN = IAM_example("inputFile.xlsx"),
#'   t_init = 1984, nbStep = 12, folderFleet = IAM_example("fleets"),
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
#' IAM_argum_1984 <- IAM.input2args(IAM_input_1984)
#' IAM_argum_1984@arguments$Eco$dr <- 0.04
#' IAM_argum_1984@arguments$Eco$perscCalc <- as.integer(1)
#' }
#'
#' data(IAM_input_1984)
#' data(IAM_argum_1984)
#'
#' @name IAM_dataset
#' @rdname IAM_dataset
#' @source
#' Based on Florence Briton dataset.
"IAM_input_1984"


#' @rdname IAM_dataset
"IAM_argum_1984"
