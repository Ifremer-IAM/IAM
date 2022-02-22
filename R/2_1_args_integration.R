#' 'iamArgs' objects creator
#'
#' @param object  An \code{\link[IAM]{iamInput-class}} or \code{\link[IAM]{iamArgs-class}} object.
#' @param ... desc argument
#'
#' @name IAM.input2args
#'
#' @export
setGeneric("IAM.input2args", function(object,  ...){
  standardGeneric("IAM.input2args")
}
)

#' @param desc Object descriptor (default value : \code{as.character(NA)}).
#' If not provided, copied the description slot of object.
#'
#' @rdname IAM.input2args
#' @export
setMethod("IAM.input2args", signature("iamInput"),function(object, desc= NULL){


  ALLVarRep = c(
    "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth", "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
    "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f", "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
    "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f", "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
    "Fr_fmi", "C_efmit", "P", "Pstat"
  )

  if(is.null(desc)){ desc <- object@desc }
  # init the arg object with shiny default
  ## Create argum ####
  # TODO : this is where to modify the default for the GUI now !
  init_recru <- function(name, object){
    list(modSRactive = 1,
         typeMODsr = "Mean",
         parAmodSR = unname(object@input[[name]]$N_it0[1]),
         parBmodSR = 0, parCmodSR = 0, wnNOISEmodSR = 0, noiseTypeSR = 1,
         simuSTOCHactive = 0, typeSIMUstoch = 1
    )
  }

  sp <- as.list(object@specific$Species) ; names(sp) <- sp
  spDyn <- sp[object@specific$Q == 0]
  Recruitement <- lapply(spDyn, function(x) init_recru(name = x, object = object))
  rm(sp)

  Replicates <- list(active =0, nbIter =500,
                     SELECTvar = ALLVarRep[c(1:20, 22:38, 43:44)] )

  Scenario <- list(active = 0, ALLscenario = names(object@scenario), SELECTscen = 1)

  Eco <- list(active = 0, type = 1,
              adj = 1, ue_choice = 1, oths = 0, othsFM = 0, # useless
              perscCalc = 0,
              report = 0, # useless
              dr = 0)

  init_gest <- function(object){
    years <- object@specific$times
    fleets <- object@specific$Fleet
    spp <- c(na.omit(object@specific$Species))
    sppStat <- c(na.omit(object@specific$StaticSpp))
    Espece <- c(na.omit(c(spp,sppStat)))[1]

    tacfbar <- matrix(as.numeric(NA),nrow=2,ncol=length(years),
                      dimnames=list(c("TAC","Fbar"),years))
    mfm <- object@input$Fleet$nbv_f_m ; mfm[!is.na(mfm)] <- 1

    TACOTH <- matrix(as.numeric(NA),nrow=length(sppStat)+length(spp),
                     ncol=length(years),dimnames=list(c(spp,sppStat),years))
    if (nrow(TACOTH)>1) {
      TACOTH <- as.data.frame(t(TACOTH[rownames(TACOTH)!=Espece,]))
      tabO <- lapply(TACOTH, function(x) {
        attr(x, "DimCst") <- as.integer(c(0,0,0,length(x)))
        names(x) <- years
        return(x)
        })
    } else {
      tabO <- NULL
    }
    # TODO : utiliser with pour enlever plein de verbose
    Gestion <- list(active = 0, control = "Nb vessels", target = "TAC",
                    espece = Espece, delay = 2, typeG = 0, upd = 1,
                    sup = 0, inf =0,
                    tac = tacfbar["TAC",],
                    fbar = tacfbar["Fbar",],
                    othSpSup = tabO,
                    effSup = matrix(as.numeric(NA),nrow=length(fleets),
                                    ncol=length(years),dimnames=list(fleets,years)),
                    mfm = mfm,
                    TACbyF = matrix(as.numeric(NA),nrow=length(fleets),
                                    ncol=length(years),dimnames=list(fleets,years)))

    attr(Gestion$tac, "DimCst") <- as.integer(c(0,0,0,length(Gestion$tac)))
    attr(Gestion$fbar, "DimCst") <- as.integer(c(0,0,0,length(Gestion$fbar)))
    attr(Gestion$effSup, "DimCst") <- as.integer(c(nrow(Gestion$effSup),0,0,
                                                      ncol(Gestion$effSup)))
    attr(Gestion$mfm, "DimCst") <- as.integer(c(dim(Gestion$mfm),0,0))
    attr(Gestion$TACbyF, "DimCst") <- as.integer(c(nrow(Gestion$TACbyF),0,0,
                                                      ncol(Gestion$TACbyF)))
    return(Gestion)
  }
  Gestion <- init_gest(object)

  argum <- list(Recruitment = Recruitement, Replicates = Replicates,
                Scenario = Scenario, Gestion = Gestion, Eco = Eco)
  ## Copy specific ####
  specif <- object@specific
  args <- new("iamArgs", desc = desc, arguments = argum, specific = specif)

  return(args)

})
