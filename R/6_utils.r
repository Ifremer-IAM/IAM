

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Methodes utilitaires diverses (mise en forme,...)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# old doc
# dimAbsc = numeric(1) : dim en abscisse

#' Formating method
#'
#' fonction interne d'agregation et de formatage de la donnee a representer
#' selon divers parametres
#'
#' @param object An "iamOutput" or "iamOutputRep" object. # TODO link here
#' @param elmt Name of the operating variable. chr (1).
#' @param spp Name of the considered species (relevant only for some variables). chr (1).
#' @param agg Index(es) describing dimensions on which agregation must be done,
#' referring to 'DimCst' attributes (1=fleet, 2=metier, 3=age/category, 4=time). num (<4)
#' @param headers Optionnal. Dimension that will be developped as headers in
#' the output data.frame. chr.
#' @param subs An optionnal 'subset' argument to operate a final subsetting on output dataframe. chr (1)
#' @param index # TODO argum is ?
#'
#' @return data.frame avec colonnes : value + dimensions
#'
iamFormat <- function(object, elmt = character(), spp = character(), agg = NA,
                      headers = as.character(NA), subs, index = NA) {

#on identifie la dimension d'abscisse
if (headers%in%"Scen") headers <- as.character(NA)
x <- switch(headers, f = 1, m = 2, a_c = 3, t = 4, NULL)
heads <- TRUE

#on commence par s?lectionner la variable

if (!is.null(object@output[[elmt]])) {

  vrbl <- object@output[[elmt]]
  if (!is.na(index)) vrbl <- vrbl[[index]]

} else {

  if (!is.null(object@outputSp[[elmt]][[spp]])) {

    vrbl <- object@outputSp[[elmt]][[spp]]

  } else {

    if (!is.null(object@outputSp[[elmt]][[index]][[spp]])) {

      vrbl <- object@outputSp[[elmt]][[index]][[spp]]

    } else {

      stop("can't find the specified variable in the input object!! Check 'elmt' or 'spp' parameter!!")

    }
  }
}


#il ne faut pas que 'x' soit dans 'agg'
if (!is.null(x)) {if (x%in%agg) stop("'headers' can't be in 'agg'!! Check your parameters!!")}
#il ne faut pas que 'x' soit une dimension nulle de DimCst
if (!is.null(x)) {if (attributes(vrbl)$DimCst[x]==0) stop("'headers'th dimension is not defined in specified variable!!")}
if (is.null(x)) heads <- FALSE


#nouvelles dimensions apr?s agr?gation
if (!is.null(dim(vrbl))) {

  if (!all(is.na(agg))) { #on doit op?rer l'agr?gation

    newAgg <- agg[attributes(vrbl)$DimCst[agg]>0]
    dimApply <- (1:sum(attributes(vrbl)$DimCst>0))[-((1:4)-cumsum(attributes(vrbl)$DimCst==0))[newAgg]]

    if (length(dimApply)==0) stop("can't aggregate over those indexes!!")

    vrblPerm <- apply(vrbl,dimApply,sum,na.rm=TRUE)
    #dimension de la nouvelle variable correspondant ? l'abscisse de repr?sentation
    dimTemp <- attributes(vrbl)$DimCst ; dimTemp[newAgg] <- 0 ;

    if (all(dimTemp%in%0)) stop("no dimensions left after agregation!!")

    nam <- c("f","m","a_c","t")[dimTemp>0]
    #on permutte la matrice et on formate si plus de deux dimensions
    if (!is.null(dim(vrblPerm))) {

      df <- cbind.data.frame(value=as.vector(vrblPerm),expand.grid(dimnames(vrblPerm)))
      names(df) <- c("value",nam)

    } else {

      df <- data.frame(value=as.vector(vrblPerm), absc=names(vrblPerm))
      names(df) <- c("value",nam)

    }

  } else {

    dimTemp <- attributes(vrbl)$DimCst

    nam <- c("f","m","a_c","t")[dimTemp>0]

    df <- cbind.data.frame(value=as.vector(vrbl),expand.grid(dimnames(vrbl)))
    names(df) <- c("value",nam)


  }

} else {    #pas d'agr?gation possible

  dimTemp <- attributes(vrbl)$DimCst
  nam <- c("f","m","a_c","t")[dimTemp>0]

  df <- data.frame(value=as.vector(vrbl), absc=names(vrbl))
  names(df) <- c("value",nam)
}

  #if (length(subset)>0)
  if (!missing(subs)) df <- df[eval(substitute(subs),df),]  #on subsette


  if (heads & (ncol(df)>2) & (nrow(df)>0)) {
    if (!missing(subs)) df <- cbind.data.frame(value=df$value,as.data.frame(as.matrix(df[,2:ncol(df)]))) #on met ? jour les levels apr?s subset
    df <- df[,c("value",headers,names(df)[!names(df)%in%c("value",headers)])] #on met la colonne de headers en deuxi?me position
    mat <- tapply(df[,1],list(do.call('paste',c(df[,3:ncol(df),drop=FALSE],list(sep=":-:-:"))),df[,2]),function(x) x)
    mat1 <- do.call("rbind",lapply(rownames(mat),function(x) strsplit(x,":-:-:")[[1]]))
    df1 <- as.data.frame(mat1) ; names(df1) <- names(df)[-(1:2)]
    df2 <- as.data.frame(mat) ; rownames(df2) <- NULL
    df <- cbind.data.frame(df1,df2)
  }


  rownames(df) <- 1:nrow(df)

  return(df)

}


# Methodes associees

#' Formating method
#'
#' @param object An "iamOutput" or "iamOutputRep" object. # TODO link here
#' @param ... Further arguments :
#' \describe{
#' \item{elmt}{ Name of the operating variable. chr.}
#' \item{spp}{Name of the considered species (relevant only for some variables). chr}
#' \item{agg}{Index(es) describing dimensions on which agregation must be done,
#' referring to 'DimCst' attributes (1=fleet, 2=metier, 3=age/category, 4=time). num.}
#' \item{headers}{Optionnal. Dimension that will be developped as headers in the output data.frame. chr.}
#' \item{subs}{An optionnal 'subset' argument to operate a final subsetting on output dataframe. logic}
#' }
#'
#' @rdname IAM.format-method
#' @export
setGeneric("IAM.format", function(object, ...){ # IAM.format gen ####
	standardGeneric("IAM.format")
	}
)

#' @rdname IAM.format-method
#' @aliases IAM.format,iamOutput-method
setMethod("IAM.format", signature("iamOutput"),function(object, ...){ ## iamOutput meth ####

  iamFormat(object, ...)

})

#' @rdname IAM.format-method
#' @aliases IAM.format,iamOutputRep-method
setMethod("IAM.format", signature("iamOutputRep"),function(object, ...){ ## iamOutputRep meth ####

  lapply(1:object@arguments$Replicates$nbIter,function(x)
                  iamFormat(object, ..., index=x))

})



#::::::::::::::::::
#Examples
#::::::::::::::::::


#out <- IAM.input("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Inputs_SIAD_SEL_2.xls",t_init=2010,nbStep=21)
#
#arg <- IAM.args(out)
#
#mod <- IAM.model(arg,out)
#
#IAM.format(mod,elmt = "C", spp = "Langoustine", agg = 1:2, headers = "a_c", t%in%(2012:2016))





#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Methode composite issue de IAM.format, ajoutant le regroupement des
# scenarii --> input des m?thodes graphiques

#' Formating method for grouping various iamOutput or iamOutputRep objects
#'
#' @param listObj A list containing "iamOutput" OR "iamOutputRep" objects. Mixed lists are not allowed.
#' @param ... Further arguments :
#' \describe{
#' \item{elmt}{ Name of the operating variable. chr.}
#' \item{spp}{Name of the considered species (relevant only for some variables). chr}
#' \item{agg}{Index(es) describing dimensions on which agregation must be done,
#' referring to 'DimCst' attributes (1=fleet, 2=metier, 3=age/category, 4=time). num.}
#' \item{headers}{Optionnal. Dimension that will be developped as headers in the output data.frame. chr.}
#' \item{subs}{An optionnal 'subset' argument to operate a final subsetting on output dataframe. logic}
#' }
#'
#' @rdname IAM.unite
#' @export
setGeneric("IAM.unite", function(listObj,...){  # IAM.unite gen ####
	standardGeneric("IAM.unite")
	}
)

#' @rdname IAM.unite
#' @aliases IAM.unite,list-method
setMethod("IAM.unite", signature(listObj="list"), function(listObj,...){  ## list meth ####

	if (!(all(unlist(lapply(listObj,class))%in%"iamOutput") | all(unlist(lapply(listObj,class))%in%"iamOutputRep")))
    stop("only 'iamOutput' objects OR (not AND) 'iamOutputRep' objects allowed in input list!!")



  df <- do.call("rbind",lapply(listObj, function(x) {
                                  if (class(x)%in%"iamOutput") {
                                    return(cbind.data.frame(IAM.format(x,...),
                                                            Scen=ifelse(x@arguments$Scenario$active==0,
                                                                        "Status Quo",
                                                                        x@arguments$Scenario$ALLscenario[x@arguments$Scenario$SELECTscen])))
                                  } else {
                                    ll <- IAM.format(x,...)
                                    return(cbind.data.frame(do.call("rbind",ll),
                                                            iter=rep(1:length(ll),unlist(lapply(ll,nrow))),
                                                            Scen=ifelse(x@arguments$Scenario$active==0,
                                                                        "Status Quo",
                                                                        x@arguments$Scenario$ALLscenario[x@arguments$Scenario$SELECTscen])))
                                  }

                                  }))


    call <- match.call()
    llhead <- as.list(call)[-1]$headers[1]


    if (!is.null(llhead)) {
      if (llhead%in%"Scen") {

    if ((ncol(df)>2) & (nrow(df)>0)) {
    mat <- tapply(df[,1],list(do.call('paste',c(df[,2:(ncol(df)-1),drop=FALSE],list(sep=":-:-:"))),df[,"Scen"]),function(x) x)
    mat1 <- do.call("rbind",lapply(rownames(mat),function(x) strsplit(x,":-:-:")[[1]]))
    df1 <- as.data.frame(mat1) ; names(df1) <- names(df)[-(c(1,ncol(df)))]
    df2 <- as.data.frame(mat) ; rownames(df2) <- NULL
    df <- cbind.data.frame(df1,df2)
  }

  rownames(df) <- 1:nrow(df)

    }
  }

  return(df)

})



#:::::::::::::::
#Examples
#:::::::::::::::


#arg2 <- IAM.args(arg)
#
#mod2 <- IAM.model(arg2,out)
#
#IAM.unite(list(mod,mod2),elmt = "C", spp = "Langoustine", agg = 1:2, headers = "a_c", t%in%(2012:2025))




