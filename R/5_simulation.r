

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# M?thode de simulation, avec en entr?e deux objets de classe 'iamArgs' et 'iamInput'
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#' Running an IAM simulation.
#'
#' @param objArgs output from \code{\link[IAM]{IAM.args-methods}}
#' @param objInput output from \code{\link[IAM]{IAM.input-methods}}
#' @param ... there is so much more arg.... # TODO
#'
#' @importFrom methods new
#'
#' @useDynLib IAM, .registration = TRUE
#' @docType methods
#' @name IAM.model-methods
#' @aliases IAM.model
#' @export
setGeneric("IAM.model", function(objArgs, objInput, ...){ #  Generic model####
	standardGeneric("IAM.model")
	}
)

#' @rdname IAM.model-methods
#'
#' @param verbose Show messages to follow the process. Mainly used for debug process.
setMethod("IAM.model", signature("iamArgs","iamInput"),function(objArgs, objInput, desc=as.character(NA), mOTH=0, updateE=0,
                  TACbyF=NULL, TACtot=NULL, Ftarg=NULL, W_Ftarg=NULL, MeanRec_Ftarg=NULL,#sont g?n?r?s en interne Ztemp, SPPstatOPT, SPPspictOPT et SPPdynOPT, qui sont ins?r?s dans l'export tacCTRL
                  TACbyFoptimCTRL=list(maxIter = as.integer(7), diffZmax = 0.0001, lambda = 0.9, t_stop = 0),
                  recList=list(), recParamList=list(), #new 24/04/2018  31/05/2018
                  ParamSPMList = list(), #added 16/09/19 aleatoire pour Global Surplus Production Model
                  parBehav=list(active=as.integer(0),type=as.integer(3),FMT=NULL,MU=NULL,MUpos=as.integer(0),ALPHA=NULL),
                  parOptQuot=list(active=as.integer(0),pxQuIni=NULL, pxQuMin=NULL, pxQuMax=NULL, lambda=0.1, sdmax=0, ftol=0.0000001, itmax=500),
                  tacControl=list(tolVarTACinf=NA,tolVarTACsup=NA,corVarTACval=NA,corVarTACnby=2,Blim=NA,Bmax=NA,BlimTrigger=as.integer(0),typeMng=NA),
                  stochPrice=list(), #liste d'?l?ments nomm?s par esp?ce consid?r?e, chaque ?l?ment ?tant une liste selon le sch?ma :
                                #list(type=NA (ou 1 ou 2,...), distr=c("norm",NA,NA,NA) (ou "exp" ou...), parA=c(0,NA,NA,NA), parB=c(1,NA,NA,NA), parC=c(NA,NA,NA,NA))
                  parOQD=list(activeQR=as.integer(0),listQR=NULL,listQR_f=NULL),      #10/07/17   activeQR=0 => d?sactiv?, sinon, commence ? l'instant sp?cifi?
                  verbose = FALSE, force_t = NULL, ...){


specific <- objInput@specific
if(is.null(force_t)){force_t = specific$NbSteps}

verbose <- verbose || app_dev()
#Ajout 20/09/2018
nT <- specific$NbSteps
nF <- length(specific$Fleet)
ni <- lapply(specific$Ages,length)

#on compl?te MeanRec_Ftarg avec une potentielle esp?ce Spict
NamSpict = intersect(intersect(names(ni[ni==1]),names(Ftarg)),names(W_Ftarg))

if (length(NamSpict >0) ){
  MeanRec_Ftarg[NamSpict] = NA
}

#validation doublet arguments Ftarg/W_Ftarg
intrs <- intersect(names(Ftarg),names(W_Ftarg))
if (!is.null(Ftarg) & !is.null(W_Ftarg) & length(intrs)>0) {
  Ftarg <- Ftarg[intrs] ; W_Ftarg <- W_Ftarg[intrs] ; MeanRec_Ftarg <- MeanRec_Ftarg[intrs]
  Ftarg <- lapply(Ftarg,function(x) rep(as.double(x),length=nT))
  empty <- lapply(W_Ftarg,function(x) if ((nrow(x)!=(nF+1)) | (ncol(x)!=nT)) stop("Check your 'W_Ftarg' input !!!"))

  W_Ftarg = rapply( W_Ftarg, f=function(x) ifelse(is.na(x),0,x), how="replace" )

  #names(Ftarg) <- names(MeanRec_Ftarg) <- intrs
  TACbyF <- lapply(Ftarg,function(x) matrix(as.double(NA),nrow=nF,ncol=nT,dimnames=list(specific$Fleet,specific$times)))
  TACtot <- lapply(Ftarg,function(x) {tmp <- rep(as.double(NA),length=nT) ; names(tmp) <- specific$times ; return(tmp)})
} else {
  Ftarg <- W_Ftarg <- MeanRec_Ftarg <- NULL
}

if (!is.null(MeanRec_Ftarg)){ # ajout Florence
  MeanRec_Ftarg <- lapply(MeanRec_Ftarg,function(x){
    if (length(x)==1) {
      return(as.integer(x))    #moyenne mobile sur d?lai=n
    } else {
      if ((is.null(nrow(x)) & (length(x)==nT)) | (is.matrix(x) && (nrow(x)==4) && (ncol(x)==nT)) | (is.matrix(x) && (nrow(x)==2) && (ncol(x)==nT))){
        return(as.double(x))    #for?age recrutement XSA(vec:length=nT) ou SS3(mat:dim=4*nT) ou sex-based(mat:dim=2*nT)
      } else {
        stop("Check your 'MeanRec_Ftarg' input !!!")
      }
    }
  })
}

SPPdyn <- lengths(specific$Ages)

#Ajout 27/03/2018 ----------------
#TACbyF <- TACbyF[names(TACbyF)%in%names(TACtot)]
if ((length(TACbyF)==0) | (length(TACtot)==0)) {
 # message("Pas d'ajustement TAC opere car 'TACbyF' ou 'TACtot' est manquant!!")
 TACbyF <- TACtot <- NULL
 SPPstatOPT <- SPPspictOPT <- SPPdynOPT <- integer(0)
} else {
 #if (length(TACtot)==0) TACtot <- list()
 TACtot <- TACtot[names(TACbyF)] ; names(TACtot) <- names(TACbyF)     #TACbyF et TACtot listes de structure similaire # TODO : why second statment ?
 SPPstatOPT <- match(names(TACbyF),specific$StaticSpp) ; SPPstatOPT <- SPPstatOPT[!is.na(SPPstatOPT)] ; if (length(SPPstatOPT)==0) SPPstatOPT <- integer(0)
 SPPspictOPT <- match(names(TACbyF)[names(TACbyF)%in%names(SPPdyn[SPPdyn==1])],specific$Species) ; SPPspictOPT <- SPPspictOPT[!is.na(SPPspictOPT)]
 if (length(SPPspictOPT)==0) SPPspictOPT <- integer(0)
 SPPdynOPT <- match(names(TACbyF)[names(TACbyF)%in%names(SPPdyn[SPPdyn>1])],specific$Species) ; SPPdynOPT <- SPPdynOPT[!is.na(SPPdynOPT)]
 if (length(SPPdynOPT)==0) SPPdynOPT <- integer(0)
 for(i in 1:length(TACbyF)) attributes(TACbyF[[i]])$DimCst <- as.integer(c(nF,0,0,nT))
}

Ztemp <- lapply(specific$Species,function(x, SPPdyn) {
  if (specific$Q[x]==1){
    rep(as.numeric(0),16*SPPdyn[x])
  }  else if (specific$S[x]==1){
    rep(as.numeric(0),2*SPPdyn[x])
  }  else {
    rep(as.numeric(0),SPPdyn[x])
  }
}, SPPdyn)
names(Ztemp) <- specific$Species
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Add SR info for SS3 species 19/11/2014
objArgsIni <- objArgs
if (any(specific$Q == 1)) {
   listType <- list(modSRactive=as.integer(0),typeMODsr="Mean",parAmodSR=as.double(NA),parBmodSR=as.double(0),parCmodSR=as.double(0),
                    wnNOISEmodSR=as.double(0),noiseTypeSR=as.integer(1),simuSTOCHactive=as.integer(0),typeSIMUstoch=as.integer(1))
   namQ <- names(specific$Q)[specific$Q == 1]
   lll <- lapply(namQ,function(z) return(listType))
   names(lll) <- namQ
   newL <- c(objArgs@arguments$Recruitment,lll)
   objArgs@arguments$Recruitment <- newL[specific$Species]
} # Eof SS3 SR


#on v?rifie le formatage des ?l?ments de recList
if (length(recList)>0) {
  devRecL <- recList[specific$Species]   #devRecL -> liste de taille nbE et avec NULL si pas d'info dans recList
  for (elem in 1:length(devRecL)) {
    if (!is.null(devRecL[[elem]])) {
      if(specific$Q[elem]==1){
        if (!is.matrix(devRecL[[elem]])) {
          devRecL[[elem]] <- matrix(rep(devRecL[[elem]],length=nT*4),ncol=nT)    #vecteur ? r?pliquer dans une matrice au format convenable
        } else {
          matTMP <- matrix(as.numeric(0),ncol=nT,nrow=4)           #matrice ? int?grer dans une matrice au format convenable
          matTMP[1:min(4,nrow(devRecL[[elem]])),1:min(nT,ncol(devRecL[[elem]]))] <- as.numeric(devRecL[[elem]][1:min(4,nrow(devRecL[[elem]])),1:min(nT,ncol(devRecL[[elem]]))])
          devRecL[[elem]] <- matTMP
        }
      } else if (specific$S[elem]==1){
        if (!is.matrix(devRecL[[elem]])) {
          devRecL[[elem]] <- matrix(rep(devRecL[[elem]],length=nT*2),ncol=nT)    #vecteur ? r?pliquer dans une matrice au format convenable
        } else {
          matTMP <- matrix(as.numeric(0),ncol=nT,nrow=2)           #matrice ? int?grer dans une matrice au format convenable
          matTMP[1:min(2,nrow(devRecL[[elem]])),1:min(nT,ncol(devRecL[[elem]]))] <- as.numeric(devRecL[[elem]][1:min(2,nrow(devRecL[[elem]])),1:min(nT,ncol(devRecL[[elem]]))])
          devRecL[[elem]] <- matTMP
        }
      } else {#XSA
        devRecL[[elem]] <- rep(as.numeric(as.character(c(devRecL[[elem]],rep(rev(devRecL[[elem]])[1],100)))),length=nT)
      }
    }
  }
  recList <- devRecL #on remplace l'argument initial par l'argument format?
}


#on v?rifie le formatage des ?l?ments de recParamList
if (length(recParamList)>0) {
  devRecParamL <- recParamList[specific$Species]   #devRecParamL -> liste de taille nbE et avec NULL si pas d'info dans recParamList
  for (elem in 1:length(devRecParamL)) {
    if (!is.null(devRecParamL[[elem]])) {
      if (specific$Q[elem]==0 & specific$S[elem]==0 & length(specific$Ages[[elem]])>1) {#XSA
        del <- as.integer(as.character(specific$Ages[[elem]][1]))
        devRecParamL[[elem]] <- list(param=devRecParamL[[elem]][1:nT,!(colnames(devRecParamL[[elem]])=="type")],
                                     delay=del,
                                     type = as.integer(devRecParamL[[elem]][1:nT,"type"]))
        devRecParamL[[elem]]$param[] <- as.numeric(as.character(unlist(devRecParamL[[elem]]$param[]))) ; devRecParamL[[elem]]$param[1:max(1,del),] <- as.numeric(NA)
      } else if (specific$Q[elem]==0 & specific$S[elem]==1 & length(specific$Ages[[elem]])>1){ #sex-based
        Nze <- c(objInput@input[[elem]]$N_i0t_G1[1],objInput@input[[elem]]$N_i0t_G2[1])
        del <- as.integer(as.character(specific$Ages[[elem]][1]))
        devRecParamL[[elem]] <- list(param=devRecParamL[[elem]][1:nT,!(colnames(devRecParamL[[elem]])=="type")],
                                     delay=del,ventil=as.numeric(as.character(Nze/sum(Nze,na.rm=TRUE))),
                                     type = as.integer(devRecParamL[[elem]][1:nT,"type"]))
        devRecParamL[[elem]]$param[] <- as.numeric(as.character(unlist(devRecParamL[[elem]]$param[]))) ; devRecParamL[[elem]]$param[1:max(1,del),] <- as.numeric(NA)

      } else if (specific$Q[elem]==1) { #SS3
        Nze <- c(objInput@input[[elem]]$Ni0_S1M1,objInput@input[[elem]]$Ni0_S2M2,objInput@input[[elem]]$Ni0_S3M3,objInput@input[[elem]]$Ni0_S4M4)
        del <- as.integer(as.character(specific$Ages[[elem]][1]))
        devRecParamL[[elem]] <- list(param=devRecParamL[[elem]][1:nT,!(colnames(devRecParamL[[elem]])=="type")],
                                     delay=del,ventil=as.numeric(as.character(Nze/sum(Nze,na.rm=TRUE))),
                                     type = as.integer(devRecParamL[[elem]][1:nT,"type"]))
        devRecParamL[[elem]]$param[] <- as.numeric(as.character(unlist(devRecParamL[[elem]]$param[]))) ; devRecParamL[[elem]]$param[1:max(1,del),] <- as.numeric(NA)
      }
    }else {
      devRecParamL[elem] <- list(NULL)
    }

  }

  names(devRecParamL) <- specific$Species
  recParamList <- devRecParamL
}

#on v?rifie le formatage des ?l?ments de ParamSPMList
if (length(ParamSPMList)>0) {
  devParamSPMList <- ParamSPMList[specific$Species]   #devRecParamL -> liste de taille nbE et avec NULL si pas d'info dans recParamList
  for (elem in 1:length(devParamSPMList)) {
    if (!is.null(devParamSPMList[[elem]])) {
      devParamSPMList[[elem]][] = as.numeric(devParamSPMList[[elem]][1:nT,])
    }else {
      devParamSPMList[elem] <- list(NULL)
    }

  }

  names(devParamSPMList) <- specific$Species
  ParamSPMList <- devParamSPMList
}

#on ?tend les listes 'parOQD' ? l'ensemble des esp?ces mod?lis?es                                          #10/07/17
allSpp <- c(specific$Species,specific$StaticSpp)                                           #10/07/17
listQR_TMP <- lapply(parOQD$listQR,function(z) rep(z,length=nT)) ; names(listQR_TMP) <- names(parOQD$listQR) ; parOQD$listQR <- listQR_TMP      #10/07/17
listQR_f_TMP <- lapply(parOQD$listQR_f,function(z) if ((nrow(z)!=nF) & (ncol(z)!=nT)) return(NULL) else return(z)) ; names(listQR_f_TMP) <- names(parOQD$listQR_f) ; parOQD$listQR_f <- listQR_f_TMP    #10/07/17
newParOQD <- list(activeQR=parOQD$activeQR,listQR=parOQD$listQR[allSpp],listQR_f=parOQD$listQR_f[allSpp])  #10/07/17
names(newParOQD$listQR) <- names(newParOQD$listQR_f) <- allSpp                                             #10/07/17

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

newStochPrice <- lapply(stochPrice,function(y) return(list(type=as.integer(y$type[1]),
                                 distr=as.character(rep(c(y$distr,rep(NA,4)),length=4)),
                                 parA=as.double(rep(c(y$parA,rep(NA,4)),length=4)),
                                 parB=as.double(rep(c(y$parB,rep(NA,4)),length=4)),
                                 parC=as.double(rep(c(y$parC,rep(NA,4)),length=4)))))



if (objArgs@arguments$Scenario$active==1) {
  scenar <- objArgs@arguments$Scenario$ALLscenario[objArgs@arguments$Scenario$SELECTscen]
 } else {
  scenar <- ""
 }

if (is.null(tacControl$tolVarTACinf)) tacControl$tolVarTACinf <- NA
if (is.null(tacControl$tolVarTACsup)) tacControl$tolVarTACsup <- NA
if (is.null(tacControl$corVarTACval)) tacControl$corVarTACval <- NA
if (is.null(tacControl$corVarTACnby)) tacControl$corVarTACnby <- NA
if (is.null(tacControl$Blim)) tacControl$Blim <- NA
if (is.null(tacControl$Bmax)) tacControl$Bmax <- NA
if (is.null(tacControl$BlimTrigger)) tacControl$BlimTrigger <- as.integer(0)   #application de l'ajustement restrictif du Fmsy en fonction de la SSB ???

def_Px <- lapply(TACtot,function(x) {tmp <- rep(as.double(0.0),length=nT) ; names(tmp) <- specific$times ; return(tmp)})
if (is.null(parOptQuot$active)) parOptQuot$active <- as.integer(0)
if (is.null(parOptQuot$pxQuIni)) parOptQuot$pxQuIni <- def_Px
if (is.null(parOptQuot$pxQuMin)) parOptQuot$pxQuMin <- def_Px
if (is.null(parOptQuot$pxQuMax)) parOptQuot$pxQuMax <- NULL
if (is.null(parOptQuot$lambda)) parOptQuot$lambda <- 0.1
if (is.null(parOptQuot$sdmax)) parOptQuot$sdmax <- 0
if (is.null(parOptQuot$ftol)) parOptQuot$ftol <- 0.0000001
if (is.null(parOptQuot$itmax)) parOptQuot$itmax <- 500
def_Holdings <- lapply(parOptQuot$pxQuIni,function(x) {tmp <- rep(as.double(0.0),length=nT) ; names(tmp) <- specific$times ; return(tmp)})
if (is.null(parOptQuot$holdings)) parOptQuot$holdings <- def_Holdings



if (is.null(parBehav$active)) parBehav$active <- as.integer(0)
if (is.null(parBehav$type)) parBehav$type <- as.integer(3)
if (is.null(parBehav$MUpos)) parBehav$MUpos <- as.integer(0)

Rectyp <- unlist(lapply(objArgs@arguments$Recruitment,function(x) x$simuSTOCHactive * x$typeSIMUstoch))

mOth <- rep(mOTH,length=length(specific$Species)) # ; mOth[match(objArgs@arguments$Gestion$espece,specific$Species)] <- mOTH

TRGT <- match(objArgs@arguments$Gestion$target,c("TAC","Fbar","TAC->Fbar"))
if (objArgs@arguments$Gestion$target%in%"biomasse") TRGT <- 999

# browser()

if(verbose) cat('\n ---- C++ node begin ----\n')
out <-  .Call("IAM", objInput@input, objInput@specific, objInput@stochastic, objInput@scenario[[scenar]],
                    RecType1=as.integer(Rectyp==1), RecType2=as.integer(Rectyp==2), RecType3=as.integer(Rectyp==3),
                    Scenarii = as.integer(objArgs@arguments$Scenario$active), Bootstrp = as.integer(objArgs@arguments$Replicates$active),
                    nbBoot = as.integer(objArgs@arguments$Replicates$nbIter),
                    GestInd = as.integer(objArgs@arguments$Gestion$active),
                    mOth = as.double(mOth),
                    bounds = as.double(c(objArgs@arguments$Gestion$inf,objArgs@arguments$Gestion$sup)),
                    TAC = TACtot, TACglob = as.double(objArgs@arguments$Gestion$tac), FBAR = as.double(objArgs@arguments$Gestion$fbar),
                    effSup = as.double(objArgs@arguments$Gestion$effSup),
                    GestParam = as.integer(c(eTemp = match(objArgs@arguments$Gestion$espece,c(specific$Species,specific$StaticSpp))-1,
                                 var = match(objArgs@arguments$Gestion$control,c("Nb trips","Nb vessels")),
                                 trgt = TRGT,
                                 delay = objArgs@arguments$Gestion$delay,
                                 upd = objArgs@arguments$Gestion$upd, typeG = objArgs@arguments$Gestion$typeG)), # eointeger
                    persCalc = as.integer(c(perscCalc = objArgs@arguments$Eco$perscCalc)),
                    dr = as.double(objArgs@arguments$Eco$dr),
                    SRind = as.integer(unlist(lapply(objArgs@arguments$Recruitment,function(x) x$modSRactive))),
                    listSR = lapply(objArgs@arguments$Recruitment,function(x) as.double(c(rep(x$parAmodSR,length=nT),rep(x$parBmodSR,length=nT),
                            rep(x$parCmodSR,length=nT),rep(x$wnNOISEmodSR,length=nT),rep(x$noiseTypeSR,length=nT)))),#modif MM 27/08/2013 : permet de definir un jeu de parametres SR par annee en vectorisant (a la main) chaque composante
                    TypeSR = lapply(objArgs@arguments$Recruitment,function(x)
                                as.integer(match(x$typeMODsr,c("Mean","Hockey-Stick","Beverton-Holt","Ricker","Shepherd","Quadratic-HS","Smooth-HS")))),
                    mFM = as.double(objArgs@arguments$Gestion$mfm),
                    TACbyF = TACbyF,
                    Ftarg = Ftarg, W_Ftarg = W_Ftarg, MeanRec_Ftarg = MeanRec_Ftarg,
                    parBHV = parBehav,
                    parQEX = list(active=as.integer(parOptQuot$active),pxQuIni=parOptQuot$pxQuIni, pxQuMin=parOptQuot$pxQuMin,
                          pxQuMax=parOptQuot$pxQuMax, lambda=as.double(parOptQuot$lambda),sdmax=as.double(parOptQuot$sdmax), ftol=as.double(parOptQuot$ftol), itmax = as.integer(parOptQuot$itmax),
                         holdings = parOptQuot$holdings),                           #fonctionne en conjugaison avec TACbyF
                    tacCTRL = list(tolVarTACinf=as.double(tacControl$tolVarTACinf),tolVarTACsup=as.double(tacControl$tolVarTACsup),
                          corVarTACval=as.double(tacControl$corVarTACval),corVarTACnby=as.integer(tacControl$corVarTACnby),
                          Blim=as.double(tacControl$Blim),Bmax=as.double(tacControl$Bmax),BlimTrigger=as.integer(tacControl$BlimTrigger),typeMng=as.integer(tacControl$typeMng),
                          maxIter=as.integer(TACbyFoptimCTRL$maxIter),diffZmax=as.double(TACbyFoptimCTRL$diffZmax),lambda=as.double(TACbyFoptimCTRL$lambda),t_stop=as.integer(TACbyFoptimCTRL$t_stop),
                          Ztemp=Ztemp, SPPstatOPT=SPPstatOPT, SPPspictOPT=SPPspictOPT, SPPdynOPT=SPPdynOPT, recList=recList, recParamList=recParamList, ParamSPMList=ParamSPMList), #for?age recrutements ins?r?
                    stochPrice = newStochPrice,       #liste d'?l?ments esp?ce (pas forc?ment toutes pr?sentes, liste vide aussi possible)
                                         #de format d?crit par la ligne de code de construction de 'newStochPrice'
                    updateE = as.integer(updateE),
                    parOQD = newParOQD,                                                  #10/07/17
                    bootVar = as.character(objArgs@arguments$Replicates$SELECTvar),
                    verbose = as.integer(verbose),
                    force_t = as.integer(force_t)
              )


if (objArgs@arguments$Replicates$active==1) {     #objet de classe 'iamOutputRep'

  if (is.na(desc)) desc <- "My iamOutputRep object"

  return(new("iamOutputRep", desc=desc, arguments=objArgsIni@arguments, specific=objInput@specific,
              outputSp = list(
                    F = out$Ffmi,
                    Fr = out$Fr_fmi,
                    Fothi = out$Foth,
                    Fbar = out$Fbar,
                    Z = out$Zeit,
                    N = out$N,
                    B = out$B,
                    SSB = out$SSB,
                    C = out$C_efmit,
                    Ctot = out$Ctot,
                    Y = out$Yfmi,
                    Ytot = out$Ytot,
                    D = out$D_efmit,
                    Li = out$L_efmit,
                    GVL_f_m_e = out$GVL_fme,
                    statY = out$Ystat,
                    statL = out$Lstat,
                    statD = out$Dstat,
                    statGVL_f_m = out$StatGVL_fme),
              output = list(nbv_f = lapply(out$Eff,function(x) x$nbv_f),
                  nbds_f = lapply(out$Eff,function(x) x$nbds_f),
                  nbv_f_m = lapply(out$Eff,function(x) x$nbv_f_m),
                  nbds_f_m = lapply(out$Eff,function(x) x$nbds_f_m),
                  GVLtot_f_m = out$GVLtot_fm,
                  GVLtot_f = out$GVLtot_f,
                  GVLav_f = out$GVLav_f,
                  vcst_f = out$vcst_f,
                  vcst_f_m = out$vcst_fm,
                  rtbs_f = out$rtbs_f,
                  rtbsAct_f = out$rtbsAct_f,
                  cs_f = out$cs_f,
                  csAct_f = out$csAct_f,
                  gva_f = out$gva_f,
                  gvaAct_f = out$gvaAct_f,
                  ccwCr_f = out$ccwCr_f,
                  wagen_f = out$wagen_f,
                  gcf_f = out$gcf_f,
                  gcfAct_f = out$gcfAct_f,
                  gp_f = out$gp_f,
                  ps_f = out$ps_f,
                  psAct_f = out$psAct_f,
                  sts_f = out$sts_f,
                  stsAct_f = out$stsAct_f)
  ))


} else {                                          #objet de classe 'iamOutput'

  if (is.na(desc)) desc <- "My iamOutput object"

    return(new("iamOutput", desc=desc, arguments=objArgsIni@arguments, specific=objInput@specific,
		           outputSp = list(
                    F = out$F,
                    Fr = out$Fr,
                    Fothi = out$Fothi,
                    Fbar = out$Fbar,
                    Z = out$Z,
                    N = out$N,
                    B = out$B,
                    SSB = out$SSB,
                    C = out$C,
                    Ctot = out$Ctot,
                    Y = out$Y,
                    Ytot = out$Ytot,
                    D = out$D,
                    Li = out$Li,
                    Lc = out$Lc,
                    Ltot = out$Ltot,
                    L_et= out$L_et,
                    L_pt = out$L_pt,
                    P = out$P,
                    GVL_f_m_e = out$E$GVL_f_m_e_out,
                     GVLcom_f_m_e = out$E$GVLcom_f_m_e_out,
                     GVLst_f_m_e = out$E$GVLst_f_m_e_out,
                    statY = out$Ystat,
                    statL = out$Lstat,
                    statD = out$Dstat,
                    statP = out$Pstat,
                    statGVL_f_m = out$E$GVL_f_m_eStat_out,
                     statGVLcom_f_m = out$E$GVLcom_f_m_eStat_out,
                     statGVLst_f_m = out$E$GVLst_f_m_eStat_out,
                    PQuot = out$PQuot,
                    TradedQ = out$TradedQ_f,
                    # F_G1 = out$F_G1,F_G2 = out$F_G2,
                    # Fr_G1 = out$Fr_G1,Fr_G2 = out$Fr_G2,
                    # Fothi_G1 = out$Fothi_G1,Fothi_G2 = out$Fothi_G2,
                    # N_G1 = out$N_G1,N_G2 = out$N_G2,
                    # Z_G1 = out$Z_G1,Z_G2 = out$Z_G2,
                    # C_G1 = out$C_G1,C_G2 = out$C_G2,
                    # Ctot_G1 = out$Ctot_G1,Ctot_G2 = out$Ctot_G2,
                    F_S1M1= out$F_S1M1,F_S1M2= out$F_S1M2,F_S1M3= out$F_S1M3,F_S1M4= out$F_S1M4,
                    F_S2M1= out$F_S2M1,F_S2M2= out$F_S2M2,F_S2M3= out$F_S2M3,F_S2M4= out$F_S2M4,
                    F_S3M1= out$F_S3M1,F_S3M2= out$F_S3M2,F_S3M3= out$F_S3M3,F_S3M4= out$F_S3M4,
                    F_S4M1= out$F_S4M1,F_S4M2= out$F_S4M2,F_S4M3= out$F_S4M3,F_S4M4= out$F_S4M4,
                    Fr_S1M1= out$Fr_S1M1,Fr_S1M2= out$Fr_S1M2,Fr_S1M3= out$Fr_S1M3,Fr_S1M4= out$Fr_S1M4,
                    Fr_S2M1= out$Fr_S2M1,Fr_S2M2= out$Fr_S2M2,Fr_S2M3= out$Fr_S2M3,Fr_S2M4= out$Fr_S2M4,
                    Fr_S3M1= out$Fr_S3M1,Fr_S3M2= out$Fr_S3M2,Fr_S3M3= out$Fr_S3M3,Fr_S3M4= out$Fr_S3M4,
                    Fr_S4M1= out$Fr_S4M1,Fr_S4M2= out$Fr_S4M2,Fr_S4M3= out$Fr_S4M3,Fr_S4M4= out$Fr_S4M4,
                    Z_S1M1= out$Z_S1M1,Z_S1M2= out$Z_S1M2,Z_S1M3= out$Z_S1M3,Z_S1M4= out$Z_S1M4,
                    Z_S2M1= out$Z_S2M1,Z_S2M2= out$Z_S2M2,Z_S2M3= out$Z_S2M3,Z_S2M4= out$Z_S2M4,
                    Z_S3M1= out$Z_S3M1,Z_S3M2= out$Z_S3M2,Z_S3M3= out$Z_S3M3,Z_S3M4= out$Z_S3M4,
                    Z_S4M1= out$Z_S4M1,Z_S4M2= out$Z_S4M2,Z_S4M3= out$Z_S4M3,Z_S4M4= out$Z_S4M4,
                    N_S1M1= out$N_S1M1,N_S1M2= out$N_S1M2,N_S1M3= out$N_S1M3,N_S1M4= out$N_S1M4,
                    N_S2M1= out$N_S2M1,N_S2M2= out$N_S2M2,N_S2M3= out$N_S2M3,N_S2M4= out$N_S2M4,
                    N_S3M1= out$N_S3M1,N_S3M2= out$N_S3M2,N_S3M3= out$N_S3M3,N_S3M4= out$N_S3M4,
                    N_S4M1= out$N_S4M1,N_S4M2= out$N_S4M2,N_S4M3= out$N_S4M3,N_S4M4= out$N_S4M4,
                    DD_efmi= out$DD_efmi,
                    DD_efmc= out$DD_efmc,
                    LD_efmi= out$LD_efmi,
                    LD_efmc= out$LD_efmc,
                    statDD_efm= out$statDD_efm,
                    statLD_efm= out$statLD_efm,
                    statLDst_efm= out$statLDst_efm,
                    statLDor_efm= out$statLDor_efm,
                    oqD_ef= out$oqD_ef,
                    oqD_e= out$oqD_e,
                    oqDstat_ef= out$oqDstat_ef,
                    TACtot = out$TACtot,
                    TACbyF = out$TACbyF,
                    PQuot_conv = out$PQuot_conv,
                    diffLQ_conv = out$diffLQ_conv),
                output = list(
                  # typeGest = out$typeGest,
                  nbv_f = out$Eff$nbv_f,
                  effort1_f = out$Eff$effort1_f,
                  effort2_f = out$Eff$effort2_f,
                  nbv_f_m = out$Eff$nbv_f_m,
                  effort1_f_m = out$Eff$effort1_f_m,
                  effort2_f_m = out$Eff$effort2_f_m,
                  #Lbio_f = out$E$Lbio_f,
                  GVLtot_f_m = out$E$GVLtot_f_m_out,
                  GVLav_f_m = out$E$GVLav_f_m_out,
                  GVLtot_f = out$E$GVLtot_f_out,
                  GVLav_f = out$E$GVLav_f_out,
                  #GVLoths_f = out$GVLoths_f,
                  NGVLav_f_m = out$E$NGVLav_f_m_out,
                  NGVLav_f = out$E$NGVLav_f_out,
                  ET_f_m = out$E$ET_f_m_out,
                  cnb_f_m = out$E$cnb_f_m_out,
                  cnb_f = out$E$cnb_f_out,
                  #vcst_f_m = out$E$vcst_f_m,
                  #vcst_f = out$E$vcst_f,
                  rtbs_f_m = out$E$rtbs_f_m_out,
                  rtbs_f = out$E$rtbs_f_out,
                  rtbsAct_f = out$E$rtbsAct_f_out,
                  cshrT_f_m = out$E$cshrT_f_m_out,
                  cshrT_f = out$E$cshrT_f_out,
                  ncshr_f = out$E$ncshr_f_out,
                  ocl_f = out$E$ocl_f_out,
                  cs_f = out$E$cs_f_out,
                  csAct_f = out$E$csAct_f_out,
                  csTot_f = out$E$csTot_f_out,
                  gva_f = out$E$gva_f_out,
                  gvaAct_f = out$E$gvaAct_f_out,
                  gvamargin_f = out$E$gvamargin_f_out,
                  gva_FTE_f = out$E$gva_FTE_f_out,
                  ccw_f = out$E$ccw_f_out,
                  ccwCr_f = out$E$ccwCr_f_out,
                  wageg_f = out$E$wageg_f_out,
                  wagen_f = out$E$wagen_f_out,
                  wageg_FTE_f = out$E$wageg_FTE_f_out,
                  wageg_h_f = out$E$wageg_h_f_out,
                  gp_f = out$E$gp_f_out,
                  gpAct_f = out$E$gpAct_f_out,
                  gpmargin_f = out$E$gpmargin_f_out,
                  ncf_f = out$E$ncf_f_out,
                  np_f = out$E$np_f_out,
                  npmargin_f = out$E$npmargin_f_out,
                  prof_f = out$E$prof_f_out,
                  npmargin_trend_f = out$E$npmargin_trend_f_out,
                  ssTot_f = out$E$ssTot_f_out,
                  ps_f = out$E$ps_f_out,
                  psAct_f = out$E$psAct_f_out,
                  sts_f = out$E$sts_f_out,
                  stsAct_f = out$E$stsAct_f_out,
                  BER_f = out$E$BER_f_out,
                  CR_BER_f = out$E$CR_BER_f_out,
                  fuelEff_f = out$E$fuelEff_f_out,
                  ratio_fvol_gva_f = out$E$ratio_fvol_gva_f_out,
                  ratio_gp_gva_f = out$E$ratio_gp_gva_f_out,
                  ratio_GVL_K_f = out$E$ratio_GVL_K_f_out,
                  ratio_gp_K_f = out$E$ratio_gp_K_f_out,
                  RoFTA_f = out$E$RoFTA_f_out,
                  ROI_f = out$E$ROI_f_out,
                  ratio_np_K_f = out$E$ratio_np_K_f_out,
                  ratio_GVL_cnb_ue_f = out$E$ratio_GVL_cnb_ue_f_out,
                  YTOT_fm = out$YTOT_fm,
                  reconcilSPP = out$reconcilSPP,
                  quotaExp_f = out$E$QuotaExp_f_out,
                  allocEff_f_m = out$allocEff_fm,
                  GoFish = out$GoFish)
  ))
}

})



#:::::::::::::::::::::::::::::::
#Examples
#:::::::::::::::::::::::::::::::


#out <- IAM.input("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Inputs_SIAD_SEL_2.xls",t_init=2010,nbStep=21)
#
#arg <- IAM.args(out)
#
#mod <- IAM.model(arg,out)
#
