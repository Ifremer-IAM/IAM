#fonction de conversion des inputs au niveau m?tier Eco (inclus impl?mentation de l'allocation de mortalite par peche

#' convertInput
#'
#' Fonction de ventilation du jeu de donnee
#'
#' @param inp # TODO
#' @param Fq_fmi # TODO
#' @param Fg_fmi # TODO
#' @param verbose TRUE will print information about the function process.
#' Implemented for debug purpose only.
#'
#' @importFrom abind adrop
#'
convertInput <- function(inp,Fq_fmi=NULL, Fg_fmi=NULL, verbose = FALSE) {

  namF <- inp@specific$Fleet ; nF <- length(namF)
  namM <- inp@specific$Metier ; nM <- length(namM)
  llF <- list() # TODO : initialise with good length and names

  #1ere etape : ventilation de la mortalite (ATTENTION : ici, la ventilation sur indices non communs n'est pas envisag?e)
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (i in inp@specific$Species) {
    if(verbose) cat("--", i, ' :')
    if (inp@specific$Q[i]==0 & inp@specific$S[i]==0) {
      if(verbose) cat(' is XSA')
      namI <- inp@specific$Ages[[i]] ; nI <- length(namI)

      Fini <- inp@input[[i]]$F_fmi
      Fmi <- array(NA, dim=c(nM,nI), dimnames=list(namM,namI))
      Ffmi <- array(NA, dim=c(nF,nM,nI), dimnames=list(namF,namM,namI))

      Cmi <- inp@input[[i]]$C_mi
      Ci <- inp@input[[i]]$C_i
      Ymi <- inp@input[[i]]$Y_mi
      Yi <- inp@input[[i]]$Y_i
      Lref <- inp@input[[i]]$Lref_f_e
      Lrefm <- inp@input[[i]]$Lref_f_m_e
      FM <- inp@input[[i]]$fm

      if (attributes(Fini)$DimCst[1]>0 & attributes(Fini)$DimCst[2]>0) {
        Ffmi[] <- Fini[]
      } else {
        #ventilation m?tier
        if (attributes(Fini)$DimCst[2]==0) {                                         #ie F = Fi (cas 1, 2)
          if (!all(is.na(Cmi)) & !all(is.na(Ci)) & all(attributes(Cmi)$DimCst[2:3]>0) & attributes(Ci)$DimCst[3]>0) {   #ie C_mi renseign? avec composante m?tier et ?ge, et C_i renseign? avec composante ?ge -> cas 1 ou 2
            if (attributes(Cmi)$DimCst[1]>0) {                                        #ie C_mi = Cfmi  (cas 2)
              Ffmi[] <- Cmi[] # TODO : is it usefull ? add some time.
              aggCmi <- apply(Cmi,3,sum,na.rm=TRUE) ; Ci[aggCmi>Ci] <- aggCmi[aggCmi>Ci]  #on remplace dans Ci les valeurs agr?g?es issues de C_mi sup?rieures
              Ffmi <- Ffmi*rep(Fini/Ci,each=nF*nM)
            } else {                                                                  #ie C_mi = Cmi (cas 1)
              Fmi[] <- Cmi[] # TODO : useless line, write Fmi as Cmi below
              aggCmi <- apply(Cmi,2,sum,na.rm=TRUE) # TODO : colSums
              Ci[aggCmi>Ci] <- aggCmi[aggCmi>Ci]
              Fmi <- Fmi*rep(Fini/Ci,each=nM)
              # TODO : replace by Fmi*matrix(Fini/Ci, ncol = nI, nrow = nM, byrow = TRUE)
            }
          } else { if (!all(is.na(Ymi)) & !all(is.na(Yi)) & all(attributes(Ymi)$DimCst[2:3]>0) & attributes(Yi)$DimCst[3]>0) {   #ie Y_mi renseign? avec composante m?tier et ?ge, et Y_i renseign? avec composante ?ge -> cas 1 ou 2
            if (attributes(Ymi)$DimCst[1]>0) {                                #ie Y_mi = Yfmi  (cas 2)
              Ffmi[] <- Ymi[]
              aggYmi <- apply(Ymi,3,sum,na.rm=TRUE) ; Yi[aggYmi>Yi] <- aggYmi[aggYmi>Yi]
              Ffmi <- Ffmi*rep(Fini/Yi,each=nF*nM)
            } else {                                                          #ie Y_mi = Ymi (cas 1)
              Fmi[] <- Ymi[]
              aggYmi <- apply(Ymi,2,sum,na.rm=TRUE) ; Yi[aggYmi>Yi] <- aggYmi[aggYmi>Yi]
              Fmi <- Fmi*rep(Fini/Yi,each=nM)
            }
          } else {                   #
            Fmi[] <- Fini[]          #  added 18/07/16
          }                          #
          }
        } else {
          Fmi[] <- Fini[] # TODO : remonter ce cas dans le premier if.
        }

        #on poursuit avec la ventilation flottille si Ffmi non compl?t? (on suppose ? ce stade que Fmi a ?t? compl?t?)

        if (all(is.na(Ffmi))) {                                                       #ie cas 3 ou 4
          if (!all(is.na(Cmi)) & !all(is.na(Ci)) & all(attributes(Cmi)$DimCst[1:3]>0) & all(attributes(Ci)$DimCst[2:3]>0)) {   #ie cas 4
            Ffmi[] <- Cmi[]
            aggCmi <- apply(Cmi,2:3,sum,na.rm=TRUE) ; Ci[aggCmi>Ci] <- aggCmi[aggCmi>Ci]
            Ffmi <- Ffmi*rep(Fmi/Ci,each=nF)
          } else {
            if (!all(is.na(Ymi)) & !all(is.na(Yi)) & all(attributes(Ymi)$DimCst[1:3]>0) & all(attributes(Yi)$DimCst[2:3]>0)) {   #ie cas 4
              Ffmi[] <- Ymi[]
              aggYmi <- apply(Ymi,2:3,sum,na.rm=TRUE) ; Yi[aggYmi>Yi] <- aggYmi[aggYmi>Yi]
              Ffmi <- Ffmi*rep(Fmi/Yi,each=nF)
            } else {                                                                # il faut alors utiliser les donn?es de d?barquements des feuillets Eco, redistribu?es par m?tierBio via la matrice fm
              if (!all(is.na(Lref)) & !all(is.na(FM)) & attributes(Ymi)$DimCst[1]==0 & attributes(Ymi)$DimCst[2]>0) { #cas 3 avec Ctot_m calcul? ? partir de Y_mi
                CtotM <- apply(Ymi,1,sum,na.rm=TRUE) # TODO : rowSums
                Cfm <- FM*as.vector(Lref)
                aggC <- apply(Cfm,2,sum,na.rm=TRUE) ; CtotM[aggC>CtotM] <- aggC[aggC>CtotM]
                Ffmi[] <- Cfm/rep(CtotM,each=nF)
                Ffmi <- Ffmi*rep(Fmi,each=nF)
              } else {
                if (!all(is.na(Lref)) & !all(is.na(FM)) & attributes(Yi)$DimCst[1]==0 & attributes(Yi)$DimCst[2]>0) { #cas 3 avec Ctot_m calcul? ? partir de Y_i
                  CtotM <- apply(Yi,1,sum,na.rm=TRUE)
                  Cfm <- FM*as.vector(Lref)
                  aggC <- apply(Cfm,2,sum,na.rm=TRUE) ; CtotM[aggC>CtotM] <- aggC[aggC>CtotM]
                  Ffmi[] <- Cfm/rep(CtotM,each=nF)
                  Ffmi <- Ffmi*rep(Fmi,each=nF)
                } else {    # cas ultime : on utilise les d?barquements f*m de r?f?rence dans les feuillets Fleet sur Fi et Ytot calcul? ? partir de Yi
                  if (attributes(Fini)$DimCst[1]==0 & attributes(Fini)$DimCst[2]==0 & !all(is.na(Lrefm)) & !all(is.na(Yi))) {
                    Ffmi[] <- (Lrefm/sum(Yi,na.rm=TRUE))%o%Fini            # fait le 26/07/2013
                  }
                }
              }
            }
          }
        }
      }

      # ici, Ffmi devrait etre dispo
      if(verbose) cat(' check')
      if (all(is.na(Ffmi))) stop("wrong or missing data for F allocation! Check C_mi, C_i, Y_mi, Y_i, or Lref_f_m inputs!!")

      llF[[i]] <- Ffmi

    } else if (inp@specific$Q[i]==1 & inp@specific$S[i]==0){
      if(verbose) cat(' is SS3')
      llF[[i]] <- adrop(Fq_fmi[[i]][1,1,,,,drop=FALSE],1:2)
    } else if (inp@specific$Q[i]==0 & inp@specific$S[i]==1){
      if(verbose) cat(' is SEX')
      llF[[i]] <- adrop(Fg_fmi[[i]][1,,,,drop=FALSE],1)
    }
  }

  # modif
  #2eme etape : on redefinit chaque variable BIO suivant le niveau metier 'ECO' ? l'aide de la matrice MM
  # -> ne concerne que la variable d_i  <-- peut ?voluer avec l'?volution de 'sr' et l'int?gration des variables de prix au niveau m?tier
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (i in inp@specific$Species) {

    namI <- inp@specific$Ages[[i]] ; nI <- length(namI)
    namME <- inp@specific$MetierEco ; nME <- length(namME)

    MM <- inp@input[[i]]$mm
    if (is.null(MM)) MM <- NA
    if (length(MM)==1) tabMM <- NA else {
      tabMM <- cbind.data.frame(expand.grid(dimnames(MM)),value2=as.vector(MM))
      names(tabMM) <- c("fm","mEco","val2") }

    #conversion des donn?es de mortalit?s --> on utilise les valeurs brutes
    tabF <- cbind.data.frame(expand.grid(dimnames(llF[[i]])),value1=as.vector(llF[[i]])) ; names(tabF) <- c("f","m","a","val1")
    tabF$fm <- paste(tabF$f,tabF$m,sep="__")
    if (all(is.na(tabMM))) TABF <- cbind(tabF,mEco=tabF$m,val2=as.vector(MM)) else TABF <- merge(tabMM,tabF,all=TRUE)
    TABF$val <- TABF$val1*TABF$val2 ; TABF <- TABF[!is.na(TABF$val),]
    TABF$f <- factor(as.character(TABF$f),levels=namF)
    TABF$mEco <- factor(as.character(TABF$mEco),levels=namME)
    TABF$a <- factor(as.character(TABF$a),levels=namI)
    FF <- with(TABF,tapply(val,list(f,mEco,a),function(x) x))
    attributes(FF)$DimCst <- as.integer(c(nF,nME,nI,0))

    if (attributes(inp@input[[i]]$F_fmi)$DimCst[1]>0 & attributes(inp@input[[i]]$F_fmi)$DimCst[2]>0) {
      inp@input[[i]]$F_i <- apply(inp@input[[i]]$F_fmi,3,sum,na.rm=TRUE)
      attributes(inp@input[[i]]$F_i)$DimCst <- as.integer(c(0,0,nI,0))
    } else {
      if (attributes(inp@input[[i]]$F_fmi)$DimCst[2]>0) {
        inp@input[[i]]$F_i <- apply(inp@input[[i]]$F_fmi,2,sum,na.rm=TRUE)
        attributes(inp@input[[i]]$F_i)$DimCst <- as.integer(c(0,0,nI,0))
      } else {
        inp@input[[i]]$F_i <- inp@input[[i]]$F_fmi
      }
    }

    inp@input[[i]]$F_fmi <- FF

    #conversion des variables ratios --> on applique les m?mes valeurs pour tous les m?tiers correspondants (on les suppose d?finie aux ?ges)
    if(inp@specific$Q[i]==0 & inp@specific$S[i]==0) {
      di <- inp@input[[i]]$d_i
      if (attributes(di)$DimCst[2]>0) {
        if (attributes(di)$DimCst[1]==0) {
          vec_di <- llF[[i]] ; vec_di[] <- NA #matrice au format fmi
          vec_di[] <- rep(di,each=nF)
        } else {
          vec_di <- di
        }

        tabD <- cbind.data.frame(expand.grid(dimnames(vec_di)),value1=as.vector(vec_di)) ; names(tabD) <- c("f","m","a","val1")
        tabD$fm <- paste(tabD$f,tabD$m,sep="__")
        if (all(is.na(tabMM))) TABD <- cbind(tabD,mEco=tabD$m,val2=MM) else TABD <- merge(tabMM,tabD,all=TRUE)
        #TABD$val <- TABD$val1*TABD$val2  change Florence 05/2019
        TABD$val <- TABD$val1; TABD <- TABD[!is.na(TABD$val),]
        TABD$f <- factor(as.character(TABD$f),levels=namF)
        TABD$mEco <- factor(as.character(TABD$mEco),levels=namME)
        TABD$a <- factor(as.character(TABD$a),levels=namI)
        FD <- with(TABD,tapply(val,list(f,mEco,a),function(x) x))
        attributes(FD)$DimCst <- as.integer(c(nF,nME,nI,0))

        inp@input[[i]]$d_i <- FD
      }
    } else if(FALSE){
      di_G1 <- inp@input[[i]]$d_i_G1 #G1
      if (attributes(di_G1)$DimCst[2]>0) {
        if (attributes(di_G1)$DimCst[1]==0) {
          vec_di <- llF[[i]] ; vec_di[] <- NA #matrice au format fmi
          vec_di[] <- rep(di_G1,each=nF)
        } else {
          vec_di <- di_G1
        }

        tabD <- cbind.data.frame(expand.grid(dimnames(vec_di)),value1=as.vector(vec_di)) ; names(tabD) <- c("f","m","a","val1")
        tabD$fm <- paste(tabD$f,tabD$m,sep="__")
        if (all(is.na(tabMM))) TABD <- cbind(tabD,mEco=tabD$m,val2=MM) else TABD <- merge(tabMM,tabD,all=TRUE)
        TABD$val <- TABD$val1*TABD$val2 ; TABD <- TABD[!is.na(TABD$val),]
        TABD$f <- factor(as.character(TABD$f),levels=namF)
        TABD$mEco <- factor(as.character(TABD$mEco),levels=namME)
        TABD$a <- factor(as.character(TABD$a),levels=namI)
        FD <- with(TABD,tapply(val,list(f,mEco,a),function(x) x))
        attributes(FD)$DimCst <- as.integer(c(nF,nME,nI,0))

        inp@input[[i]]$d_i_G1 <- FD
      }

      di_G2 <- inp@input[[i]]$d_i_G2 #G2
      if (attributes(di_G2)$DimCst[2]>0) {
        if (attributes(di_G2)$DimCst[1]==0) {
          vec_di <- llF[[i]] ; vec_di[] <- NA #matrice au format fmi
          vec_di[] <- rep(di_G2,each=nF)
        } else {
          vec_di <- di_G2
        }

        tabD <- cbind.data.frame(expand.grid(dimnames(vec_di)),value1=as.vector(vec_di)) ; names(tabD) <- c("f","m","a","val1")
        tabD$fm <- paste(tabD$f,tabD$m,sep="__")
        if (all(is.na(tabMM))) TABD <- cbind(tabD,mEco=tabD$m,val2=MM) else TABD <- merge(tabMM,tabD,all=TRUE)
        TABD$val <- TABD$val1*TABD$val2 ; TABD <- TABD[!is.na(TABD$val),]
        TABD$f <- factor(as.character(TABD$f),levels=namF)
        TABD$mEco <- factor(as.character(TABD$mEco),levels=namME)
        TABD$a <- factor(as.character(TABD$a),levels=namI)
        FD <- with(TABD,tapply(val,list(f,mEco,a),function(x) x))
        attributes(FD)$DimCst <- as.integer(c(nF,nME,nI,0))

        inp@input[[i]]$d_i_G2 <- FD
      }
    }
  }

  return(inp)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# standFormat
#
# reformat a dataframe into a standart formated matrix
#
# @param DF a data frame
# @param nbStep Number of year to use in the model
#
standFormat <- function(DF,nbStep,modF,modM,modI,modC,alk,as.na=NULL) {

  if (is.null(ncol(DF))) {
    Mat <- as.numeric(as.character(DF))
  } else {

    if (ncol(DF)==1) {
      Mat <- as.numeric(as.character(DF$value))
      attributes(Mat)$DimCst <- as.integer(rep(0,4))
    } else {
      #  if ("t" %in% names(DF)) DF <- expand.time(DF,nbStep)

      if ("f" %in% names(DF) & !all(is.na(modF))) {
        DF$f <- factor(as.character(DF$f),levels=modF)
        dimL.f <- length(modF)
      } else {
        dimL.f <- 0
      }

      if ("m" %in% names(DF) & !all(is.na(modM))) {
        DF$m <- factor(as.character(DF$m),levels=modM)
        dimL.m <- length(modM)
      } else {
        dimL.m <- 0
      }

      if ("i" %in% names(DF)  & !all(is.na(modI))) {
        DF$i <- factor(as.character(DF$i),levels=modI)
        dimL.i <- length(modI)
      } else {
        dimL.i <- 0
      }

      if ("l" %in% names(DF)  & !all(is.na(DF))) {
        DF$l <- factor(as.character(DF$l),levels=paste0("l__",dimnames(alk)[[1]]))
        dimL.l <- 1
      } else {
        dimL.l <- 0
      }

      if ("c" %in% names(DF)  & !all(is.na(DF))) {
        DF$c <- factor(as.character(DF$c),levels=modC)
        dimL.c <- length(modC)
      } else {
        dimL.c <- 0
      }

      if ("t" %in% names(DF)) {ve <- as.character(DF$t)
      unve <- unique(ve)
      DF$t <- factor(ve,levels=unve[order(as.numeric(substr(unve,4,1000)))])
      }


      dimL <- c(f=dimL.f, m=dimL.m,             # e=length(unique(DF$e)),
                i=dimL.i, t=length(unique(DF$t)), c=dimL.c)

      eval(parse('',text=paste0("Mat <- suppressWarnings(with(DF,tapply(as.numeric(as.character(value)),list(",
                               paste(c("l"[dimL.l>0],"f"[dimL["f"]>0],"m"[dimL["m"]>0],#"e"[dimL["e"]>0],  #on organisera les modules bio par esp?ce,
                                       "i"[dimL["i"]>0],"c"[dimL["c"]>0],"t"[dimL["t"]>0]),collapse=","),  #et regroup?s en liste apr?s traitement
                               "),function(x) x)))")))

      if (!is.null(as.na)) Mat[is.na(Mat)] <- as.na
      #on enl?ve des en-t?tes les pr?fixes indicateurs
      dimnames(Mat) <- lapply(dimnames(Mat), function(x) sapply(1:length(x),function(y) substring(x[y],4,nchar(x[y]))))

      #on applique la cl? taille-?ge si besoin
      if ((!is.null(alk)) & (dimL.l>0)) {
        Mat <- t(alk/apply(alk,1,sum,na.rm=TRUE))%*%Mat
        Mat <- aperm(Mat, match(c("f"[dimL["f"]>0],"m"[dimL["m"]>0],"i"[dimL.l>0],"t"[dimL["t"]>0]),
                                c("i"[dimL.l>0],"f"[dimL["f"]>0],"m"[dimL["m"]>0],"t"[dimL["t"]>0])))
        dimL[3] = length(modI)
      }

      attributes(Mat)$DimCst <- as.integer(dimL[1:4]) ; if (dimL.c>0) {
        if (dimL.i>0) {
          attributes(Mat)$DimCst <- NULL
        } else {
          attributes(Mat)$DimCst[3] <- as.integer(dimL.c)
        }
      }
    }
  }
  return(Mat)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

twoDto1D <- function(df,dim="2D") {

  if (dim=="2D") {
    df <- as.matrix(df) ; df[is.na(df)] <- ""
    nColHead <- sum(df[1,]=="")
    dfTemp <- expand.grid(apply(df[-1,1:nColHead,drop=FALSE],1,paste,collapse=":@:@:"),df[1,-(1:nColHead)])

    if (nColHead==1){
      DF <- cbind.data.frame(as.character(dfTemp[,1]), as.character(dfTemp[,2]), as.numeric(df[-1,-1]))

    } else {
      DF <- cbind.data.frame(do.call("rbind",lapply(as.character(dfTemp[,1]), function(x) strsplit(x,":@:@:")[[1]])),
                             as.character(dfTemp[,2]), as.numeric(df[-1,-(1:nColHead)]))
    }
  } else {
    DF <- df
  }

  names(DF) <- c(substring(apply(DF,1,as.character)[1:(ncol(DF)-1)],1,1),"value")
  return(DF)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#extrapole les valeurs aux temps non d?crits pour chaque pas de temps
#(attention, n'est valable que si pr?sence de champs "t" et "value" dans df

expand.time <- function(df,t_init,nbStep=1,scenario=FALSE){

  if (is.null(nbStep)) stop("nbStep parameter is NULL!!")
  #modalit?s de temps
  occ = paste0("t__",seq(t_init,length=nbStep))
  #les modalit?s qui ?volueront au cours du temps
  mod = unique(df[,-match(c("t","value"),names(df)),drop=FALSE])
  DF <- NULL ; TAB = NULL
  for (i in 1:nbStep) {
    if (ncol(df)>2){

      #on cr?? la portion de table qu'il faudra remplir
      port <- cbind.data.frame(mod,data.frame(t=occ[i]))
      #on merge avec la partie de la table en commun
      tab <- merge(port,df,all.x=TRUE)
      #s'il reste des NA, on va chercher les valeurs de l'instant (t-1) -> table TAB
      if (any(is.na(tab$value)))
      {
        if (scenario) {
          #tab$value[is.na(tab$value)] <- 1
        } else {

          tab$value[is.na(tab$value)] <- merge(
            tab[is.na(tab$value),-match(c("t","value"),names(tab)),drop=FALSE],TAB,all.x=TRUE)$value
        }}

    } else { #ie seulement les champs 'temps' et 'valeurs'
      tab = merge(df,data.frame(t=occ[i]),all.y=TRUE)
      if (is.na(tab$value)) tab$value <- TAB$value
    }
    TAB <- tab
    DF <- rbind(DF,TAB)
  }
  return(DF)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#on aura besoin des objets d?finis
#source("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Input_object.r")

#' @importFrom methods new
reformat <- function(x, slotN="stockInput") {
  n <- names(new(slotN)@input)
  ll <- x[n]
  names(ll) <- n
  ll <- lapply(ll,function(y) {if (is.null(y)) return(as.numeric(NA)) else {if (all(is.na(y))) return(as.numeric(NA)) else return(y)}})
  return(lapply(ll,function(y){if (length(y)==1 & is.null(names(y)) & is.null(dimnames(y)) & is.numeric(y)) attributes(y)$DimCst <- as.integer(c(0,0,0,0))   #added 18/07/16
  return(y)}))
}    #si on veut plut?t des NAs, on remplace 'll)}' par 'lapply(ll,function(y) if (is.null(y)) NA else y))}'


#' read.Pflex
#'
#' Function to import a Price_flexibility sheet from the input file.
#'
#' @param file excell file (.xlsx) imported by IAM.input.
#' Must contain a Price_flexibility sheet.Character string.
#' @param nam_stock vector of character string with stock names.
#' Format is "Stock__X", X being an abreviation for a precise stock.
#' @param nam_stock_bis # TODO what is this format
#'
#' # TODO add example file.
#'
#' @importFrom reshape2 acast
#'
#' @author Florence Briton 05/2019
read.Pflex <- function(file, nam_stock, nam_stock_bis){
  namList <- gsub("Stock__","",nam_stock)
  indDYN <- length(nam_stock)>0

  PFlex <- read.sheet(file = file, sheet = "Price_flexibility")

  vec <- as.vector(PFlex)
  vec <- vec[vec!=""]
  MOD=NULL
  MOD[[1]] <- c("NONE",gsub('p__',"",unique(vec[sapply(vec,function(y) substring(y,1,3)=='p__')]))) #modalites produits marche
  Spp_dyn=if (indDYN) as.character(namList) else character(0)
  StaticSpp=as.character(nam_stock_bis)
  MOD[[2]] <- c(Spp_dyn,StaticSpp) #modalites especes

  indEmpt <- suppressWarnings(apply(PFlex,1,function(x) min(unlist(sapply(c("v__","e__","p__"),grep,x)))))
  invisible(sapply(1:nrow(PFlex),function(x) if (is.finite(indEmpt[x])) {if (indEmpt[x]>1) PFlex[x,1:(indEmpt[x]-1)] <<- ""}))

  #on filtre tout ce qui n'est ni numerique, ni parametre
  #conversion en numerique
  num <- apply(suppressWarnings(apply(PFlex,1,as.numeric)),1,as.character)

  #on ajoute les parametres
  indic <- substring(PFlex,1,3)%in%c("v__","e__","p__")
  num[indic] <- PFlex[indic]

  indicRow <- apply(PFlex,1,function(x) any(substring(x,1,3)%in%c("v__","e__","p__")))
  indicTbl <- cumsum(apply(num,1,function(x) all(is.na(x))))
  #on separe les tables (sauts de lignes)
  sepTabl <- split(as.data.frame(num)[indicRow,],indicTbl[indicRow])
  tbl <- lapply(sepTabl,function(x) x[,!apply(x,2,function(y) all(is.na(y)))])

  #il faut maintenant filtrer toutes les anomalies de format
  #on distingue pour commencer les tables 1D des tables 2D
  tbl2Dind <- lapply(tbl,function(x) !(substring(as.character(x[1,1]),1,3)%in%c("v__","e__","p__")))
  tbl2D <- tbl[(1:length(tbl))[unlist(tbl2Dind)]]

  if (length(tbl2D)>0) {tbl2D <- lapply(tbl2D,function(x) x[,apply(x,2,function(y) any(substring(as.matrix(y),1,3)%in%c("v__","e__","p__")))])
  invisible(lapply(1:length(tbl2D),function(x) tbl2D[[x]][tbl2D[[x]]==-1] <<- as.numeric(NA)))}

  #on peut maintenant separer les variables
  #pour cela, il faut tout mettre sous forme 1D
  if (length(tbl2D)>0) tbl2 <- lapply(tbl2D,twoDto1D,"2D") else tbl2 <- NULL
  ListPflex <- tbl2
  ListPflex <- lapply(ListPflex,function(x) split(x[,-match("v",names(x)),drop=FALSE],as.character(x[,match("v",names(x))])))
  namL <- gsub("v__","",unlist(lapply(ListPflex,names)))
  ListPflex <- unlist(ListPflex,recursive=FALSE,use.names = FALSE)
  names(ListPflex) <- namL

  #il faut regrouper les tables de meme variable decoupees en plusieurs parties
  Nam <- unique(names(ListPflex))
  ListPflex <- lapply(Nam,function(x) {tt <- do.call("rbind",ListPflex[names(ListPflex)%in%x])
                                       rownames(tt) <- NULL
                                       return(tt)})
  names(ListPflex) <- Nam

  #il faut etendre la matrice ep a toutes les especes
  combi = expand.grid(e=paste0('e__',MOD[[2]]),p=paste0('p__',MOD[[1]]))
  ep = ListPflex$ep
  ep=merge(combi,ep,all.x = TRUE)
  ep[is.na(ep)]=0
  ListPflex$ep = ep

  combi2 = expand.grid(p=paste0('p__',MOD[[1]]),p.1=paste0('p__',MOD[[1]]))
  beta_pp = ListPflex$beta_pp
  beta_pp = merge(combi2, beta_pp, all.x=TRUE)
  beta_pp[is.na(beta_pp)]=0
  ListPflex$beta_pp = beta_pp

  ListPflex$beta_pp$p = gsub(x=as.character(ListPflex$beta_pp$p),pattern='p__',replacement='')
  ListPflex$beta_pp$p.1 = gsub(x=as.character(ListPflex$beta_pp$p.1),pattern='p__',replacement='')
  ListPflex$ep$p = gsub(x=as.character(ListPflex$ep$p),pattern='p__',replacement='')
  ListPflex$ep$e = gsub(x=as.character(ListPflex$ep$e),pattern='e__',replacement='')

  ListPflex$beta_pp$p <- factor(as.character(ListPflex$beta_pp$p),levels=MOD[[1]])
  ListPflex$beta_pp$p.1 <- factor(as.character(ListPflex$beta_pp$p.1),levels=MOD[[1]])
  ListPflex$ep$p <- factor(as.character(ListPflex$ep$p),levels=MOD[[1]])
  ListPflex$ep$e <- factor(as.character(ListPflex$ep$e),levels=MOD[[2]])
  ListPflex$ep$value = as.integer(ListPflex$ep$value)
  dim.p <- length(MOD[[1]])
  dim.e <- length(MOD[[2]])
  dimL <- c(p=dim.p, e=dim.e)

  ListPflex$beta_pp = acast(ListPflex$beta_pp,p~p.1,value.var='value')
  attributes(ListPflex$beta_pp)$DimCst = as.integer(c(dimL['p'],dimL['p']))
  ListPflex$ep = acast(ListPflex$ep,e~p,value.var='value')
  attributes(ListPflex$ep)$DimCst = as.integer(c(dimL['e'],dimL['p']))

  ListPflex$ep[which(rowSums(ListPflex$ep)==0),'NONE'] = as.integer(1) #on assigne les especes pour lesquelles pas de produit a produit NONE

  ListPflex$modE = as.character(MOD[[2]])
  ListPflex$modP = as.character(MOD[[1]])

  return(ListPflex)
}


#' read.Scenarii
#'
#' @param file excell file (.xlsx) imported by IAM.input.
#' Must contain a Scenarii sheet.Character string.
#'
read.Scenar <- function(file){

  scenar <- read.sheet(file, "Scenarii")

  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")

  #on ne prend pas en compte les 100 premieres lignes (Attention : format fixe a respecter)
  scenar <- scenar[101:nrow(scenar),]
  # dans le cas du fichier med il y a un soucis car une valeur dans la colonne 13.

  #il faut maintenant tenir compte des scenarios couples ('... & ...') : on duplique afin de n'avoir qu'un scenario par ligne
  repVec <- apply(scenar,1,function(y) length(gregexpr(" & ",as.character(y[1]))[[1]]))
  count <- apply(scenar,1,function(y) grepl(" & ",as.character(y[1])))
  repVec[count] <- repVec[count] + 1
  newSc <- strsplit(as.vector(scenar[,1])," & ")
  newSc <- lapply(newSc,function(x) if (length(x)==0) "" else x)
  scenar <- scenar[rep(1:nrow(scenar),repVec),]
  scenar[,1] <- unlist(newSc)
  scenar[scenar[,1]!="",1] <- paste0("s__",scenar[scenar[,1]!="",1])

  scenar1 <- scenar[,1] ; scenar2 <- scenar[,2:ncol(scenar)]

  indEmpt <- suppressWarnings(apply(scenar2,1,function(x) min(unlist(sapply(prefix,grep,x)))))

  #tout ce qui se trouve avant une modalit? de variable est pass? ? ""
  invisible(sapply(1:nrow(scenar2),function(x) if (is.finite(indEmpt[x])) {if (indEmpt[x]>1) scenar2[x,1:(indEmpt[x]-1)] <<- ""}))

  #on filtre tout ce qui n'est ni num?rique, ni param?tre
  #conversion en num?rique
  num <- apply(suppressWarnings(apply(scenar2,1,as.numeric)),1,as.character)

  #on ajoute les param?tres
  indic <- substring(scenar2,1,3) %in% prefix          #attention : depuis ajout openxlsx, scenar2 --> as.matrix(scenar2)
  num[indic] <- scenar2[indic] ; num[is.na(num)] <- ""  ; scenar <- cbind(scenar1,num)


  #on s?pare les tables
  indicRow <- apply(scenar,1,function(x) any(substring(x,1,3) %in% prefix))
  scenar[!indicRow,] <- rep("",ncol(scenar))
  indicTbl <- cumsum(apply(scenar,1,function(x) all(x=="")))

  sepSc <- split(as.data.frame(scenar)[indicRow,],indicTbl[indicRow])

  ##on colle le pr?fixe ? la colonne sc?nario
  #sepSc <- lapply(sepSc, function(x) {x[x[,1]!="",1] <- paste0("s__",x[x[,1]!="",1])
  #                           return(x)})

  #on distingue pour commencer les tables 1D des tables 2D
  tbl2DindS <- lapply(sepSc,function(x) !substring(as.character(x[1,1]),1,3) %in% c("s__",prefix))
  tbl2DS <- sepSc[(1:length(sepSc))[unlist(tbl2DindS)]]
  tbl1DS <- sepSc[(1:length(sepSc))[!unlist(tbl2DindS)]]

  #r?gles des tables 1d :
  #une seule colonne de num?riques
  if (length(tbl1DS)>0) tbl1DS <- lapply(tbl1DS,function(x) x[,1:((1:ncol(x))[!substring(as.matrix(x[1,]),1,3)%in%c("s__",prefix)][1])])

  #regles des tables 2d :
  #une colonne de numerique doit etre precedee d'une variable
  if (length(tbl2DS)>0) tbl2DS <- lapply(tbl2DS,function(x) x[,apply(x,2,function(y) any(substring(as.matrix(y),1,3)%in%c("s__",prefix)))])

  #on peut maintenant s?parer les variables
  #pour cela, il faut tout mettre sous forme 1D
  if (length(tbl2DS)>0) tbl2S <- lapply(tbl2DS,twoDto1D,"2D") else tbl2S <- NULL
  if (length(tbl1DS)>0) tbl1S <- lapply(tbl1DS,twoDto1D,"1D") else tbl1S <- NULL


  ListS <- c(tbl1S,tbl2S)
  return(ListS)
}

#' read.sheet
#'
#' Import a sheet from an excell file.
#'
#' @param file excell file (.xlsx) imported by IAM.input.
#' Must contain a sheet named after the second parameter.Character string.
#' @param sheet sheet name inside the file imported by IAM.input.
#' Character string.
#'
#' @details
#' Replace "," with "." and all columns into character.
#' Replace NA values with empty string ("")
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom methods rbind2
#'
read.sheet <- function(file, sheet){
  sheet <- read.xlsx(file,sheet=sheet,rowNames=FALSE,colNames=FALSE,skipEmptyRows = FALSE,skipEmptyCols = FALSE)
  sheet[] <- lapply(sheet, function(x) gsub(",",".",as.character(x)))
  sheet <- as.matrix(rbind2("",sheet))
  sheet[is.na(sheet)] <- ""

  return(sheet)
}


#' Reapeated code that I don't understand yet
#'
#' Idea is DRY (Don't Repeat Yourself)
#'
#' @param result a readed sheet from a file, modified before
#' @param indEmpt something I need to understand # TODO
#'
result_filtre <- function(result, indEmpt){

  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")

  #tout ce qui se trouve avant une modalit? de variable est pass? ? ""
  invisible(sapply(1:nrow(result),function(x) if (is.finite(indEmpt[x])) {if (indEmpt[x]>1) result[x,1:(indEmpt[x]-1)] <<- ""}))

  #on filtre tout ce qui n'est ni num?rique, ni param?tre
  #conversion en num?rique

  num <- apply(suppressWarnings(apply(result,1,as.numeric)),1,as.character)

  #on ajoute les param?tres
  indic <- substring(result,1,3) %in% prefix
  num[indic] <- result[indic]

  indicRow <- apply(result,1,function(x) any(substring(x,1,3) %in% prefix))
  indicTbl <- cumsum(apply(num,1,function(x) all(is.na(x))))
  #on s?pare les tables (sauts de lignes)
  sepTabl <- split(as.data.frame(num)[indicRow,],indicTbl[indicRow])
  tbl <- lapply(sepTabl,function(x) x[,!apply(x,2,function(y) all(is.na(y)))])

  #il faut maintenant filtrer toutes les anomalies de format

  #on distingue pour commencer les tables 1D des tables 2D
  tbl2Dind <- lapply(tbl,function(x) !substring(as.character(x[1,1]),1,3) %in% prefix)
  tbl2D <- tbl[(1:length(tbl))[unlist(tbl2Dind)]]
  tbl1D <- tbl[(1:length(tbl))[!unlist(tbl2Dind)]]

  #r?gles des tables 1d :
  #une seule colonne de num?riques

  if (length(tbl1D)>0) {
    tbl1D <- lapply(tbl1D,function(x) x[,1:((1:ncol(x))[!substring(as.matrix(x[1,]),1,3) %in% prefix][1])])
    invisible(lapply(1:length(tbl1D),function(x) tbl1D[[x]][tbl1D[[x]]==-1] <<- as.numeric(NA))) # TODO : quel est le but de cette ligne ? pkoi value == -1 ?
  }

  #r?gles des tables 2d :
  #une colonne de num?rique doit ?tre pr?c?d?e d'une variable

  if (length(tbl2D)>0) {
    tbl2D <- lapply(tbl2D,function(x) x[,apply(x,2,function(y) any(substring(as.matrix(y),1,3) %in% prefix))])
    invisible(lapply(1:length(tbl2D),function(x) tbl2D[[x]][tbl2D[[x]]==-1] <<- as.numeric(NA)))
  }
  #on peut maintenant s?parer les variables
  #pour cela, il faut tout mettre sous forme 1D
  if (length(tbl2D)>0) tbl2 <- lapply(tbl2D,twoDto1D,"2D") else tbl2 <- NULL
  if (length(tbl1D)>0) tbl1 <- lapply(tbl1D,twoDto1D,"1D") else tbl1 <- NULL

  List <- c(tbl1,tbl2)

  return(List)
}


#' init_listHisto
#'
#' initialise list Historique and list Input for a specific stock
#'
#' @param List a result table filtered before
#' @param t_init Initial time of the simulation
#' @param t_hist_max first time of the modelisation
#' @param nbStep number of step for the modelisation
#'
init_listHisto <- function(List, t_init, t_hist_max, nbStep){
  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")
  #on en fait maintenant des objets standards accompagn?s de leur attribut 'DimCst' pour les inputs, et on laisse sous forme de DF pour l'historique
  #il faut consid?rer l'historique... (t<=t_init)
  listHisto <- List[!grepl("s__",names(List))]
  invisible(sapply(prefix,
                   function(y) listHisto <<- lapply(listHisto,function(x) as.data.frame(gsub(y,"",as.matrix(x))))))
  listHisto <- lapply(listHisto, function(x) {if (ncol(x)==1) {
    return(as.numeric(as.character(x$value)))
  } else {
    rownames(x) <- 1:nrow(x)
    if ("t"%in%names(x)) {    #si l'occurence n'est pas pr?sente, on prend toute la table
      rp <- match(as.character(t_hist_max),as.character(x$t))
      if (!is.na(rp)) {
        x <- x[unique(sort(c(1:match(as.character(t_hist_max),as.character(x$t)),
                             (1:nrow(x))[x$t%in%t_hist_max]))),]
      }}
    x$value <- as.numeric(as.character(x$value))
    return(x)
  }})
  #... et les param?tres d'entr?e (t>=t_init)
  listInput <- List[!grepl("s__",names(List))]

  listInput <- lapply(listInput, function(x) {if (ncol(x)==1) {
    return(x)
  } else {
    if ("t"%in%names(x)) {
      #il faut distinguer ce qui va servir ? calculer la valeur initiale (tab), et ce qui sert pour les projections (proj)
      ind <- grep("t__t__",as.character(x$t))
      indic <- length(ind)>0
      #occ <- unique(as.character(x$t)) ; occ <- occ[length(occ)]
      if (indic) tab <- x[ind,] else
        tab <- x[x$t%in%paste0("t__",t_init),]
      if (indic) {

        if (max(ind)<nrow(x)) {
          proj <- x[(max(ind)+1):nrow(x),]
        } else {
          proj <- NULL }  #pas de donn?e de projection

      } else {

        if (match(paste0("t__",t_init),rev(as.character(x$t)))==1) {
          proj <- NULL
        } else {
          proj <- x[(nrow(x)+2-match(paste0("t__",t_init),rev(as.character(x$t)))):nrow(x),] }

      }

      if (ncol(tab)==2) {

        if (is.null(proj))
          return(mean(as.numeric(as.character(tab$value)),na.rm=TRUE))
        else {
          intTab <- rbind.data.frame(tab[nrow(tab),],proj)
          intTab$value[1] <- mean(as.numeric(as.character(tab$value)),na.rm=TRUE)
          intTab$t <- gsub("t__t__","t__",intTab$t)
          return(expand.time(intTab,t_init,nbStep))
        }

      } else {

        nams <- names(tab) ; nams <- nams[-match(c("t","value"),nams)]
        eval(parse('',text=paste0("TAB <- with(tab,aggregate(as.numeric(as.character(value)),list(",
                                 paste(nams,collapse=","),"),mean,na.rm=TRUE))")))

        if (is.null(proj)) {
          names(TAB) <- c(nams,"value")
          return(TAB)
        } else {
          TAB$newField <- paste0("t__",t_init)
          names(TAB) <- c(nams,"value","t")
          return(expand.time(rbind.data.frame(TAB[,names(proj)],proj),t_init,nbStep))
        }
      }
    } else {
      return(x)
    }
  }# end of lapply
  })
  return(list(listHisto, listInput))
}


#' @importFrom openxlsx getSheetNames read.xlsx
#' @importFrom methods rbind2 new
#' @importFrom utils read.table
read.input <- function(file, t_init, nbStep, t_hist_max = t_init,
                       desc = "My input", folderFleet = NULL, verbose = FALSE ) {
  # if(!is.null(getOption("dev"))){ # will be triggered if option(dev = TRUE)
  #   rm(list = ls())
  #   library(openxlsx)
  #   file <- "inst/extdata/inputFile_simpl.xlsx"
  #   t_init <- 2020
  #   t_hist_max <- 2020
  #   nbStep <-  5
  #   desc <- "input tryhard"
  #   folderFleet = NULL
  # }
  ## Read sheets names ####
  tbls <- getSheetNames(file)
  nam_stock <- tbls[grep("stock",tolower(tbls))]
  n_stock <- length(nam_stock) # maybe add n_stock>0 variable here.
  if (is.null(folderFleet)) {
    namF <- tbls[grepl("^f__.*", tbls)]
  } else {
    namF <- sort(list.files(folderFleet))
    namF <- gsub(".csv","",namF[grepl("^f__.*", namF)])
  }

  namList <- gsub("Stock__","",nam_stock) # extract species names
  if(verbose) cat('Species :', namList, '\n')

  #LL <- list(historique=list(),input=list(),scenario=list()) ; LL$historique <- LL$input <- vector("list", length(namList))
  #names(LL$historique) <- names(LL$input) <- namList

  #modalites flottilles et metiers (bio et eco)
  modF <- NULL
  modMbio <- NULL
  modMeco <- NULL

  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")

  #yyy <- try(wb <- loadWorkbook(file))
  #if (sum(attributes(yyy)$class%in%"try-error")==1) {
  # index.openxlsx <- TRUE      #l'importation en utilisant XLConnect a ?chou?
  # print("'openxlsx' package will be used")
  #} else {
  # print("'XLConnect' package will be used")
  #}

  #if (index.openxlsx) {
  # detach("package:XLConnect", character.only = TRUE)
  # require(openxlsx)  #on travaille d?sormais avec openxlsx
  #}

  if (n_stock>0) {
    for (k in nam_stock) {

      result <- read.xlsx(file,sheet = k, rowNames = FALSE, colNames = FALSE,
                          skipEmptyRows = FALSE, skipEmptyCols = FALSE)
      result[] <- lapply(result, function(x) gsub(",",".",as.character(x)))
      result <- as.matrix(rbind2("",result))

      result[is.na(result)] <- ""

      #on commence par analyser les modalit?s de chaque type de variable
      vec <- as.vector(result)
      vec <- vec[vec!=""]
      #MOD <- lapply(c("f__","m__"),function(x) {gsub(x,"",unique(vec[grepl(x,vec)]))})    #plus lent
      MOD <- lapply(c("f__","m__"),function(x) {gsub(x,"",unique(vec[sapply(vec,function(y) substring(y,1,3)==x)]))})
      modF <- unique(c(modF,MOD[[1]])) ; modMbio <- unique(c(modMbio,MOD[[2]]))
    }
  }
  rm(k, vec, result, MOD)

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Read Fleet ####
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  Fleet <- NULL

  if (!is.null(folderFleet)) {
    for (k in 1:length(namF)) {
      if (k==1) {

        Fleet <- as.matrix(read.table(file.path(folderFleet,paste0(namF[k],".csv")),sep=";",quote="\""))
        Fleet <- Fleet[,-match("annee",Fleet[1,])][,1:7]
      } else {

        FLtemp <- as.matrix(read.table(file.path(folderFleet,paste0(namF[k],".csv")),sep=";",quote="\""))
        Fleet <- rbind2(Fleet,FLtemp[,-match("annee",FLtemp[1,])][-1,1:7])
      }
    }

  } else {
    for (k in 1:length(namF)) {
      if (k==1) {

        Fleet <- read.xlsx(file,sheet=namF[k],rowNames=FALSE,colNames=FALSE,skipEmptyRows = FALSE,skipEmptyCols = FALSE)[,1:7]
        Fleet[] <- lapply(Fleet, function(x) gsub(",",".",as.character(x)))
        Fleet <- as.matrix(Fleet)

        Fleet[is.na(Fleet)] <- ""
      } else {

        FLEETtmp <- read.xlsx(file,sheet=namF[k],rowNames=FALSE,colNames=FALSE,skipEmptyRows = FALSE,skipEmptyCols = FALSE)[-1,1:7]
        FLEETtmp[] <- lapply(FLEETtmp, function(x) gsub(",",".",as.character(x)))
        FLEETtmp <- as.matrix(FLEETtmp)

        FLEETtmp[is.na(FLEETtmp)] <- ""

        Fleet <- rbind2(Fleet,FLEETtmp)
      }
    }
    rm(FLEETtmp)
  }
  rm(k)

  Fleet[is.na(Fleet)] <- ""
  #on en profite pour finaliser modF et creer modMeco
  vec <- as.vector(Fleet)
  vec <- vec[vec!=""]
  MOD <- lapply(c("f__","m__"),function(x) {gsub(x,"",unique(vec[grepl(x,vec)]))})
  modF <- unique(c(modF,MOD[[1]])) ; modMeco <- unique(c(modMeco,MOD[[2]]))
  rm(vec)

  if (length(modMbio)==0) modMbio <- modMeco

  Fleet <- as.data.frame(Fleet[,c(1,4:7)])
  if (!is.null(folderFleet)) Fleet <- Fleet[-1,]
  Fleet[,5] <- suppressWarnings(as.numeric(as.character(Fleet[,5])))
  Fleet[,1] <- gsub("v__","",Fleet[,1])
  Fleet[,4] <- gsub("e__","",Fleet[,4])
  names(Fleet) <- c("v","f","m","e","value")

  #on distingue ce qui se decline par espece --> a integrer dans les parametres stocks
  Fstock <- Fleet[Fleet[,4]!="",]
  Fleet <- Fleet[Fleet[,4]=="",] ### End Fleet ####
  # ici ressort Fleet, Fstock, MOD, modF, modMbio, modMeco
  rm(folderFleet)
  if(verbose) cat('Fleet done \n')
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ## Init LL ####
  # require Fstock, namList
  nam_stock_bis <- unique(Fstock$e[!Fstock$e%in%c("espece",namList)])
  LL <- list(historique = list(), input = list(), scenario = list())
  LL$historique <- LL$input <- vector("list", length(namList)+length(nam_stock_bis))
  names(LL$historique) <- names(LL$input) <- c(namList,nam_stock_bis)


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##Scenarii ####
  ListS <- read.Scenar(file)
  if(verbose) cat('read.Scenar \n')
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ## Market sheet ####
  Market <- read.xlsx(file, sheet = "Market", rowNames = FALSE, colNames = FALSE,
                      skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  Market[] <- lapply(Market, function(x) gsub(",",".",as.character(x)))
  Market <- as.matrix(rbind2("",Market))

  Market[is.na(Market)] <- ""
  if(verbose) cat('Market read \n')
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Stock Param ####
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  iCATtab <- NULL

  if (n_stock>0) {
    ### Dyna Sp ####
    for (k in 1:n_stock) {
      if(verbose) cat('begin ', namList[k], ' -')

      result <- read.sheet(file = file, sheet = nam_stock[k])

      #on va ajouter la table market
      MarketSp <- Market[Market[,6]%in%paste0("e__",namList[k]),c(1,4:5,7:8),drop=FALSE]
      #on s?pare les tables par variables en intercalant une (ou 2) ligne vide
      tabicat <- do.call("rbind",lapply(c("v__OD_e","v__theta_e","v__Pst_e","v__P_fmce","v__Q_fmce","v__alpha_fmce","v__beta_fmce","v__gamma_fmce"),
                                        function(x) MarketSp[c(NA,grep(x,MarketSp[,1]),NA),]))
      tabicat[,5][tabicat[,5]==""] <- "-1"

      #on commence par analyser les modalit?s de chaque type de variable (age et taille)
      vec <- as.vector(result)
      vec <- vec[vec!=""]
      #MOD <- lapply(c("i__","l__","c__"),function(x) {gsub(x,"",unique(vec[grepl(x,vec)]))})
      MOD <- lapply(c("i__","l__"),function(x) {gsub(x,"",unique(vec[sapply(vec,function(y) substring(y,1,3)==x)]))})

      if (k==1) {            #on insere les variables 'fm' et 'icat' et marche

        FM <- read.sheet(file = file, sheet = "fm_matrix")
        MM <- read.sheet(file = file, sheet = "mm_matrix")
        ICAT <- read.sheet(file = file, sheet = "icat_matrix")

        #on transforme un peu les deux matrices pour qu'elles aient le meme nombre de colonnes
        ncolMax <- max(ncol(result),ncol(FM),ncol(ICAT),ncol(tabicat),ncol(MM))
        result <- rbind2(rbind2(rbind2(rbind2(eval(parse('',text=paste0("cbind(",paste(c("result",rep("\"\"",ncolMax-ncol(result))),collapse=","),")"))),
                                              eval(parse('',text=paste0("cbind(",paste(c("FM",rep("\"\"",ncolMax-ifelse(length(FM)>0,ncol(FM),0))),collapse=","),")")))),
                                       eval(parse('',text=paste0("cbind(",paste(c("MM",rep("\"\"",ncolMax-ifelse(length(MM)>0,ncol(MM),0))),collapse=","),")")))),
                                eval(parse('',text=paste0("cbind(",paste(c("ICAT",rep("\"\"",ncolMax-ifelse(length(ICAT)>0,ncol(ICAT),0))),collapse=","),")")))),
                         eval(parse('',text=paste0("cbind(",paste(c("tabicat",rep("\"\"",ncolMax-ifelse(length(tabicat)>0,ncol(tabicat),0))),collapse=","),")"))))
      } else {  #on insere seulement les variables marche

        ncolMax <- max(ncol(result),ncol(tabicat))
        result <- rbind2(eval(parse('',text=paste0("cbind(",paste(c("result",rep("\"\"",ncolMax-ncol(result))),collapse=","),")"))),
                         eval(parse('',text=paste0("cbind(",paste(c("tabicat",rep("\"\"",ncolMax-ncol(tabicat))),collapse=","),")"))))

      }
      if(verbose) cat(' matrix')
      result[is.na(result)] <- "NA"

      #on termine en analysant les modalit?s de la derni?re variable (cat?gorie)
      MOD[[3]] <- gsub("c__","",unique(MarketSp[,4][sapply(MarketSp[,4],function(y) substring(y,1,3)=="c__")]))


      #on commence par extraire le tableau de codage des variables
      indEmpt <- suppressWarnings(apply(result,1,function(x) min(unlist(sapply(prefix,grep,x)))))
      recode <- result[5:39,1:4] #recode <- result[5:(match(TRUE,is.finite(indEmpt))-1),1:4] # qu'es-ce que ca fait la en fait ?
      #on compl?te les recodages non sp?cifi?s
      recode[recode[,2]%in%c("","NA"),2] <- recode[recode[,2]%in%c("","NA"),1]
      recode[is.na(recode[,2]),2] <- recode[is.na(recode[,2]),1]
      #table de recodage des variables
      rec <- as.data.frame(recode[2:(match("",recode[,1])-1),])
      names(rec) <- recode[1,]

      # Regrouped this part in a independant function maxime 30/09/21
      List <- result_filtre(result = result, indEmpt = indEmpt)

      #on va l?g?rement retoucher la table mm pour injecter dans la colonne value la place de l'indice m?tier_eco correspondant
      indMM <- (1:length(List))[unlist(lapply(List,function(x) ("v__mm"%in%x$v)))]
      if (length(indMM)>0) MMmodif <- List[[indMM]]
      #for (i in indMM) {
      #  tempMM <- MMmodif <- List[[i]]                                         #modif MM  23/02/2012
      #  tempMM$value[!is.na(tempMM$value)] <- match(gsub("m__","",tempMM[,5]),modMeco)[!is.na(tempMM$value)]
      #  tempMM <- tempMM[,c(1:4,6)]
      #  List[[i]] <- tempMM[apply(tempMM,1,function(x) !any(is.na(x))),]
      #  }            #modif 17/04/2013

      if(verbose) cat(' insert_mat')
      #si k==1, on n'oublie pas d'extraire les donn?es "fm" et "mm" pour les ins?rer avec les donn?es flottille par esp?ce
      if (k==1) {
        TAB_FM <- do.call("rbind",lapply(List,function(x) if ("v__fm"%in%x$v) {x$e <- gsub("e__","",x$e) ; x$v <- gsub("v__","",x$v) ; return(x[x$v%in%"fm",names(Fstock)]) } else NULL))

        TAB_MM <- do.call("rbind",lapply(List,function(x) if ("v__mm"%in%x$v) {x$e <- gsub("e__","",x$e) ; x$v <- gsub("v__","",x$v) ; return(x[x$v%in%"mm",]) } else NULL))

        #List <- lapply(List,function(x) if ("v__fm"%in%x$v) return(NULL) else x)
        iCATtab <- do.call("rbind",lapply(List,function(x) if ("v__icat"%in%x$v) return(x[x$v%in%"v__icat",])  else NULL))


        LLL <- length(List)
        invisible(sapply(LLL:1,function(x) if ("v__icat"%in%List[[x]]$v) List[[x]] <<- NULL)) #on efface les ?l?ments de v_icat
        Fstock <- rbind(Fstock,TAB_FM)#,TAB_MM)
      }

      List <- c(List,list(iCATtab[iCATtab$e%in%paste0("e__",namList[k]),-2]))

      #on peut int?grer ici les tables de sc?narios

      #test pour savoir quelles tables int?grer dans la liste d?j? construite
      testS <- lapply(ListS,function(x) {
        tst <- FALSE # TODO : single if statement here
        if ((as.character(x$v[1])%in%paste0("v__",as.character(rec$Alias))) & !"e"%in%names(x)) {
          tst <- TRUE
        } else {
          if ("e"%in%names(x)) {
            if (namList[k]%in%gsub("e__","",as.character(x$e))) tst <- TRUE
          }
        }
        return(tst)
        })

      #test de donn?e flottille (? n'op?rer que lors de la premi?re it?ration)
      if (k==1) testF <- lapply(ListS,function(x) { (!as.character(x$v[1])%in%paste0("v__",as.character(rec$Alias))) &
          (!"e"%in%names(x)) })

      ListStemp <- lapply(ListS,function(x) {x$v <- paste0(x$v,x$s) ; return(x)})
      keep <- ListStemp[unlist(testS)]
      #on retire la colonne "esp?ce" apr?s avoir filtr? sur l'esp?ce
      keep <- lapply(keep,function(x) x[x$e%in%paste0("e__",namList[k]),])
      keep <- lapply(keep,function(x) if ("e"%in%names(x)) x[,-match("e",names(x))])
      if (k==1) keepF <-  ListStemp[unlist(testF)]
      if (k==1) {if (length(keepF)>0) {keepF <- lapply(keepF,function(x) {x$v <- paste0(x$v,"f__") ; return(x)})
      keep <- c(keep,keepF)}
      } #on balise les infos Fleet

      if (length(keep)>0) List <- c(List,keep)


      List <- lapply(List,function(x) split(x[,-match("v",names(x)),drop=FALSE],as.character(x[,match("v",names(x))])))
      namL <- gsub("v__","",unlist(lapply(List,names)))
      List <- unlist(List,recursive=FALSE,use.names = FALSE)
      names(List) <- namL

      #il faut regrouper les tables de m?me variable d?coup?es en plusieurs parties
      Nam <- unique(names(List))
      List <- lapply(Nam,function(x) {tt <- do.call("rbind",List[names(List)%in%x])
      rownames(tt) <- NULL
      return(tt)})
      names(List) <- Nam

      if(verbose) cat(' histo')
      #on en fait maintenant des objets standards accompagn?s de leur attribut 'DimCst' pour les inputs, et on laisse sous forme de DF pour l'historique
      #il faut consid?rer l'historique... (t<=t_init)
      res <- init_listHisto(List, t_init, t_hist_max, nbStep)
      listHisto <- res[[1]] ; listInput <- res[[2]]
      rm(res)

      #on recode les noms de variables conform?ment ? 'rec' (on ajoute les variables SS3 et Sex-based pilotables par le module scenario)
      SS3nam_N <- paste0("Ni0_S",1:4)
      SS3nam_F <- paste0("Ffmi_",as.vector(t(outer(paste0("S",1:4),paste0("M",1:4),paste0))))
      SS3nam_Flanwt <- paste0("FLWfmi_",as.vector(t(outer(paste0("S",1:4),paste0("M",1:4),paste0))))
      SS3nam_Fdiswt <- paste0("FDWfmi_",as.vector(t(outer(paste0("S",1:4),paste0("M",1:4),paste0))))

      Sexnam_F <- paste0("F_i_G",1:2)
      Sexnam_Nt0 <- paste0("N_it0_G",1:2)
      Sexnam_Ni0 <- paste0("N_i0t_G",1:2)
      Sexnam_mat <- paste0("mat_i_G",1:2)
      Sexnam_M <- paste0("M_i_G",1:2)
      Sexnam_wS <- paste0("wStock_i_G",1:2)
      Sexnam_wL <- paste0("wL_i_G",1:2)
      Sexnam_wD <- paste0("wD_i_G",1:2)
      Sexnam_C <- paste0("C_i_G",1:2)
      Sexnam_doth <- paste0("doth_i_G",1:2)
      Sexnam_d <- paste0("d_i_G",1:2)
      Sexnam_Fbar <- paste0("Fbar_G",1:2)
      # TODO : why twice variable here ?
      renam <- c(as.character(rec$Variable),as.character(rec$Variable),c("OD_e","theta_e","Pst_e"),SS3nam_N,SS3nam_F,SS3nam_Flanwt,SS3nam_Fdiswt,
                 Sexnam_F,Sexnam_Nt0,Sexnam_Ni0,Sexnam_mat,Sexnam_M,Sexnam_wS,Sexnam_wL,Sexnam_wD,Sexnam_C,Sexnam_doth,Sexnam_d,Sexnam_Fbar) ; # TODO : why ??
      # TODO : why not use renam simply ?
      names(renam) <- c(as.character(rec$Alias),as.character(rec$Variable),c("OD_e","theta_e","Pst_e"),SS3nam_N,SS3nam_F,SS3nam_Flanwt,SS3nam_Fdiswt,
                        Sexnam_F,Sexnam_Nt0,Sexnam_Ni0,Sexnam_mat,Sexnam_M,Sexnam_wS,Sexnam_wL,Sexnam_wD,Sexnam_C,Sexnam_doth,Sexnam_d,Sexnam_Fbar)
      # TODO : same as before...use renam
      renam <- renam[!duplicated(names(renam))]
      names(listHisto) <- renam[names(listHisto)] ; names(listInput) <- renam[names(listInput)] # TODO : ...used to rm element not in name
      #et on applique le multiplicateur ? chaque variable dans les deux listes
      rec <- rec[suppressWarnings(!is.na(as.numeric(as.character(rec$Multi)))),]
      invisible(sapply(1:nrow(rec),function(x) if (as.character(rec$Variable)[x]%in%names(listHisto))
        try(listHisto[[as.character(rec$Variable)[x]]]$value <<-
              listHisto[[as.character(rec$Variable)[x]]]$value*as.numeric(as.character(rec$Multi))[x],silent=TRUE)))

      #il ne reste plus qu'? ajouter ? listInput les param?tres par esp?ce issus des fichiers flottilles
      Fle <- Fstock[Fstock[,4]%in%namList[k],] ; n <- unique(Fle[,1])
      Fle <- lapply(n,function(x) {df <- as.data.frame(Fle[Fle[,1]%in%x,c(2,3,5)]); rownames(df) <- 1:nrow(df); return(df)})
      Fle <- lapply(1:length(n),function(x) {df <- Fle[[x]] ; if (all(is.na(df[,3]))) df[,1:2] <- "" ; return(df)})
      Fle <- lapply(Fle,function(x) x[,c(apply(x[,1:(ncol(x)-1)],2,function(y) !all(y=="")),TRUE)])   #on ejecte les colonnes vides
      #on g?re les constantes
      Fle <- lapply(Fle,function(x) if (is.null(dim(x))) return(x[1]) else return(x))
      names(Fle) <- n
      listInput <- c(listInput,Fle)

      ALK <- NULL
      #on commence par analyser la cl? taille-?ge si elle existe
      if ("alk"%in%names(listInput)) {

        if (!all(is.na(listInput$alk))) {

          if (all(c("l","i","value")%in%names(listInput$alk))) {

            ALK <- suppressWarnings(with(listInput$alk,tapply(as.numeric(as.character(value)),list(as.character(l),
                                                                                                   factor(as.character(i),levels=paste0("i__",MOD[[1]]))),function(x) x)))
            dimnames(ALK) <- lapply(dimnames(ALK), function(x) sapply(1:length(x),function(y) substring(x[y],4,nchar(x[y]))))

            listInput <- listInput[-match("alk",names(listInput))]

          }
        }
      }

      #on ajoute mm ? partir de TAB_MM
      tabMMtemp <- TAB_MM[TAB_MM$e%in%namList[k],c(2,4,5,6)]
      if (is.null(tabMMtemp) || nrow(tabMMtemp)==0) {
        tabMMtemp <- NULL
      } else {
        colnames(tabMMtemp) <- c("f","mBio","mEco","value")
        rownames(tabMMtemp) <- 1:nrow(tabMMtemp)
      }

      #mod_i <- unique(unlist(lapply(listInput,function(x) if (length(x)<2) return(NULL) else if ("i"%in%names(x)) return(as.character(x$i)) else return(NULL) )))

      listInputBio <- lapply(listInput[!names(listInput)%in%c("GVLref_f_m_e","Lref_f_m_e","dd1_f_m_e","dd2_f_m_e","P_fmce","Q_fmce","P_fme","Q_fme","mm")],
                             standFormat,nbStep,paste0("f__",modF),paste0("m__",modMbio),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),ALK)
      listInputEco <- lapply(listInput[c("GVLref_f_m_e","Lref_f_m_e","dd1_f_m_e","dd2_f_m_e","P_fmce","Q_fmce","P_fme","Q_fme")],
                             standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),ALK)
      listInput <- c(listInputBio,list(mm=tabMMtemp),listInputEco)

      invisible(sapply(1:nrow(rec),function(x) if (as.character(rec$Variable)[x]%in%names(listInput))
        try(listInput[[as.character(rec$Variable)[x]]] <<-
              listInput[[as.character(rec$Variable)[x]]]*as.numeric(as.character(rec$Multi))[x],silent=TRUE)))

      #on ajoute les occurences ?ges, tailles et cat?gories
      listInput$modI <- MOD[[1]] ; listInput$modL <- MOD[[2]] ; listInput$modC <- MOD[[3]] ; listInput$alk <- ALK

      listScenar <- List[grepl("s__",names(List))]
      listScenar <- lapply(listScenar,function(x) if ("t"%in%names(x)) expand.time(x,t_init,nbStep,TRUE) else return(x))

      indEc <- apply(do.call("rbind",lapply(c("nbv_f_m","cnb_f_m","nbds_f_m","effort1_f_m","effort2_f_m","Lref_f_m","Lref_f_m_e","Lref_f_m_e","Lref_f_m_e","GVLref_f_m",
                                              "GVLref_f_m_e","GVLref_f_m_e","GVLref_f_m_e","gc_f_m","nbh_f_m","nbtrip_f_m","fc_f_m","vf_f_m",
                                              "ovc_f_m","oilc_f_m","bc_f_m","foc_f_m","icec_f_m","cshr_f_m"),function(x) grepl(x,names(listScenar)))),2,any)
      if(verbose) cat(' eco')
      listScenarBio <- lapply(listScenar[!indEc],standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),ALK,NA)
      listScenarEco <- lapply(listScenar[indEc],standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),ALK,NA)
      listScenar <- c(listScenarBio,listScenarEco)

      #on ajoute l'attribut 'intervention'
      for (nn in names(listScenar)) {
        if (length(grep("__x__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(1)  #1 -> multiplication
        if (length(grep("__+__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(2)  #2 -> addition
        if (length(grep("__o__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(3)  #3 -> remplacement
        if (length(grep("__x__",nn))==0 & length(grep("__+__",nn))==0 & length(grep("__o__",nn))==0) attributes(listScenar[[nn]])$type <- as.integer(0)  #0 -> par d?faut (multiplication ??)
      }

      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__x__","",NN))
      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__+__","",NN))
      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__o__","",NN))

      #if (k==1) {disti <- grepl("f__",names(listScenar)) ; )
      namS <- names(listScenar) ; namV <- sapply(namS,function(x) strsplit(x,"s__")[[1]][1])
      namS <- sapply(namS,function(x) strsplit(x,"s__")[[1]][2])
      #### on ne proc?de au recodage que sur les variables Stock (except?e les variables internes)    <<----- variables internes ? mettre ? jour ici  <<<--------
      testSt <- grepl("f__",namS) | grepl("Foth_i",namV) | grepl("F_fmi",namV) | grepl("Froth_i",namV) | grepl("FDWToth_i",namV) | grepl("FRWToth_i",namV) | grepl("Ni0_",namV) | grepl("N_i0t",namV)

      if (length(listScenar)>0) {
        names(listScenar)[!testSt] <- paste(paste(renam[namV],namS,sep="s__"),namList[k],sep="e__")[!testSt]
        names(listScenar)[testSt] <- paste(names(listScenar)[testSt],namList[k],sep="e__")
      }

      # TODO : use a function here !
      #il ne reste plus qu'? transformer les captures en nombres en captures en poids   (on fera de m?me pour la partie historique)
      if (length(listInput$Y_mi)==0) {
        if (length(listInput$C_mi)!=0) {
          listInput$Y_mi <- listInput$C_mi
          indProd <- attributes(listInput$C_mi)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_mi <- listInput$Y_mi*rep(listInput$wL_i,each=prod(indProd))/1000
        }
      }

      if (length(listInput$Y_mi_G1)==0) {
        if (length(listInput$C_mi_G1)!=0) {
          listInput$Y_mi_G1 <- listInput$C_mi_G1
          indProd <- attributes(listInput$C_mi_G1)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_mi_G1 <- listInput$Y_mi_G1*rep(listInput$wL_i_G1,each=prod(indProd))/1000
        }
      }

      if (length(listInput$Y_mi_G2)==0) {
        if (length(listInput$C_mi_G2)!=0) {
          listInput$Y_mi_G2 <- listInput$C_mi_G2
          indProd <- attributes(listInput$C_mi_G2)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_mi_G2 <- listInput$Y_mi_G2*rep(listInput$wL_i_G2,each=prod(indProd))/1000
        }
      }

      if (length(listInput$Y_i)==0) {
        if (length(listInput$C_i)!=0) {
          listInput$Y_i <- listInput$C_i
          indProd <- attributes(listInput$C_i)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_i <- listInput$Y_i*rep(listInput$wL_i,each=prod(indProd))/1000
        }
      }
      if (length(listInput$Y_i_G1)==0) {
        if (length(listInput$C_i_G1)!=0) {
          listInput$Y_i_G1 <- listInput$C_i_G1
          indProd <- attributes(listInput$C_i_G1)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_i_G1 <- listInput$Y_i_G1*rep(listInput$wL_i_G1,each=prod(indProd))/1000
        }
      }

      if (length(listInput$Y_i_G2)==0) {
        if (length(listInput$C_i_G2)!=0) {
          listInput$Y_i_G2 <- listInput$C_i_G2
          indProd <- attributes(listInput$C_i_G2)$DimCst[1:2]
          indProd[indProd==0] <- 1
          listInput$Y_i_G2 <- listInput$Y_i_G2*rep(listInput$wL_i_G2,each=prod(indProd))/1000
        }
      }

      listInput[is.na(names(listInput))] <- NULL # TODO : why not the same with listHisto

      LL$historique[[k]] <- listHisto
      LL$input[[k]] <- listInput
      LL$scenario <- c(LL$scenario,listScenar)

      if(verbose)cat(namList[k], ' done\n')
    }

  } else {
    ### No dyna Sp ####
    #on rattrape le coup sur listScenar si pas d'espece dynamique

    #on filtre de ListS les sc?narios especes, ils seront pris en charge ensuite dans la boucle des esp?ces statiques
    ListStmp <- ListS[!unlist(lapply(ListS,function(x) "e"%in%names(x)))]
    List <- lapply(ListStmp,function(x) {x$v <- paste0(x$v,x$s) ; return(x)})
    List <- lapply(List,function(x) split(x[,-match("v",names(x)),drop=FALSE],as.character(x[,match("v",names(x))])))
    namL <- gsub("v__","",unlist(lapply(List,names)))
    List <- unlist(List,recursive=FALSE,use.names = FALSE)
    names(List) <- namL
    Nam <- unique(names(List))
    List <- lapply(Nam,function(x) {tt <- do.call("rbind",List[names(List)%in%x])
    rownames(tt) <- NULL
    return(tt)})
    names(List) <- Nam
    listScenar <- List[grepl("s__",names(List))]
    listScenar <- lapply(listScenar,function(x) if ("t"%in%names(x)) expand.time(x,t_init,nbStep,TRUE) else return(x))
    listScenar <- lapply(listScenar,standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),NA,NA,NA,NA)

    for (nn in names(listScenar)) {
      if (length(grep("__x__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(1)  #1 -> multiplication
      if (length(grep("__+__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(2)  #2 -> addition
      if (length(grep("__o__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(3)  #3 -> remplacement
      if (length(grep("__x__",nn))==0 & length(grep("__+__",nn))==0 & length(grep("__o__",nn))==0) attributes(listScenar[[nn]])$type <- as.integer(0)  #0 -> par d?faut (multiplication ??)
    }

    names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__x__","",NN))
    names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__+__","",NN))
    names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__o__","",NN))

    names(listScenar) <- paste0(names(listScenar),"e__Fleet")

    LL$scenario <- c(LL$scenario,listScenar)

  }

  ### ? ####

  if (length(nam_stock_bis)>0) {

    for (k in 1:length(nam_stock_bis)) {

      nam <- nam_stock_bis[k]

      #on va ajouter la table market
      MarketSp <- Market[Market[,6]%in%paste0("e__",nam),c(1,4:5,8),drop=FALSE]
      #on s?pare les tables par variables en intercalant une (ou 2) ligne vide
      tabicat <- do.call("rbind",lapply(c("v__OD_e","v__theta_e","v__Pst_e","v__P_fme","v__Q_fme","v__alpha_fme","v__beta_fme","v__gamma_fme"),
                                        function(x) MarketSp[c(NA,grep(x,MarketSp[,1]),NA),]))
      tabicat[,4][tabicat[,4]==""] <- "-1"

      result <- tabicat
      result[is.na(result)] <- "NA"

      #on commence par extraire le tableau de codage des variables
      indEmpt <- suppressWarnings(apply(result,1,function(x) min(unlist(sapply(prefix,grep,x)))))
      # Regrouped this part in a independant function
      List <- result_filtre(result = result, indEmpt = indEmpt)

      #on peut int?grer ici les tables de sc?narios

      #test pour savoir quelles tables int?grer dans la liste d?j? construite
      testS <- lapply(ListS,function(x) { tst <- FALSE
      if ("e"%in%names(x)) {
        if (nam%in%gsub("e__","",as.character(x$e))) tst <- TRUE
      }
      return(tst)})

      ListStemp <- lapply(ListS,function(x) {x$v <- paste0(x$v,x$s) ; return(x)})
      keep <- ListStemp[unlist(testS)]
      #on retire la colonne "esp?ce" apr?s avoir filtr? sur l'esp?ce
      keep <- lapply(keep,function(x) x[x$e%in%paste0("e__",nam),])
      keep <- lapply(keep,function(x) if ("e"%in%names(x)) x[,-match("e",names(x))])

      if (length(keep)>0) List <- c(List,keep)



      List <- lapply(List,function(x) split(x[,-match("v",names(x)),drop=FALSE],as.character(x[,match("v",names(x))])))
      namL <- gsub("v__","",unlist(lapply(List,names)))
      List <- unlist(List,recursive=FALSE,use.names = FALSE)
      names(List) <- namL

      #il faut regrouper les tables de m?me variable d?coup?es en plusieurs parties
      Nam <- unique(names(List))
      List <- lapply(Nam,function(x) {tt <- do.call("rbind",List[names(List)%in%x])
      rownames(tt) <- NULL
      return(tt)})
      names(List) <- Nam


      #on en fait maintenant des objets standards accompagn?s de leur attribut 'DimCst' pour les inputs, et on laisse sous forme de DF pour l'historique
      #il faut consid?rer l'historique... (t<=t_init)
      res <- init_listHisto(List, t_init, t_hist_max, nbStep)
      listHisto <- res[[1]] ; listInput <- res[[2]]
      rm(res)

      #il ne reste plus qu'? ajouter ? listInput les param?tres par esp?ce issus des fichiers flottilles
      Fle <- Fstock[Fstock[,4]%in%nam,] ; n <- unique(Fle[,1])
      Fle <- lapply(n,function(x) {df <- as.data.frame(Fle[Fle[,1]%in%x,c(2,3,5)]); rownames(df) <- 1:nrow(df); return(df)})
      Fle <- lapply(1:length(n),function(x) {df <- Fle[[x]] ; if (all(is.na(df[,3]))) df[,1:2] <- "" ; return(df)})
      Fle <- lapply(Fle,function(x) x[,c(apply(x[,1:(ncol(x)-1)],2,function(y) !all(y=="")),TRUE)])   #on ejecte les colonnes vides
      #on g?re les constantes
      Fle <- lapply(Fle,function(x) if (is.null(dim(x))) return(x[1]) else return(x))
      names(Fle) <- n
      listInput <- c(listInput,Fle)

      listInput <- lapply(listInput,standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),alk=NULL)

      listScenar <- List[grepl("s__",names(List))]
      listScenar <- lapply(listScenar,function(x) if ("t"%in%names(x)) expand.time(x,t_init,nbStep,TRUE) else return(x))

      listScenar <- lapply(listScenar,standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),paste0("i__",MOD[[1]]),paste0("c__",MOD[[3]]),alk=NULL,NA)

      #on ajoute l'attribut 'intervention'
      for (nn in names(listScenar)) {
        if (length(grep("__x__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(1)  #1 -> multiplication
        if (length(grep("__+__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(2)  #2 -> addition
        if (length(grep("__o__",nn))>0) attributes(listScenar[[nn]])$type <- as.integer(3)  #3 -> remplacement
        if (length(grep("__x__",nn))==0 & length(grep("__+__",nn))==0 & length(grep("__o__",nn))==0) attributes(listScenar[[nn]])$type <- as.integer(0)  #0 -> par d?faut (multiplication ??)
      }

      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__x__","",NN))
      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__+__","",NN))
      names(listScenar) <- sapply(names(listScenar),function(NN) gsub("__o__","",NN))

      #if (k==1) {disti <- grepl("f__",names(listScenar)) ; )
      namS <- names(listScenar) ; namV <- sapply(namS,function(x) strsplit(x,"s__")[[1]][1])
      namS <- sapply(namS,function(x) strsplit(x,"s__")[[1]][2])

      if (length(listScenar)>0) names(listScenar) <- paste(names(listScenar),nam,sep="e__")

      listInput[is.na(names(listInput))] <- NULL

      LL$historique[[n_stock+k]] <- listHisto
      LL$input[[n_stock+k]] <- listInput
      LL$scenario <- c(LL$scenario,listScenar)

    }
  }
  ### End Stock ####
  #il manque un element Fleet ? historique --> ? voir
  LL$historique$Fleet <- list() # TODO why here ?

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Fleet ####
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # require Fleet, nbStep, modF, modMeco and LL$input but only to create fleet element
  n <- unique(Fleet[,1])
  Fleet <- lapply(n,function(x) {df <- as.data.frame(Fleet[Fleet[,1]%in%x,c(2,3,5)]); rownames(df) <- 1:nrow(df); return(df)})
  Fleet <- lapply(1:length(n),function(x) {df <- Fleet[[x]] ; if (all(is.na(df[,3]))) df[,1:2] <- "" ; return(df)})
  Fleet <- lapply(Fleet,function(x) x[,c(apply(x[,1:(ncol(x)-1)],2,function(y) !all(y=="")),TRUE)])   # rm empty cols
  #on gere les constantes
  Fleet <- lapply(Fleet,function(x) if (is.null(dim(x))) return(x[1]) else return(x))
  names(Fleet) <- n
  LL$input$Fleet <- lapply(Fleet,standFormat,nbStep,paste0("f__",modF),paste0("m__",modMeco),"","",NULL) # reformat
  rm(Fleet, n)

  #on calcule les valeurs totales a partir des valeurs moyennes sur les champs "nbact_f_tot","nbds_f_tot","nbdf_f_tot",
  #"prodtot_f_tot","GR_f_tot","nbh_f_tot","nbds_f_m_tot","nbdf_f_m_tot","prodtot_f_m_tot","GR_f_m_tot","nbh_f_m_tot"
  namFtot <- c("effort_f_tot","nbds_f_tot","GVLref_f_tot","nbh_f_tot","effort_f_m_tot","nbds_f_m_tot","GVLref_f_m_tot","nbh_f_m_tot")
  Ftot <- LL$input$Fleet[gsub("_tot","",namFtot)]
  Ftot[[1]] <- LL$input$Fleet$effort1_f * LL$input$Fleet$effort2_f
  Ftot[[5]] <- LL$input$Fleet$effort1_f_m * LL$input$Fleet$effort2_f_m
  if (is.null(LL$input$Fleet$nbv_f)) LL$input$Fleet$nbv_f <- NA
  if (is.null(LL$input$Fleet$nbv_f_m)) LL$input$Fleet$nbv_f_m <- NA
  Ftot <- lapply(1:length(namFtot),function(x) if (x<5) Ftot[[x]]*as.vector(LL$input$Fleet$nbv_f) else Ftot[[x]]*as.vector(LL$input$Fleet$nbv_f_m))
  names(Ftot) <- namFtot
  LL$input$Fleet <- c(LL$input$Fleet,Ftot)
  l=lapply(LL$input[-length(LL$input)],function(x) x$Lref_f_m_e)
  l = lapply(l, function(x){ x[which(is.na(x))] = 0; return(x)})
  LL$input$Fleet$Yothsue_f_m <- (LL$input$Fleet$Lref_f_m - Reduce("+",l))/Ftot[[5]]
  LL$input$Fleet$Yothsue_f_m[which(is.na(LL$input$Fleet$Yothsue_f_m))] <- 0
  LL$input$Fleet$tripLgthIniMax_f_m <- LL$input$Fleet$tripLgth_f_m * as.vector(LL$input$Fleet$H_f*LL$input$Fleet$nbTrip_f*LL$input$Fleet$nbv_f / LL$input$Fleet$Lref_f)
  rm(Ftot,l, namFtot)
  # End with LL$input$Fleet
  ### END Fleet ####
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  # on va traiter les inputs sc?narios pour les organiser de la m?me mani?re que dans input
  SC <- LL$scenario
  # on recree une table des occurences des noms
  step1 <- t(sapply(names(SC),function(x) strsplit(x,"e__")[[1]])) ; rownames(step1) <- 1:nrow(step1)
  repF <- grepl("f__",step1[,1]) ; step1[,1] <- gsub("f__","",step1[,1]) ; step1[repF,2] <- "Fleet"
  step2 <- t(sapply(step1[,1],function(x) strsplit(x,"s__")[[1]])) ; rownames(step2) <- 1:nrow(step2)
  tabSc <- data.frame(v=step2[,1],s=step2[,2],e=step1[,2],ind=1:nrow(step2))

  gg <- lapply(split(tabSc,tabSc$s),function(x) split(x,x$e,drop=TRUE))
  gg <- lapply(gg,function(x) lapply(x,function(y) split(y,y$v,drop=TRUE)))
  LL$scenario <- lapply(gg,function(x) lapply(x,function(y) lapply(y,function(z) SC[[z$ind]])))

  ## Formating ####
  # il ne reste plus qu'a mettre en forme

  LL$historique <- c(lapply(LL$historique[1:n_stock],reformat),
                     if (length(nam_stock_bis)>0) lapply(LL$historique[(1:length(nam_stock_bis))+n_stock],reformat,"staticStockInput") else NULL,
                     lapply(LL$historique[n_stock+length(nam_stock_bis)+1],reformat,"fleetInput"))

  LL$input <- c(if (n_stock>0) lapply(LL$input[1:n_stock],reformat) else NULL,
                if (length(nam_stock_bis)>0) lapply(LL$input[(1:length(nam_stock_bis))+n_stock],reformat,"staticStockInput") else NULL,
                lapply(LL$input[n_stock+length(nam_stock_bis)+1],reformat,"fleetInput"))


  # on laisse l'element "scenario" tel quel pour le moment

  #on remplit la partie "stochastique" avec les variables issues des valeurs historiques de recrutement (seulement ?a pour le moment)
  STO <- list()

  stoch <- read.xlsx(file,sheet="Stochasticite_Sensibilite",rowNames=FALSE,colNames=FALSE,skipEmptyRows = FALSE,skipEmptyCols = FALSE)
  stoch[] <- lapply(stoch, function(x) gsub(",",".",as.character(x)))
  stoch <- as.matrix(stoch)

  stoch[is.na(stoch)] <- ""

  indexSt <- seq(from=grep("Samples : recruitment" ,stoch[,1]),to=grep("END Samples" ,stoch[,1])-1)
  indexSp <- grep("e__",stoch[,1])
  tabSto <- stoch[indexSp[indexSp%in%indexSt],]
  tabSto[,1] <- gsub("e__","",tabSto[,1])
  STO$RecHist <- lapply(as.character(namList),function(x) {tab <- tabSto[tabSto[,1]%in%x,,drop=FALSE] ;
  return(rep(as.numeric(tab[,3])*as.numeric(tab[,4]),as.numeric(tab[,5])))})

  names(STO$RecHist) <- as.character(namList)
  STO$GeoMeanRec <- lapply(STO$RecHist,function(x) if (length(x)>0) prod(x,na.rm=TRUE)^(1/sum(!is.na(x))) else numeric(0))
  STO$RecResiduals <- mapply(function(x,y) x-y,STO$RecHist,STO$GeoMeanRec,SIMPLIFY=FALSE)

  allSp <- tapply(tabSto[,1],list(tabSto[,2]),function(x) all(namList%in%x))       #ann?es pour lesquelles on a un historique pour toutes les eps?ces
  allSp <- names(allSp)[allSp]
  STO$RecHistLink <- lapply(as.character(namList),function(x) {tab <- tabSto[tabSto[,1]%in%x,,drop=FALSE] ; tab <- tab[match(allSp,tab[,2]),,drop=FALSE];
  return(rep(as.numeric(tab[,3])*as.numeric(tab[,4]),as.numeric(tab[,5])))})
  names(STO$RecHistLink) <- as.character(namList)
  STO$GeoMeanRecLink <- lapply(STO$RecHistLink,function(x) if (length(x)>0) prod(x,na.rm=TRUE)^(1/sum(!is.na(x))) else numeric(0))
  STO$RecResidualsLink <- mapply(function(x,y) x-y,STO$RecHistLink,STO$GeoMeanRecLink,SIMPLIFY=FALSE)


  #on remplit aussi, toujours pour le recrutement, la loi de distribution d?sir?e avec ses param?tres pour chaque esp?ce

  indexSt2 <- seq(from=grep("Random-variate : recruitment" ,stoch[,1]),to=grep("END Random-variate" ,stoch[,1])-1)
  tabSto2 <- stoch[indexSt2,]
  tabSto2 <- tabSto2[grep("e__",tabSto2[,5]),,drop=FALSE]
  tabSto2[,5] <- gsub("e__","",tabSto2[,5])

  ff <- lapply(1:4,function(x) {
    if (length(as.character(namList))>0) {
      comp <- lapply(as.character(namList),function(y) {tab <- tabSto2[tabSto2[,5]%in%y,,drop=FALSE] ;
      if (length(tab)>0) {
        if (tab[1,x]!="") {
          if (x%in%(2:4)) as.numeric(tab[1,x]) else as.character(tab[1,x])
        } else {
          if (x%in%(2:4)) as.numeric(NA) else ""
        }} else {
          if (x%in%(2:4)) as.numeric(NA) else ""}})
      names(comp) <- as.character(namList)
      return(comp)
    } else return(NA)
  })
  names(ff) <- c("RecDist","RecDistPar1","RecDistPar2","RecDistPar3")
  STO <- c(STO,ff)
  #on retourne alors l'objet bemInput renseign? (manquent encore les param?tres stochastiques et d'optimisation


  #modif MM 23/02/2012
  #? partir de MMmodif, on modifie mm pour chaque esp?ce
  indDYN <- (n_stock>0)
  if (indDYN) {
    if (exists("MMmodif")) {
      names(MMmodif)[4] <- "mBio"
      MMmodif$f <- gsub("f__","",MMmodif$f)
      MMmodif$e <- gsub("e__","",MMmodif$e)
      MMmodif$mBio <- gsub("m__","",MMmodif$mBio)
      MMmodif$m <- gsub("m__","",MMmodif$m)
      #tabMM <- with(MMmodif,tapply(value,list(paste(f,mBio,sep="__"),m,e),function(x) x[1]))
      #modif MM 27/08/2013 : ligne au dessus engendre un plantage si 1 seule occurence m?tier (g?n?ralisation du drop)
      tabMM <- lapply(unique(MMmodif$e),function(z) with(MMmodif[MMmodif$e%in%z,],tapply(value,list(paste(f,mBio,sep="__"),m),function(x) x[1])))
      names(tabMM) <- unique(MMmodif$e)
      for (spp in as.character(namList)) LL$input[[spp]]$mm <- tabMM[[spp]]
    }
  }

  if("Price_flexibility" %in% tbls){ ## Price flexibility  ####
    ListPflex <- read.Pflex(file = file, nam_stock = nam_stock,
                            nam_stock_bis = nam_stock_bis)
    LL$input[['Market']] <- reformat(ListPflex,"marketInput")
  }


  Qvec <- as.integer(rep(0,length(namList))) ; names(Qvec) <- as.character(namList)
  Svec <- as.integer(rep(0,length(namList))) ; names(Svec) <- as.character(namList)

  ## Return ####
  return(new("iamInput", desc=desc,
             specific=list(Species=if (indDYN) as.character(namList) else character(0),
                           StaticSpp=as.character(nam_stock_bis),
                           AllSpp = c(if (indDYN) as.character(namList) else character(0),as.character(nam_stock_bis)),
                           Fleet=modF, Metier=modMbio, MetierEco=modMeco,
                           Ages=if (indDYN) lapply(LL$input,function(x) x$modI)[namList] else list(),
                           Cat=if (indDYN) lapply(LL$input,function(x) x$modC)[namList] else list(),
                           t_init=t_init,
                           NbSteps=as.integer(nbStep),
                           times=as.integer(as.character(seq(t_init,by=1,length=nbStep))),
                           Q=if(indDYN) Qvec else integer(),#initialise Q
                           S=if(indDYN) Svec else integer()),#initialise S
             historical = LL$historique, input = LL$input,
             scenario = LL$scenario, stochastic = STO))

}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#M?thodes d'importation des param?tres d'input

# Generic IAM.input ####


# only accpet .xlsx files, will stop on xls files.
# maybe possibility to load both with regex

#' 'iamInput' objects creator
#'
#' @param fileIN .xls standard IAM input file
#' @param fileSPEC deprecated argument.
#' @param fileSCEN deprecated argument.
#' @param fileSTOCH deprecated argument.
#' @param ... Further arguments used when the input is a .xls. Described below.
#'
#' @docType methods
#' @name IAM.input-methods
#' @aliases IAM.input
#' @export
setGeneric("IAM.input", function(fileIN, fileSPEC, fileSCEN, fileSTOCH, ...){
  standardGeneric("IAM.input")
}
)

# why do we use a method here, it can be a single function that create a iamInput class object...

# a partir d'un fichier .xls

#' @rdname IAM.input-methods
#'
#' @param t_init Only for first signature. Initial year.
#' @param nbStep Only for first signature. Number of timesteps (including initial year).
#' @param t_hist_max Only for first signature. Last year considered for 'historical' slot.
#' @param desc Object descriptor (default value : "My Input").
#' @param folderFleet Folder containing fleets input sheets (Optionnal. Default value : NULL).
#' @param Fq_i List containing SS3 parameters arrays for quarterly stock dynamic simulation.
#' One element per considered species, all names must match with input species.
#' Fishing mortality per season, morph, and age.
#' (Optionnal. Default value : NULL. Example : list(hake=array(...,dim=c(4,4,nAge))).
#' @param Fq_fmi See above. Fishing mortality per season, morph, fleet, "economic" metiers and age
#' (Optionnal. Default value : NULL).
#' @param Ni0q See above. Recruits numbers per season (Optionnal. Default value : NULL).
#'
#'
setMethod("IAM.input", signature("character", "missing", "missing", "missing"),
                   function(fileIN, t_init, nbStep=20, t_hist_max=t_init, desc="My Input", folderFleet=NULL,
                                    Fq_i=NULL,iniFq_i=NULL,Fq_fmi=NULL,iniFq_fmi=NULL,
                                    FqLwt_i=NULL,iniFqLwt_i=NULL,FqLwt_fmi=NULL,iniFqLwt_fmi=NULL,
                                    FqDwt_i=NULL,iniFqDwt_i=NULL,FqDwt_fmi=NULL,iniFqDwt_fmi=NULL,
                                    Nt0s1q=NULL,Ni0q=NULL,iniNt0q=NULL,matwt=NULL,
                                    Fg_fmi=NULL,  dg_fmi=NULL, verbose = FALSE, ...) {
                                                                               #Fq_i matrice season*morph*age
                                                                               #Fq_fmi matrice season*morph*flottille*metier*age
                                                                               #Nt0s1q matrice morph*age  (Effectifs initiaux de projection (saison 1 de t=1) par cohorte)
                                                                               #iniNt0q matrice season*morph*age  (Effectifs initiaux (t=0))
                                                                               #Ni0q vecteur season (Recrutement par saison)
                                                                               #matwt poids moyen pond?r? pour calcul de la SSB (morph*age)


if (!substring(fileIN,nchar(fileIN)-4,nchar(fileIN))%in%".xlsx") stop("'fileIN' must be an .xlsx file!!")
# if (!sub("(.*)([[:punct:]]xlsx*)", "\\2",fileIN)%in% c(".xlsx", ".xls")) {
#   stop("'fileIN' must be an .xlsx or .xls file!!")
# }
#require(abind)
if(verbose) cat('Highway to hell \n')
out <- read.input(normalizePath(fileIN),t_init=t_init,nbStep=nbStep,t_hist_max=t_hist_max,desc=desc,folderFleet=folderFleet, verbose = verbose)
if(verbose) cat('read.input done \n')

## SS3 init ####
# TODO : Replace this with a more simpler common check.
# TODO : maybe group this in a class !
nmQ <- names(Fq_i) ; nmQ <- names(Fq_fmi)[names(Fq_fmi)%in%nmQ] ; nmQ <- names(Nt0s1q)[names(Nt0s1q)%in%nmQ]
nmQ <- names(Ni0q)[names(Ni0q)%in%nmQ] ; nmQ <- names(FqLwt_i)[names(FqLwt_i)%in%nmQ] ; nmQ <- names(FqLwt_fmi)[names(FqLwt_fmi)%in%nmQ]
nmQ <- names(FqDwt_i)[names(FqDwt_i)%in%nmQ] ; nmQ <- names(FqDwt_fmi)[names(FqDwt_fmi)%in%nmQ]
nmQ <- names(iniFq_i)[names(iniFq_i)%in%nmQ] ; nmQ <- names(iniFq_fmi)[names(iniFq_fmi)%in%nmQ]
nmQ <- names(iniFqLwt_i)[names(iniFqLwt_i)%in%nmQ] ; nmQ <- names(iniFqLwt_fmi)[names(iniFqLwt_fmi)%in%nmQ]
nmQ <- names(iniFqDwt_i)[names(iniFqDwt_i)%in%nmQ] ; nmQ <- names(iniFqDwt_fmi)[names(iniFqDwt_fmi)%in%nmQ]
nmQ <- names(iniNt0q)[names(iniNt0q)%in%nmQ] ; nmQ <- names(matwt)[names(matwt)%in%nmQ]
out@specific$Q[out@specific$Species%in%nmQ] <- as.integer(1)

## SEX init ####
if(!is.null(Fg_fmi) & !is.null(dg_fmi)){
  nmS <- names(Fg_fmi) ; nmS <- names(dg_fmi)[names(dg_fmi)%in%nmS]
  out@specific$S[out@specific$Species%in%nmS] <- as.integer(1)
} else if(is.null(Fg_fmi) & is.null(dg_fmi)){
  nmS <- NULL
} else {
  stop('To add sex to IAM, you need both Fg_fmi and dg_fmi.')
}
if(verbose) cat('Initiated SS3 and sex \n')

if (all(is.na(out@specific$Species))){
  OUT <- out
}else{
  OUT <- convertInput(out,Fq_fmi=Fq_fmi, Fg_fmi=Fg_fmi, verbose = verbose)
}
if(verbose) cat('convertInput done \n')

#sorting
ODpar_list <- unlist(lapply(OUT@input,function(x) x$OD_e))
if (all(ODpar_list%in%0)) OD <- 0 else OD <- min(ODpar_list,na.rm=TRUE)
OUT@input$Fleet$sorting <- as.integer(as.character(OD))

    # Fq_i <- lapply(OUT@specific$Species[1],function(x)
    #                array(0.1,dim=c(4,4,length(OUT@specific$Ages[[x]])),
    #                         dimnames=list(paste0("S",1:4),paste0("M",1:4),OUT@specific$Ages[[x]])))
    #
    # Fq_fmi <- lapply(OUT@specific$Species[1],function(x)
    #                array(0.1,dim=c(4,4,length(OUT@specific$Fleet),length(OUT@specific$MetierEco),length(OUT@specific$Ages[[x]])),
    #                         dimnames=list(paste0("S",1:4),paste0("M",1:4),OUT@specific$Fleet,
    #                                      OUT@specific$MetierEco,OUT@specific$Ages[[x]])))
    #
    # Ni0q <- lapply(OUT@specific$Species[1],function(x) array(10183000/4,dim=4,dimnames=list(paste0("S",1:4))))
    #
    # names(Fq_i) <- names(Fq_fmi) <- names(Ni0q) <- OUT@specific$Species[1]

#gestion des indicateurs trimestre
         # mise en commun : seules les esp?ces renseign?es pour les 3 variables sont conserv?es
FL <- OUT@specific$Fleet
ME <- OUT@specific$MetierEco

## edit ss3 ####
if (length(nmQ)>0) {
   for (i in nmQ) {
      lFq_i <- lFq_fmi <- lFothq_i <- lFqLwt_i <- lFqLwt_fmi <- lFothqLwt_i <- lFqDwt_i <- lFqDwt_fmi <- lFothqDwt_i <- lNt0s1q <- lNi0q <- list()
      liniFq_i <- liniFq_fmi <- liniFothq_i <- liniFqLwt_i <- liniFqLwt_fmi <- liniFothqLwt_i <- liniFqDwt_i <- liniFqDwt_fmi <- liniFothqDwt_i <- liniNt0q <- lmatwt <- list()
      AG <- OUT@specific$Ages[[i]]
      Nt0s1q[[i]][2:4,AG[1]] <- 0 ; Nt0s1q[[i]][1,AG[1]] <- Ni0q[[i]][1]      #pas d'individus d'?ge 0 et de morph>1 ? l'?tat initial
      for (season in 1:4) {
          temp <- as.double(Ni0q[[i]][season]) ; attributes(temp)$DimCst <- as.integer(c(0,0,0,0))
          lNi0q[[paste0("Ni0_S",season,"M",season)]] <- temp
        for (morph in 1:4) {
          if (season==1) {
             temp <- adrop(Nt0s1q[[i]][morph,AG,drop=FALSE],1) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
             lNt0s1q[[paste0("Nt0_S1M",morph)]] <- temp
             temp <- adrop(matwt[[i]][morph,AG,drop=FALSE],1) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
             lmatwt[[paste0("matwt_M",morph)]] <- temp
          }
          temp <- adrop(Fq_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFq_i[[paste0("Fi_S",season,"M",morph)]] <- temp
          temp <- adrop(Fq_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          lFq_fmi[[paste0("Ffmi_S",season,"M",morph)]] <- temp
          temp <- adrop(Fq_i[[i]][season,morph,AG,drop=FALSE]-apply(Fq_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFothq_i[[paste0("Fothi_S",season,"M",morph)]] <- temp

          temp <- adrop(iniFq_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFq_i[[paste0("iniFi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFq_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          liniFq_fmi[[paste0("iniFfmi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFq_i[[i]][season,morph,AG,drop=FALSE]-apply(iniFq_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFothq_i[[paste0("iniFothi_S",season,"M",morph)]] <- temp

          temp <- adrop(FqLwt_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFqLwt_i[[paste0("FLWi_S",season,"M",morph)]] <- temp
          temp <- adrop(FqLwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          lFqLwt_fmi[[paste0("FLWfmi_S",season,"M",morph)]] <- temp
          temp <- adrop(FqLwt_i[[i]][season,morph,AG,drop=FALSE]-apply(FqLwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFothqLwt_i[[paste0("FLWothi_S",season,"M",morph)]] <- temp

          temp <- adrop(iniFqLwt_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFqLwt_i[[paste0("iniFLWi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFqLwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          liniFqLwt_fmi[[paste0("iniFLWfmi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFqLwt_i[[i]][season,morph,AG,drop=FALSE]-apply(iniFqLwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFothqLwt_i[[paste0("iniFLWothi_S",season,"M",morph)]] <- temp

          temp <- adrop(FqDwt_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFqDwt_i[[paste0("FDWi_S",season,"M",morph)]] <- temp
          temp <- adrop(FqDwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          lFqDwt_fmi[[paste0("FDWfmi_S",season,"M",morph)]] <- temp
          temp <- adrop(FqDwt_i[[i]][season,morph,AG,drop=FALSE]-apply(FqDwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          lFothqDwt_i[[paste0("FDWothi_S",season,"M",morph)]] <- temp

          temp <- adrop(iniFqDwt_i[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFqDwt_i[[paste0("iniFDWi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFqDwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
          liniFqDwt_fmi[[paste0("iniFDWfmi_S",season,"M",morph)]] <- temp
          temp <- adrop(iniFqDwt_i[[i]][season,morph,AG,drop=FALSE]-apply(iniFqDwt_fmi[[i]][season,morph,FL,ME,AG,drop=FALSE],c(1,2,5),sum,na.rm=TRUE),1:2)
          temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniFothqDwt_i[[paste0("iniFDWothi_S",season,"M",morph)]] <- temp

          temp <- adrop(iniNt0q[[i]][season,morph,AG,drop=FALSE],1:2) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(0,0,length(AG),0))
          liniNt0q[[paste0("iniNt0q_S",season,"M",morph)]] <- temp
        }
          #OUT@specific$Q[OUT@specific$Species%in%i] <- as.integer(1)
      }
      OUT@input[[i]] <- c(OUT@input[[i]],liniNt0q,lNi0q,lNt0s1q,lmatwt,liniFq_i,lFq_i,liniFq_fmi,lFq_fmi,liniFothq_i,lFothq_i,
                     liniFqLwt_i,lFqLwt_i,liniFqLwt_fmi,lFqLwt_fmi,liniFothqLwt_i,lFothqLwt_i,liniFqDwt_i,lFqDwt_i,
                     liniFqDwt_fmi,lFqDwt_fmi,liniFothqDwt_i,lFothqDwt_i)
      OUT@input[[i]]$N_it0 <- as.double(NA) ; attributes(OUT@input[[i]]$N_it0)$DimCst <- as.integer(c(0,0,0,0))
      OUT@input[[i]]$N_i0t <- as.double(NA) ; attributes(OUT@input[[i]]$N_i0t)$DimCst <- as.integer(c(0,0,0,0))
      OUT@input[[i]]$F_i <- as.double(NA) ; attributes(OUT@input[[i]]$F_i)$DimCst <- as.integer(c(0,0,0,0))
      OUT@input[[i]]$F_fmi <- as.double(NA) ; attributes(OUT@input[[i]]$F_fmi)$DimCst <- as.integer(c(0,0,0,0))
    }
}
## edit sex ####
if (length(nmS)>0) {
  for (i in nmS) {
    lFg_fmi <- ldg_fmi <- list()
    AG <- OUT@specific$Ages[[i]]

    for (gender in as.character(1:2)) {
      temp <- adrop(Fg_fmi[[i]][gender,FL,ME,AG,drop=FALSE],1) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
      lFg_fmi[[paste0("F_fmi_G",gender)]] <- temp

      temp <- adrop(dg_fmi[[i]][gender,FL,ME,AG,drop=FALSE],1) ; temp[] <- as.double(temp) ; attributes(temp)$DimCst <- as.integer(c(length(FL),length(ME),length(AG),0))
      ldg_fmi[[paste0("d_i_G",gender)]] <- temp
    }
    OUT@input[[i]] <- c(OUT@input[[i]],lFg_fmi)
    OUT@input[[i]]$d_i_G1 <- ldg_fmi$d_i_G1
    OUT@input[[i]]$d_i_G2 <- ldg_fmi$d_i_G2
    OUT@input[[i]]$N_it0 <- as.double(NA) ; attributes(OUT@input[[i]]$N_it0)$DimCst <- as.integer(c(0,0,0,0))
    OUT@input[[i]]$N_i0t <- as.double(NA) ; attributes(OUT@input[[i]]$N_i0t)$DimCst <- as.integer(c(0,0,0,0))
    OUT@input[[i]]$F_i <- as.double(NA) ; attributes(OUT@input[[i]]$F_i)$DimCst <- as.integer(c(0,0,0,0))
    OUT@input[[i]]$F_fmi <- as.double(NA) ; attributes(OUT@input[[i]]$F_fmi)$DimCst <- as.integer(c(0,0,0,0))
  }
}

return(OUT)

})




# outdated method for txt files

  # ? partir de fichiers .txt
#' @importFrom methods new
#' @aliases IAM.input,character,character,missing,missing-method
#' @rdname IAM.input-methods
setMethod("IAM.input", signature("character", "character", "missing", "missing"),
                   function(fileIN, fileSPEC, desc="My Input", ...){

if (substring(fileIN,nchar(fileIN)-3,nchar(fileIN))!=".txt") stop("'fileIN' must be a .txt file!!")
if (substring(fileSPEC,nchar(fileSPEC)-3,nchar(fileSPEC))!=".txt") stop("'fileSPEC' must be a .txt file!!")

stop("This depend on deprecated C++ function 'Fun' ")
# specific <- suppressWarnings(.Call("Fun",normalizePath(fileSPEC),NULL))
# input <- suppressWarnings(.Call("Fun",normalizePath(fileIN),specific))

# out <- new("iamInput",desc=desc,specific=specific,historical=list(),input=input,scenario=list(),stochastic=list())
# return(convertInput(out))

})


#' @importFrom methods new
#' @aliases IAM.input,character,character,character,missing-method
#' @rdname IAM.input-methods
setMethod("IAM.input", signature("character", "character", "character", "missing"),
                   function(fileIN, fileSPEC, fileSCEN, desc="My Input", ...){

if (substring(fileIN,nchar(fileIN)-3,nchar(fileIN))!=".txt") stop("'fileIN' must be a .txt file!!")
if (substring(fileSPEC,nchar(fileSPEC)-3,nchar(fileSPEC))!=".txt") stop("'fileSPEC' must be a .txt file!!")
if (substring(fileSCEN,nchar(fileSCEN)-3,nchar(fileSCEN))!=".txt") stop("'fileSCEN' must be a .txt file!!")

stop("This depend on deprecated C++ function 'Fun' ")
# specific <- suppressWarnings(.Call("Fun",normalizePath(fileSPEC),NULL))
# input <- suppressWarnings(.Call("Fun",normalizePath(fileIN),specific))
# scenario <- suppressWarnings(.Call("Fun",normalizePath(fileSCEN),specific))

# out <- new("iamInput",desc=desc,specific=specific,historical=list(),input=input,scenario=scenario,stochastic=list())
# return(convertInput(out))

})


#' @importFrom methods new
#' @aliases IAM.input,character,character,missing,character-method
#' @rdname IAM.input-methods
setMethod("IAM.input", signature("character", "character", "missing", "character"),
                   function(fileIN, fileSPEC, fileSTOCH, desc="My Input", ...){

if (substring(fileIN,nchar(fileIN)-3,nchar(fileIN))!=".txt") stop("'fileIN' must be a .txt file!!")
if (substring(fileSPEC,nchar(fileSPEC)-3,nchar(fileSPEC))!=".txt") stop("'fileSPEC' must be a .txt file!!")
if (substring(fileSTOCH,nchar(fileSTOCH)-3,nchar(fileSTOCH))!=".txt") stop("'fileSTOCH' must be a .txt file!!")

stop("This depend on deprecated C++ function 'Fun' ")
# specific <- suppressWarnings(.Call("Fun",normalizePath(fileSPEC),NULL))
# input <- suppressWarnings(.Call("Fun",normalizePath(fileIN),specific))
# stochastic <- suppressWarnings(.Call("Fun",normalizePath(fileSTOCH),specific))

# out <- new("iamInput",desc=desc,specific=specific,historical=list(),input=input,scenario=list(),stochastic=stochastic)
# return(convertInput(out))

})


#' @importFrom methods new
#' @aliases IAM.input,character,character,character,character-method
#' @rdname IAM.input-methods
setMethod("IAM.input", signature("character", "character", "character", "character"),
                   function(fileIN, fileSPEC, fileSCEN, fileSTOCH, desc="My Input", ...){

if (substring(fileIN,nchar(fileIN)-3,nchar(fileIN))!=".txt") stop("'fileIN' must be a .txt file!!")
if (substring(fileSPEC,nchar(fileSPEC)-3,nchar(fileSPEC))!=".txt") stop("'fileSPEC' must be a .txt file!!")
if (substring(fileSCEN,nchar(fileSCEN)-3,nchar(fileSCEN))!=".txt") stop("'fileSCEN' must be a .txt file!!")
if (substring(fileSTOCH,nchar(fileSTOCH)-3,nchar(fileSTOCH))!=".txt") stop("'fileSTOCH' must be a .txt file!!")

stop("This depend on deprecated C++ function 'Fun' ")
# specific <- suppressWarnings(.Call("Fun",normalizePath(fileSPEC),NULL))
# input <- suppressWarnings(.Call("Fun",normalizePath(fileIN),specific))
# scenario <- suppressWarnings(.Call("Fun",normalizePath(fileSCEN),specific))
# stochastic <- suppressWarnings(.Call("Fun",normalizePath(fileSTOCH),specific))

# out <- new("iamInput",desc=desc,specific=specific,historical=list(),input=input,scenario=scenario,stochastic=stochastic)
# return(convertInput(out))

})



#::::::::::::::
# Examples ####
#::::::::::::::


#out <- IAM.input("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Inputs_SIAD_SEL_2.xls",t_init=2010,nbStep=21)
#
#out <- IAM.input("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/input.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/specific.txt")
#
#out <- IAM.input("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/input.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/specific.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/scenario.txt")
#
#out <- IAM.input("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/input.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/specific.txt",
#		fileSTOCH="C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/stochastic.txt")
#
#out <- IAM.input("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/input.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/specific.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/scenario.txt",
#		"C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expInput/stochastic.txt")
#
#str(out)
#
