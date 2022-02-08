# iamInput ####

# val.iamInput <- function(object){
#
# #? remplir avec les tests de validite ? appliquer aux donnees de parametrage
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }

#' Class "iamInput"
#'
#' # TODO
#'
#' @slot desc Short description of the object. chr.
#' @slot specific structural dimension
#' \describe{
#'   \item{Species}{Stock vector with 3 letter abreviation of species name and sotck number. chr vector}
#'   \item{StaticSpp}{Static species. These ones are landed but are not included in dynamic structure.
#'   This information is extracted from the Fleet sheet. chr vector}
#'   \item{AllSpp}{All the stocks names concatenated from Species and StaticSpp. chr vector}
#'   \item{Fleet}{Fleet types, defined thanks to vessel length, area and times. A vesserl can have multiple Metier. chr vector}
#'   \item{Metier}{Typology depending on fishing gear, vessel length, area and times of fisheries. chr vector}
#'   \item{MetiersEco}{Economic metier, depending on fishing gear. chr vector}
#'   \item{Ages}{Ages structure for each dynamic species. List of chr vector}
#'   \item{Cat}{Economic type for each dynamic species. List of chr vector}
#'   \item{t_init}{Initial year of the model. int}
#'   \item{NbSteps}{Number of years for the model to be used. Initial year is included as step 0. int vector}
#'   \item{times}{Years used in the model. Defined with t_init and NbSteps}
#'   \item{Q}{Status of species. 0 mean XSA dynamic 1 mean SS3 dynamic}
#'   \item{S}{Status of species. 0 mean XSA dynamic 1 mean SEX dynamic. A species can't be SS3 and SEX.}
#' }
#' @slot historical # TODO descr juste pour les graphes.pas utilise dans le modele C++
#' @slot input # TODO descr
#' @slot scenario List of scenarii with each containing one element per species
#' and one supplementary element for Fleets. # TODO explain why is there a Fleet element here ?
#' # TODO describe the format a species element
#' @slot stochastic # TODO description is not possible with med style input.
#' @slot optimization # TODO description is not possible with med style input.
#'
#' @details The \code{iamInput} class has methods defined for creation of
#' a new object, usage in \code{iamArgs} initiation and in \code{IAM.Model}.
#'
#' @examples
#' showClass("iamInput")
#'
#' @author Mathieu Merzereaud
#' @name iamInput-class
#' @rdname iamInput-class
#' @export
setClass("iamInput", ## Class ####
	representation(
		desc="character",
		specific = "list",
		historical="list",
		input="list",
		scenario="list",
		stochastic="list",
		optimization="list"
	),
	prototype(
		desc="iamInput",
		specific = list(Species=character(), # list of species names
                    StaticSpp=character(),
    	              Fleet=character(),
                  	Metier=character(),
                  	MetierEco=character(),
                    Ages=list(),
                   	Cat=list(),
                    t_init=double(),
                    NbSteps=integer(),
                   	times=integer(),
                    Q=integer(),
		                S =integer()),
		historical=list(),
		input=list(),
		scenario=list(),
		stochastic=list(),
		optimization=list()
  ) #,
	# validity=val.iamInput
)


# StockInput ####
# not used at all...why ?

# val.stockInput <- function(object){
#
# #? remplir avec les tests de validite ? appliquer aux donnees de parametrage
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }

#' stockInput Class
#'
#' # TODO descr
#'
#' @slot stock # TODO descr
#' @slot input # TODO descr
#'
#' @details Used by \code{reformat} function. # TODO find in what extent
#'
#' @name stockInput-class
#' @rdname stockInput-class
setClass("stockInput", ## Class ####
         representation(
           stock="character",
           input="list"
         ),
         prototype(
           stock="My stock",
           input=list(
             #modalites
             modI = NA,
             modL = NA,
             modC = NA,
             #BIO
             icat	= NA,     #Cl? categories commerciales/?ges
             alk	= NA,       #Cl? tailles-?ges (non temporelle pour le moment, et non spatialis?e)
             fm = NA,        #proportion par metier des valeurs de ventilation pour chaque flottille (par espece)
             mm = NA,        #matrice de red?finition du niveau metier entre le module bio/marche et le module ?co (f_m.bio_e -> m.?co)
             M_i	= NA,
             M_i_G1	= NA,
             M_i_G2	= NA,#Mortalite naturelle
             mat_i = NA,
             mat_i_G1 = NA,
             mat_i_G2 = NA,#Ogive de maturit?
             wStock_i	= NA,
             wStock_i_G1	= NA,
             wStock_i_G2	= NA,#Poids individuels moyens dans le stock (kg)
             wL_i	= NA,
             wL_i_G1	= NA,
             wL_i_G2	= NA,#Poids individuels moyens dans les debarquements (kg)
             wD_i	= NA,
             wD_i_G1	= NA,
             wD_i_G2	= NA,#Poids individuels moyens dans les rejets (kg)
             N_it0	= NA,
             N_it0_G1	= NA,
             N_it0_G2	= NA,#Effectifs population aux ?ges ? l'instant initial
             N_i0t	= NA,
             N_i0t_G1	= NA,
             N_i0t_G2	= NA,#Effectifs population ? l'?ge 0
             #Ni0_S1	= NA, Ni0_S2	= NA, Ni0_S3	= NA, Ni0_S4	= NA,   #effectifs population ? l'?ge 0 pour chaque morph (ie saison)
             F_i = NA,
             F_i_G1 = NA,
             F_i_G2 = NA,#Mortalite par p?che aux ?ges
             F_fmi	= NA,
             #F_fmi_G1	= NA,
             #F_fmi_G2	= NA,#Mortalite par p?che aux ?ges ventil?e
             #Mortalite par p?che aux ?ges pour chaque saison et chaque morph
             #Fi_S1M1 = NA, Fi_S1M2 = NA, Fi_S1M3 = NA, Fi_S1M4 = NA, Fi_S2M1 = NA, Fi_S2M2 = NA, Fi_S2M3 = NA, Fi_S2M4 = NA,
             #Fi_S3M1 = NA, Fi_S3M2 = NA, Fi_S3M3 = NA, Fi_S3M4 = NA, Fi_S4M1 = NA, Fi_S4M2 = NA, Fi_S4M3 = NA, Fi_S4M4 = NA,
             #Mortalite par p?che flottille-metier-?ges pour chaque saison et chaque morph
             #Ffmi_S1M1 = NA, Ffmi_S1M2 = NA, Ffmi_S1M3 = NA, Ffmi_S1M4 = NA, Ffmi_S2M1 = NA, Ffmi_S2M2 = NA, Ffmi_S2M3 = NA, Ffmi_S2M4 = NA,
             #Ffmi_S3M1 = NA, Ffmi_S3M2 = NA, Ffmi_S3M3 = NA, Ffmi_S3M4 = NA, Ffmi_S4M1 = NA, Ffmi_S4M2 = NA, Ffmi_S4M3 = NA, Ffmi_S4M4 = NA,
             #Mortalite par p?che "autres" aux ?ges pour chaque saison et chaque morph
             #Fothi_S1M1 = NA, Fothi_S1M2 = NA, Fothi_S1M3 = NA, Fothi_S1M4 = NA, Fothi_S2M1 = NA, Fothi_S2M2 = NA, Fothi_S2M3 = NA, Fothi_S2M4 = NA,
             #Fothi_S3M1 = NA, Fothi_S3M2 = NA, Fothi_S3M3 = NA, Fothi_S3M4 = NA, Fothi_S4M1 = NA, Fothi_S4M2 = NA, Fothi_S4M3 = NA, Fothi_S4M4 = NA,
             B_i	= NA,
             B_i_G1	= NA,
             B_i_G2	= NA, #Biomasse aux ?ges (t)
             Y_mi = NA,
             Y_mi_G1 = NA,
             Y_mi_G2 = NA,#Capture totale par metier et par a/l en poids pour ventilation de la Mortalite par p?che	(t)
             C_mi = NA,
             C_mi_G1 = NA,
             C_mi_G2 = NA,#Capture totale par metier et par a/l en nombres pour ventilation de la Mortalite par p?che
             Y_i = NA,
             Y_i_G1 = NA,
             Y_i_G2 = NA,#Capture totale par a/l en poids pour ventilation de la Mortalite par p?che	(t)
             C_i = NA,
             C_i_G1 = NA,
             C_i_G2 = NA,#Capture totale par a/l en nombres pour ventilation de la Mortalite par p?che
             d_i= NA,        #Proportion des captures totales rejetees "flottilles modelisees"
             d_i_G1=NA,
             d_i_G2=NA,
             doth_i= NA,     #Proportion des captures totales rejetees "autres flottilles"
             doth_i_G1= NA,
             doth_i_G2= NA,
             #d_fmi_G1 = NA,
             #d_fmi_G2 = NA,
             dd1_f_m_e = NA, #taux de rejets exemption en % de la capture totale de l'espece
             dd2_f_m_e = NA, #taux de rejets exemption en % de la capture totale
             sr = NA,        #Taux de survie des rejets
             #      SelRef = NA,	  #Facteur de s?lectivit? de reference (PSo)
             r = NA,         #SPiCT : Intrinsic growth rate : growth, recruitment, natural mortality
             K = NA,         #SPiCT : Carrying capacity (or equilibrium biomass or virgin stock biomass)
             n = NA,         #SPiCT : Parameter determining the shape of the production curve
             sigmaF = NA,    #SPiCT : Standard deviation of F
             sigmaB = NA,    #SPiCT : Standard deviation of B
             P_fmce = NA,    #Prix moyen par cat?gorie (euros)
             Pst_e = NA,     #Prix farine
             OD_e = NA,      #obligation de debarquement ? (oui(1)/non(0))
             theta_e = NA,   #multiplicateur de prix pour les rejets debarques (0<=...<=1)
             alpha_fmce = NA,     #Coefficient mod?le de prix
             beta_fmce = NA,      #Coefficient mod?le de prix
             gamma_fmce = NA,     #Coefficient mod?le de prix
             TAC = NA,
             Fbar = NA,
             Fbar_G1 = NA,
             Fbar_G2 = NA,
             FmaxTarget = NA,
             #ACTIVITE
             Lref_f_e = NA,  #Quantit? moyenne debarquee par navire d'une flottille par an en tonnes par espece
             Lref_f_m_e = NA,#Quantit? moyenne debarquee par navire d'une flottille-metier par an en tonnes par espece
             GVLref_f_e = NA,	  #Valeur moyenne debarquee par navire d'une flottille par an en milliers d'euro par espece
             GVLref_f_m_e = NA 	#Valeur moyenne debarquee par navire d'une flottille-metier par an en milliers d'euro par espece
           )
         ) #,
         # validity=val.stockInput
)


# staticStockInput ####

# val.staticStockInput <- function(object){
#
# #? remplir avec les tests de validite ? appliquer aux donnees de parametrage
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }


#' staticStockInput Class
#'
#' # TODO
#'
#' @slot stock # TODO descr
#' @slot input # TODO descr
#'
#' @details Used by \code{reformat} function. # TODO find in what extent
#'
#' @name staticStockInput-class
#' @rdname staticStockInput-class
setClass("staticStockInput", ## Class ####v
	representation(
		stock="character",
		input="list"
	),
	prototype(
		stock="My stock",
		input=list(
      LPUE_f_m_e=NA,  #debarquements moyens par unit? d'effort pour les especes non modelisees (t/nbds)
      d_f_m_e=NA,     #Proportion des captures totales rejetees pour les especes non modelisees sur les flottilles modelisees
      dd1_f_m_e = NA, #taux de rejets exemption en % de la capture totale de l'espece
      dd2_f_m_e = NA, #taux de rejets exemption en % de la capture totale
      dst_f_m_e = NA, #taux de rejets debarques sous-taille en % du tonnage de rejets debarques de l'espece
      P_fme = NA,     #Prix moyen especes non modelisees (euros)
      Pst_e = NA,     #Prix farine
      OD_e = NA,      #obligation de debarquement ? (oui(1)/non(0))
      theta_e = NA,   #multiplicateur de prix pour les rejets debarques (0<=...<=1)
      alpha_fme = NA,     #Coefficient mod?le de prix
      beta_fme = NA,      #Coefficient mod?le de prix
      gamma_fme = NA,
      Lref_f_e = NA,
      Lref_f_m_e = NA,
      GVLref_f_e = NA,
      GVLref_f_m_e = NA
    )
  )# ,
	# validity=val.staticStockInput
)


# fleetInput ####

# val.fleetInput <- function(object){
#
# #? remplir avec les tests de validite ? appliquer au donnees de parametrage
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }

#' fleetInput class
#'
#' # TODO
#'
#' @slot stock # TODO descr
#' @slot input # TODO descr
#'
#' @details Used by \code{reformat} function. # TODO find in what extent
#'
#' @name fleetInput-class
#' @rdname fleetInput-class
setClass("fleetInput", ## Class ####
	representation(
		stock="character",
		input="list"
	),
	prototype(
		stock="My fleet",
		input=list(
		  #modalit?s
		  modF = NA,
		  modMbio = NA,
		  modMeco = NA,
		  #variables
      sorting = NA,
		  Lref_f = NA,             #debarquements de reference par flottille
		  Lref_f_m = NA,           #debarquements de reference par flottille - metier
      GVLref_f = NA,           #CA moyen initial par navire d'une flottille
      GVLref_f_m = NA,         #CA moyen initial par navire d'une flottille - metier
      Yothsue_f_m = NA,
      nbv_f = NA,              #Nombre de navires par flottille
      nbv_f_m = NA,            #Nombre de navire par flottille - metier
      lc_f_m = NA,             #Taxes de debarquement (% CA ?co)
      lcd_f_m = NA,            #Taxes de debarquement relatives aux rejets debarques sous-taille (% CA ?co)
      gc_f = NA,               #cout total engins par navire d'une flottille
      gc_f_m = NA,             #cout total engins par navire d'une flottille - metier
      nbds_f = NA,             #Nombre moyen de Jours de Mer par navire d'une flottille par an
      nbds_f_m = NA,           #Nombre moyen de Jours de Mer par navire d'une flottille - metier par an
      nbh_f = NA,              #Nombre d'heures moteur par navire d'une flottille par an
      nbh_f_m = NA,            #Nombre d'heures moteur par navire d'une flottille - metier par an
      nbTrip_f = NA,           #Nombre de mar?es annuel par navire d'une flottille
      nbTrip_f_m = NA,         #Nombre de mar?es annuel par navire d'une flottile - metier
      tripLgth_f = NA,         #dur?e moyenne d'une mar?e par navire d'une flottille
      tripLgth_f_m = NA,       #dur?e moyenne d'une mar?e par navire d'une flottile - metier
      tripLgthIniMax_f_m = NA,
      effort1_f = NA,
      effort1_f_m = NA,
      effort2_f = NA,
      effort2_f_m = NA,
      effort1max_nbds_f = NA,
      effort1max_nbTrip_f = NA,
      effort1max_f = NA,
      H_f = NA,
      fc_f = NA,               #couts du carburant par navire d'une flottille
      fc_f_m = NA,             #couts du carburant par navire d'une flottille - metier
      vf_f = NA,               #Prix du carburant par navire d'une flottille
      vf_f_m = NA,             #Prix du carburant par navire d'une flottille - metier
      ovc_f = NA,              #Autres couts variables par navire d'une flottille
      ovc_f_m = NA,            #Autres couts variables par navire d'une flottille - metier
      oilc_f = NA,             #couts d'huile par navire d'une flottille
      oilc_f_m = NA,           #couts d'huile par navire d'une flottille - metier
      bc_f = NA,               #couts d'app?ts par navire d'une flottille
      bc_f_m = NA,             #couts d'app?ts par navire d'une flottille - metier
      foc_f = NA,              #couts de vivres par navire d'une flottille
      foc_f_m = NA,            #couts de vivres par navire d'une flottille - metier
      cnb_f = NA,              #Effectif moyen par navire d'une flottille
      cnb_f_m = NA,            #Effectif moyen par navire d'une flottille - metier
      icec_f = NA,             #Couts de glace par navire d'une flottille
      icec_f_m = NA,           #Couts de glace par navire d'une flottille - metier
      cshr_f = NA,             #Part ?quipage (ratio du RAP) par navire d'une flottille
      cshr_f_m = NA,           #Part ?quipage (ratio du RAP) par navire d'une flottille - metier
      eec_f = NA,              #Cotisations salariales par navire d'une flottille
      mwhg_f = NA,             #Salaire brut horaire minimum national
      mwh_f = NA,              #Salaire net horaire minimum national
      altwh_f = NA,            #Salaire net horaire alternatif
      rep_f = NA,              #couts entretien et r?paration
      onvc_f = NA,             #Autres couts Non Variables
      insp_f = NA,             #Primes d'assurance
      ownc_f = NA,             #Autres d?penses d'armement
      mngc_f = NA,             #Cotisation Centre de gestion
      licc_f = NA,             #cout Total Licences
      comc_f = NA,             #Taxes Comit?s
      finc_f = NA,             #cout du capital
      dep_f = NA,              #Amortissement total
      ic_f = NA,               #Int?r?t
      K_f = NA,                #Valeur d'assurance
      vc_f = NA,               #cout total gr?ement
      persc_f = NA,            #couts de personnel
      ecc_f = NA,              #Cotisations patronales
      pl_f = NA,               #Cong?s pay?s
      ovcDCF_f = NA,	         #Autres couts variables DCF par navire d'une flottille
      ovcDCF_f_m = NA,	       #Autres couts variables DCF par navire d'une flottille - metier
      fixc_f = NA,             #couts fixes (DCF)
      FTE_f = NA,
      FTE_f_m = NA,
      inv_f = NA,
      effort_f_tot=NA,
      nbds_f_tot = NA,         #Nombre total de Jours de Mer par navire par an et par flottille
      GVLref_f_tot = NA,       #Valeur totale debarquee par navire par an et par flottille en milliers d'euros (CA nav)
      nbh_f_tot = NA,          #Nombre total d'heures moteur par navire par an et par flottille
      effort_f_m_tot=NA,
      nbds_f_m_tot = NA,       #Nombre total de jour de mer par navire par an et par flottille - metier
      GVLref_f_m_tot = NA,     #Valeur totale debarquee par navire par an et par flottille - metier en milliers d'euro
      nbh_f_m_tot = NA         #Nombre total d'heure moteur par navire par an  et par flottille - metier
    )
  )# ,
	# validity=val.fleetInput
)


# marketInput ####

# val.marketInput <- function(object){
#
#   #? remplir avec les tests de validite ? appliquer au donnees de parametrage
#   #if ... stop(...)
#
#   return(TRUE)
#
# }

#' marketInput Class
#'
#' # TODO
#'
#' @slot input # TODO descr
#'
#' @details Used by \code{reformat} function. # TODO find in what extent
#'
#' @name marketInput-class
#' @rdname marketInput-class
setClass("marketInput", ## Class ####
         representation(
           input="list"
         ),
         prototype(
           input=list(
             #modalit?s
             modE = NA,
             modP = NA,
             #variables
             ep = NA,                #correspondance espece/espece marche
             beta_pp = NA            #cross price flexibilities
           )
         )# ,
         # validity=val.marketInput
)


# iamArgs ####

#' Checking validity of iamArgs Class
#'
#' @param object \code{\link[IAM]{iamArgs-class}} object
#'
#' @importFrom stats na.omit
val.iamArgs <- function(object){

  arg <- object@arguments
  spe <- object@specific

  #Recruitments
  if (length(arg$Recruitment)!=sum(spe$Q%in%0)) stop("wrong 'Recruitment' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,length)==9))) stop("missing argument in a 'Recruitment' element of iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) x$modSRactive)%in%(0:1)))) stop("wrong 'modSRactive' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) x$typeMODsr)%in%c("Mean","Hockey-Stick","Beverton-Holt","Ricker","Shepherd","Quadratic-HS","Smooth-HS")))) stop("wrong 'typeMODsr' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) is.numeric(x$parAmodSR))))) stop("wrong 'modSRactive' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) is.numeric(x$parAmodSR))))) stop("wrong 'parAmodSR' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) is.numeric(x$parBmodSR))))) stop("wrong 'parBmodSR' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) is.numeric(x$parCmodSR))))) stop("wrong 'parCmodSR' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) is.numeric(x$wnNOISEmodSR))))) stop("wrong 'wnNOISEmodSR' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) x$simuSTOCHactive)%in%(0:1)))) stop("wrong 'simuSTOCHactive' argument in iamArgs object!!")
  if (!all(unlist(lapply(arg$Recruitment,function(x) x$typeSIMUstoch)%in%(1:3)))) stop("wrong 'typeSIMUstoch' argument in iamArgs object!!")

  #Replicates
  if (!arg$Replicates$active%in%(0:1)) stop("wrong 'Rep$active' argument in iamArgs object!!")
  if (!is.numeric(arg$Replicates$nbIter)) stop("wrong 'nbIter' argument in iamArgs object!!")
  if (arg$Replicates$active==1 & length(arg$Replicates$SELECTvar)==0) stop("no specified output variables for replicates in iamArgs object!!")

  #Scenario
  ind <- is.null(arg$Scenario$ALLscenario)
  if (!arg$Scenario$active%in%(0:1)) stop("wrong 'Scen$active' argument in iamArgs object!!")
  if (arg$Scenario$active==1 & ind) warning("no available scenario. 'Scen$active' put to 0!!")
  if (!ind) {if (!is.character(arg$Scenario$ALLscenario)) stop("wrong 'ALLscenario' argument in iamArgs object!!")}
  if (!ind) {if (!arg$Scenario$SELECTscen%in%(c(1:length(arg$Scenario$ALLscenario)))) stop("wrong 'SELECTscen' argument in iamArgs object!!")}

  #Gestion
  if (!arg$Gestion$active%in%(0:1)) stop("wrong 'Gest$active' argument in iamArgs object!!")
  if (!arg$Gestion$control%in%c("Nb vessels","Nb trips")) stop("wrong 'control' argument in iamArgs object!!")
  if (!arg$Gestion$target%in%c("TAC","Fbar","TAC->Fbar","biomasse")) stop("wrong 'target' argument in iamArgs object!!")
  if (!arg$Gestion$espece%in%c(spe$Species,spe$StaticSpp)) stop("wrong 'espece' argument in iamArgs object!!")
  if (!arg$Gestion$typeG%in%(0:1)) stop("wrong 'level' argument in iamArgs object!!")
  if (!arg$Gestion$delay%in%(1:spe$NbSteps)) stop("wrong 'delay' argument in iamArgs object!!")
  if (!arg$Gestion$upd%in%(1:2)) stop("wrong 'upd' argument in iamArgs object!!")
  if (!is.numeric(arg$Gestion$sup)) stop("wrong 'sup' argument in iamArgs object!!")
  if (!is.numeric(arg$Gestion$inf)) stop("wrong 'inf' argument in iamArgs object!!")
  if (length(arg$Gestion$tac)!=length(spe$times)) stop("wrong 'tac' argument in iamArgs object!!")
  if (length(arg$Gestion$fbar)!=length(spe$times)) stop("wrong 'fbar' argument in iamArgs object!!")
  if (nrow(arg$Gestion$mfm)!=length(spe$Fleet)) stop("wrong 'mfm' argument in iamArgs object!!")
  if (ncol(arg$Gestion$mfm)!=length(spe$MetierEco)) stop("wrong 'mfm' argument in iamArgs object!!")
  if (length(arg$Gestion$othSpSup)!=length(c(na.omit(spe$Species)))+length(c(na.omit(spe$StaticSpp)))-1) stop("wrong 'othSpSup' argument in iamArgs object!!")
  if (nrow(arg$Gestion$effSup)!=length(spe$Fleet)) stop("wrong 'effSup' argument in iamArgs object!!")
  if (ncol(arg$Gestion$effSup)!=length(spe$times)) stop("wrong 'effSup' argument in iamArgs object!!")



  #Eco
  if (!arg$Eco$active%in%(0:1)) stop("wrong 'Eco$active' argument in iamArgs object!!")
  if (!arg$Eco$type%in%(1:2)) stop("wrong 'Eco$type' argument in iamArgs object!!")
  if (!arg$Eco$adj%in%(1:2)) stop("wrong 'Eco$adj' argument in iamArgs object!!")
  #if (!arg$Eco$lev%in%(1:2)) stop("wrong 'Eco$lev' argument in iamArgs object!!")
  if (!arg$Eco$ue_choice%in%(1:2)) stop("wrong 'Eco$ue_choice' argument in iamArgs object!!")
  if (!arg$Eco$oths%in%(0:1)) stop("wrong 'Eco$oths' argument in iamArgs object!!")
  if (!arg$Eco$othsFM%in%(0:1)) stop("wrong 'Eco$othsFM' argument in iamArgs object!!")
  if (!arg$Eco$perscCalc%in%(0:4)) stop("wrong 'Eco$perscCalc' argument in iamArgs object!!")
  if (!arg$Eco$report%in%(0:1)) stop("wrong 'Eco$report' argument in iamArgs object!!")
  if (!is.numeric(arg$Eco$dr)) stop("wrong 'Eco$dr' argument in iamArgs object!!")

  return(TRUE)

}

#' Class "iamArgs"
#'
#' # TODO
#'
#' @slot desc Copy of the desc slot from \code{\link[IAM]{iamInput-class}}
#' @slot arguments Arguments set in the GUI window.
#' \describe{
#'   \item{Recruitment}{Parameters for dynamic XSA species. Equation and parameters for recruitment.}
#'   \item{Replicates}{Deprecated. Parameters for replication. Module activation, number of replicates et output variables selected.}
#'   \item{Scenario}{Scenario selection. Module activation, list of scenario and selected scenario.}
#'   \item{Gestion}{Gestion module activation. Multiple parameters to select a target species,
#'   a gestion control (nbv, nbds), bounds of gestion, and if the gestion is by TAC, Fleet etc.
#'   This is in way of deprecation since TACbyF is directly implemented in IAM.model function.}
#'   \item{Eco}{Economic module, partially deprecated. Only perscCalc and dr (discount rate) are used.}
#' }
#' @slot specific Copy of the specific slot from \code{\link[IAM]{iamInput-class}}
#'
#' @details Used by \code{IAM.Args} method that use tcltk package for GUI.
#'
#' @examples
#' showClass("iamArgs")
#'
#' @author Mathieu Merzereaud
#' @name iamArgs-class
#' @rdname iamArgs-class
#' @export
setClass("iamArgs", ## Class ####
	representation(
		desc="character",
		arguments="list",
		specific="list"
	),
	prototype(
		desc="My model",
		arguments=list(
		  Recruitment = list(),
      Replicates = list(),
      Scenario = list(),
      Gestion = list(),
      Eco = list()),
		specific = list(
		  Species=character(),
		  StaticSpp=character(),
		  Fleet=character(),
		  Metier=character(),
		  MetierEco=character(),
		  Ages=list(),
		  Cat=list(),
		  t_init=double(),
		  NbSteps=integer(),
		  times=integer(),
		  Q=integer())
	),
	validity=val.iamArgs
)

