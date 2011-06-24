
#====================================================================
# D�finition de l'objet global et test de validit� (param�tres requis renseign�s,...)
#====================================================================

val.iamInput <- function(object){

#� remplir avec les tests de validit� � appliquer aux donn�es de param�trage
  #if ... stop(...)
  
	return(TRUE)

}

setClass("iamInput",
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
		specific = list(Species=character(),
    	              Fleet=character(),
                  	Metier=character(),
                  	MetierEco=character(),	
                    Ages=list(),
                   	Cat=list(),
                    t_init=double(),
                    NbSteps=integer(),
                   	times=integer()),
		historical=list(),
		input=list(),
		scenario=list(),
		stochastic=list(),
		optimization=list()
  ),
	validity=val.iamInput
)


#====================================================================
# D�finition de l'objet StockInput et test de validit� (param�tres requis renseign�s,...)
#====================================================================

val.stockInput <- function(object){

#� remplir avec les tests de validit� � appliquer aux donn�es de param�trage
  #if ... stop(...)
  
	return(TRUE)

}

setClass("stockInput",
	representation(
		stock="character",
		input="list"
	),
	prototype(
		stock="My stock",
		input=list(
		  #modalit�s
		  modI = NA,
		  modL = NA,
		  modC = NA,
		  #BIO
		  icat	= NA,     #Cl� Cat�gories commerciales/�ges	
      alk	= NA,       #Cl� tailles-�ges (non temporelle pour le moment, et non spatialis�e)	
      fm = NA,        #proportion par m�tier des valeurs de ventilation pour chaque flottille (par esp�ce)
      mm = NA,        #matrice de red�finition du niveau m�tier entre le module bio/march� et le module �co (f_m.bio_e -> m.�co) 
      M_i	= NA,	      #Mortalit� naturelle
      mat_i = NA,     #Ogive de maturit�	
      wStock_i	= NA, #Poids individuels moyens dans le stock (kg)
      wL_i	= NA,     #Poids individuels moyens dans les d�barquements (kg)	
      wD_i	= NA,     #Poids individuels moyens dans les rejets (kg)	
      N_it0	= NA,     #Effectifs population aux �ges � l'instant initial 	
      N_i0t	= NA,     #Effectifs population � l'�ge 0	
      F_fmi	= NA,	    #Mortalit� par p�che aux �ges	
      B_i	= NA,       #Biomasse aux �ges (t)	
      Y_mi = NA,      #Capture totale par m�tier et par a/l en poids pour ventilation de la mortalit� par p�che	(t)
      C_mi = NA,      #Capture totale par m�tier et par a/l en nombres pour ventilation de la mortalit� par p�che	
      Y_i = NA,       #Capture totale par a/l en poids pour ventilation de la mortalit� par p�che	(t)
      C_i = NA,       #Capture totale par a/l en nombres pour ventilation de la mortalit� par p�che	
      d_i= NA,        #Proportion des captures totales rejet�es "flottilles mod�lis�es"
      doth_i= NA,     #Proportion des captures totales rejet�es "autres flottilles"
      sr = NA,        #Taux de survie des rejets
#      SelRef = NA,	  #Facteur de s�lectivit� de r�f�rence (PSo)
      P_fmce = NA,    #Prix moyen par cat�gorie (euros)	
      alpha_fmce = NA,     #Coefficient mod�le de prix	
      beta_fmce = NA,      #Coefficient mod�le de prix	
      gamma_fmce = NA,     #Coefficient mod�le de prix
      TAC = NA,
      Fbar = NA,
      FmaxTarget = NA,
      #ACTIVITE
      Lref_f_e = NA,  #Quantit� moyenne d�barqu�e par navire d'une flottille par an en tonnes par esp�ce
      Lref_f_m_e = NA,#Quantit� moyenne d�barqu�e par navire d'une flottille-m�tier par an en tonnes par esp�ce
      GVLref_f_e = NA,	  #Valeur moyenne d�barqu�e par navire d'une flottille par an en milliers d'euro par esp�ce
      GVLref_f_m_e = NA 	#Valeur moyenne d�barqu�e par navire d'une flottille-m�tier par an en milliers d'euro par esp�ce
    )
  ),
	validity=val.stockInput
)



#====================================================================
# D�finition de l'objet Fleet et test de validit� (param�tres requis renseign�s,...)
#====================================================================

val.fleetInput <- function(object){

#� remplir avec les tests de validit� � appliquer au donn�es de param�trage
  #if ... stop(...)
  
	return(TRUE)

}

setClass("fleetInput",
	representation(
		stock="character",
		input="list"
	),
	prototype(
		stock="My fleet",
		input=list(
		  #modalit�s
		  modF = NA,
		  modMbio = NA,
		  modMeco = NA,
		  #variables 
		  Lref_f = NA,             #D�barquements de r�f�rence par flottille
		  Lref_f_m = NA,           #D�barquements de r�f�rence par flottille - m�tier
      GVLref_f = NA,           #CA moyen initial par navire d'une flottille
      GVLref_f_m = NA,         #CA moyen initial par navire d'une flottille - m�tier
      nbv_f = NA,              #Nombre de navires par flottille
      nbv_f_m = NA,            #Nombre de navire par flottille - m�tier
      lc_f = NA,               #Taxes de d�barquement (% CA �co)
      gc_f = NA,               #Co�t total engins par navire d'une flottille
      gc_f_m = NA,             #Co�t total engins par navire d'une flottille - m�tier
      nbds_f = NA,             #Nombre moyen de Jours de Mer par navire d'une flottille par an
      nbds_f_m = NA,           #Nombre moyen de Jours de Mer par navire d'une flottille - m�tier par an
      nbh_f = NA,              #Nombre d'heures moteur par navire d'une flottille par an
      nbh_f_m = NA,            #Nombre d'heures moteur par navire d'une flottille - m�tier par an
      nbtrip_f = NA,           #Nombre de mar�es annuel par navire d'une flottille
      nbtrip_f_m = NA,         #Nombre de mar�es annuel par navire d'une flottile - m�tier
      fc_f = NA,               #Co�ts du carburant par navire d'une flottille
      fc_f_m = NA,             #Co�ts du carburant par navire d'une flottille - m�tier
      vf_f = NA,               #Prix du carburant par navire d'une flottille
      vf_f_m = NA,             #Prix du carburant par navire d'une flottille - m�tier
      ovc_f = NA,              #Autres co�ts variables par navire d'une flottille
      ovc_f_m = NA,            #Autres co�ts variables par navire d'une flottille - m�tier
      oilc_f = NA,             #Co�ts d'huile par navire d'une flottille
      oilc_f_m = NA,           #Co�ts d'huile par navire d'une flottille - m�tier
      bc_f = NA,               #Co�ts d'app�ts par navire d'une flottille
      bc_f_m = NA,             #Co�ts d'app�ts par navire d'une flottille - m�tier
      foc_f = NA,              #Co�ts de vivres par navire d'une flottille
      foc_f_m = NA,            #Co�ts de vivres par navire d'une flottille - m�tier
      cnb_f = NA,              #Effectif moyen par navire d'une flottille
      cnb_f_m = NA,            #Effectif moyen par navire d'une flottille - m�tier
      icec_f = NA,             #Co�ts de glace par navire d'une flottille
      icec_f_m = NA,           #Co�ts de glace par navire d'une flottille - m�tier
      cshr_f = NA,             #Part �quipage (ratio du RAP) par navire d'une flottille
      cshr_f_m = NA,           #Part �quipage (ratio du RAP) par navire d'une flottille - m�tier
      eec_f = NA,              #Cotisations salariales par navire d'une flottille
      mwhg_f = NA,             #Salaire brut horaire minimum national
      mwh_f = NA,              #Salaire net horaire minimum national
      altwh_f = NA,            #Salaire net horaire alternatif
      rep_f = NA,              #Co�ts entretien et r�paration
      onvc_f = NA,             #Autres Co�ts Non Variables
      insp_f = NA,             #Primes d'assurance
      ownc_f = NA,             #Autres d�penses d'armement
      mngc_f = NA,             #Cotisation Centre de gestion
      licc_f = NA,             #Co�t Total Licences
      comc_f = NA,             #Taxes Comit�s
      finc_f = NA,             #Co�t du capital
      dep_f = NA,              #Amortissement total
      ic_f = NA,               #Int�r�t
      K_f = NA,                #Valeur d'assurance
      vc_f = NA,               #Co�t total gr�ement
      persc_f = NA,            #Co�ts de personnel
      ecc_f = NA,              #Cotisations patronales
      pl_f = NA,               #Cong�s pay�s
      ovcDCF_f = NA,	         #Autres co�ts variables DCF par navire d'une flottille
      fixc_f = NA,             #Co�ts fixes (DCF)

      nbds_f_tot = NA,         #Nombre total de Jours de Mer par navire par an et par flottille
      GVLref_f_tot = NA,       #Valeur totale d�barqu�e par navire par an et par flottille en milliers d'euros (CA nav)
      nbh_f_tot = NA,          #Nombre total d'heures moteur par navire par an et par flottille
      nbds_f_m_tot = NA,       #Nombre total de jour de mer par navire par an et par flottille - m�tier
      GVLref_f_m_tot = NA,     #Valeur totale d�barqu�e par navire par an et par flottille - m�tier en milliers d'euro 			
      nbh_f_m_tot = NA         #Nombre total d'heure moteur par navire par an  et par flottille - m�tier   
    )
  ),
	validity=val.fleetInput
)



#====================================================================
# D�finition de l'objet Args contenant les arguments et test de validit� (param�tres requis renseign�s,...)
#====================================================================

val.iamArgs <- function(object){

arg <- object@arguments
spe <- object@specific

#Recruitments
if (length(arg$Recruitment)!=length(spe$Species)) stop("wrong 'Recruitment' argument in iamArgs object!!") 
if (!all(unlist(lapply(arg$Recruitment,length)==8))) stop("missing argument in a 'Recruitment' element of iamArgs object!!")
if (!all(unlist(lapply(arg$Recruitment,function(x) x$modSRactive)%in%(0:1)))) stop("wrong 'modSRactive' argument in iamArgs object!!")
if (!all(unlist(lapply(arg$Recruitment,function(x) x$typeMODsr)%in%c("Mean","Hockey-Stick","Beverton-Holt","Ricker","Shepherd")))) stop("wrong 'typeMODsr' argument in iamArgs object!!")
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

#Scenario
ind <- is.null(arg$Scenario$ALLscenario)
if (!arg$Scenario$active%in%(0:1)) stop("wrong 'Scen$active' argument in iamArgs object!!")
if (arg$Scenario$active==1 & ind) warning("no available scenario. 'Scen$active' put to 0!!")
if (!ind) {if (!is.character(arg$Scenario$ALLscenario)) stop("wrong 'ALLscenario' argument in iamArgs object!!")}                           
if (!ind) {if (!arg$Scenario$SELECTscen%in%(c(1:length(arg$Scenario$ALLscenario)))) stop("wrong 'SELECTscen' argument in iamArgs object!!")}                           

#Gestion
if (!arg$Gestion$active%in%(0:1)) stop("wrong 'Gest$active' argument in iamArgs object!!")
if (!arg$Gestion$control%in%c("Nb navires","Nb jdm")) stop("wrong 'control' argument in iamArgs object!!")
if (!arg$Gestion$target%in%c("TAC","Fbar","TAC->Fbar")) stop("wrong 'target' argument in iamArgs object!!")
if (!arg$Gestion$espece%in%spe$Species) stop("wrong 'espece' argument in iamArgs object!!")
if (!arg$Gestion$delay%in%(1:spe$NbSteps)) stop("wrong 'delay' argument in iamArgs object!!")
if (!arg$Gestion$upd%in%(1:2)) stop("wrong 'upd' argument in iamArgs object!!")
if (!is.numeric(arg$Gestion$sup)) stop("wrong 'sup' argument in iamArgs object!!")
if (!is.numeric(arg$Gestion$inf)) stop("wrong 'inf' argument in iamArgs object!!")
if (length(arg$Gestion$tac)!=length(spe$times)) stop("wrong 'tac' argument in iamArgs object!!")
if (length(arg$Gestion$fbar)!=length(spe$times)) stop("wrong 'fbar' argument in iamArgs object!!")
if (length(arg$Gestion$mf)!=length(spe$Fleet)) stop("wrong 'mf' argument in iamArgs object!!")

#Eco
if (!arg$Eco$active%in%(0:1)) stop("wrong 'Eco$active' argument in iamArgs object!!")
if (!arg$Eco$type%in%(1:2)) stop("wrong 'Eco$type' argument in iamArgs object!!")
if (!arg$Eco$adj%in%(1:2)) stop("wrong 'Eco$adj' argument in iamArgs object!!")
if (!arg$Eco$lev%in%(1:2)) stop("wrong 'Eco$lev' argument in iamArgs object!!")
if (!arg$Eco$ue_choice%in%(1:2)) stop("wrong 'Eco$ue_choice' argument in iamArgs object!!")
if (!arg$Eco$oths%in%(0:1)) stop("wrong 'Eco$oths' argument in iamArgs object!!")
if (!arg$Eco$othsFM%in%(0:1)) stop("wrong 'Eco$othsFM' argument in iamArgs object!!")
if (!arg$Eco$perscCalc%in%(0:2)) stop("wrong 'Eco$perscCalc' argument in iamArgs object!!")
if (!arg$Eco$report%in%(0:1)) stop("wrong 'Eco$report' argument in iamArgs object!!")
if (!is.numeric(arg$Eco$dr)) stop("wrong 'Eco$dr' argument in iamArgs object!!")

	return(TRUE)

}

setClass("iamArgs",
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
    specific = list(Species=character(),
    	              Fleet=character(),
                  	Metier=character(),
                  	MetierEco=character(),	
                    Ages=list(),
                   	Cat=list(),
                    t_init=double(),
                    NbSteps=integer(),
                   	times=integer())
  ),
	validity=val.iamArgs
)


