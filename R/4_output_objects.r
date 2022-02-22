
# iamOuput ####
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Definition de l'objet Output global sans replicats et test de validite
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# val.iamOutput <- function(object){
#
# #? remplir avec les tests de validit? ? appliquer aux donn?es de sortie
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }

#' Class "iamOutput"
#'
#'
#' @slot desc Copy of the desc slot from \code{\link[IAM]{iamInput-class}}
#' @slot arguments Copy of the arguments slot from \code{\link[IAM]{iamArgs-class}}
#' @slot specific Copy of the specific slot from \code{\link[IAM]{iamInput-class}}
#' @slot outputSp # TODO
#' @slot output # TODO
#'
#' @examples
#' showClass("iamOutput")
#'
#'
#' @author Mathieu Merzereaud
#' @name iamOutput-class
#' @export
setClass("iamOutput",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de d?finition "esp?ce"
		output = "list"          # indicateur sans niveau de d?finition "esp?ce"
	),
	prototype(
		desc = "My output",
		arguments = list(Recruitment = list(),
                     Replicates = list(),
                     Scenario = list(),
                     Gestion = list(),
                     Eco = list()),
		specific = list(Species = character(),
                    StaticSpp=character(),
    	              Fleet = character(),
                  	Metier = character(),
                  	MetierEco = character(),
                    Ages = list(),
                   	Cat = list(),
                    t_init = double(),
                    NbSteps = integer(),
                   	times = integer(),
                    Q=integer()),
		outputSp = list(F = list(),               #mortalit? par p?che non corrig?e (-> captures)
                    Fr = list(),              #mortalit? par p?che corrig?e (-> morts)
                    Fothi = list(),           #mortalit? "autres flottilles, autres m?tiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalit? totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et m?tiers mod?lis?s
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et m?tiers mod?lis?s (composante ?ge)
                    Ytot = list(),            #captures totales en poids (composante ?ge)
                    D = list(),               #rejets totaux en poids pour les flottilles et m?tiers mod?lis?s (composante ?ge)
                    Li = list(),              #d?barquements totaux aux ?ges en poids pour les flottilles et m?tiers bio mod?lis?s  (composante ?ge)
                    Lc = list(),              #d?barquements totaux en poids par cat?gorie pour les flottilles et m?tiers bio mod?lis?s  (composante ?ge)
                    Ltot = list(),            #d?barquements totaux en poids (composante ?ge)
		                L_et = list(),            #d?barquements totaux en poids
		                L_pt = list(),            #debarquements produit marche en poids
                    P = list(),               #prix moyen par esp?ce et cat?gorie
                    GVL_f_m_e = list(),       #CA total par esp?ce, flottille et m?tier ?co
                     GVLcom_f_m_e = list(),
                     GVLst_f_m_e = list(),
                    statY = list(),           #captures par flottille, m?tier pour les esp?ces sans dynamique
                    statL = list(),           #d?barquements par flottille, m?tier pour les esp?ces sans dynamique
                    statD = list(),           #rejets par flottille, m?tier pour les esp?ces sans dynamique
                    statP = list(),           #Prix pour les esp?ces sans dynamique
                    statGVL_f_m = list(),     #CA total par esp?ce statique , flottille et m?tier ?co
                     statGVLcom_f_m = list(),
                     statGVLst_f_m = list(),
                    PQuot = list(),
		                TradedQ = list(),
                    F_S1M1= list(),F_S1M2= list(),F_S1M3= list(),F_S1M4= list(),
                    F_S2M1= list(),F_S2M2= list(),F_S2M3= list(),F_S2M4= list(),
                    F_S3M1= list(),F_S3M2= list(),F_S3M3= list(),F_S3M4= list(),
                    F_S4M1= list(),F_S4M2= list(),F_S4M3= list(),F_S4M4= list(),
                    Fr_S1M1= list(),Fr_S1M2= list(),Fr_S1M3= list(),Fr_S1M4= list(),
                    Fr_S2M1= list(),Fr_S2M2= list(),Fr_S2M3= list(),Fr_S2M4= list(),
                    Fr_S3M1= list(),Fr_S3M2= list(),Fr_S3M3= list(),Fr_S3M4= list(),
                    Fr_S4M1= list(),Fr_S4M2= list(),Fr_S4M3= list(),Fr_S4M4= list(),
                    Z_S1M1= list(),Z_S1M2= list(),Z_S1M3= list(),Z_S1M4= list(),
                    Z_S2M1= list(),Z_S2M2= list(),Z_S2M3= list(),Z_S2M4= list(),
                    Z_S3M1= list(),Z_S3M2= list(),Z_S3M3= list(),Z_S3M4= list(),
                    Z_S4M1= list(),Z_S4M2= list(),Z_S4M3= list(),Z_S4M4= list(),
                    N_S1M1= list(),N_S1M2= list(),N_S1M3= list(),N_S1M4= list(),
                    N_S2M1= list(),N_S2M2= list(),N_S2M3= list(),N_S2M4= list(),
                    N_S3M1= list(),N_S3M2= list(),N_S3M3= list(),N_S3M4= list(),
                    N_S4M1= list(),N_S4M2= list(),N_S4M3= list(),N_S4M4= list(),
		                F_G1 = list(),
		                F_G2 = list(),#mortalit? par p?che non corrig?e (-> captures)
		                Fr_G1 = list(),
		                Fr_G2 = list(),#mortalit? par p?che corrig?e (-> morts)
		                Fothi_G1 = list(),
		                Fothi_G2 = list(),#mortalit? "autres flottilles, autres m?tiers"
		                Z_G1 = list(),
		                Z_G2 = list(), #mortalit? totale (F+M)
		                N_G1 = list(),
		                N_G2 = list(),#effectifs totaux en nombre
                    DD_efmi= list(),
                    DD_efmc= list(),
                    LD_efmi= list(),
                    LD_efmc= list(),
                    statDD_efm= list(),
                    statLD_efm= list(),
                    statLDst_efm= list(),
                    statLDor_efm= list(),
                    oqD_ef= list(),
                    oqD_e= list(),
                    oqDstat_ef= list(),
                    TACtot = list(),
                    TACbyF = list(),
		                PQuot_conv = list(),
		                diffLQ_conv = list()),
    output = list(#typeGest = integer(),                #type de sc?nario de gestion appliqu?
                  nbv_f = numeric(),                   #Nb de navires par flottille
                  effort1_f = numeric(),               #1?re composante d'effort moyen par an et par flottille
                  effort2_f = numeric(),               #2?me composante d'effort moyen par an et par flottille
                  nbv_f_m = numeric(),                 #Nb de navires par flottille-m?tier
                  effort1_f_m = numeric(),             #1?re composante d'effort moyen par an et par flottille-m?tier
                  effort2_f_m = numeric(),             #2?me composante d'effort moyen par an et par flottille-m?tier
                  allocEff_f_m = numeric(),
                  #Lbio_f = numeric(),
                  GVLtot_f_m = numeric(),
                  GVLav_f_m = numeric(),
                  GVLtot_f = numeric(),
                  GVLav_f = numeric(),
                  #GVLoths_f = numeric(),
                  NGVLav_f_m = numeric(),
                  NGVLav_f = numeric(),
                  ET_f_m = numeric(),
                  cnb_f_m = numeric(),
                  cnb_f = numeric(),
                  #vcst_f_m = numeric(),
                  #vcst_f = numeric(),
                  rtbs_f_m = numeric(),
                  rtbs_f = numeric(),
                  rtbsAct_f = numeric(),
                  cshrT_f_m = numeric(),
                  cshrT_f = numeric(),
                  ncshr_f = numeric(),
                  ocl_f = numeric(),
                  cs_f = numeric(),
                  csAct_f = numeric(),
                  csTot_f = numeric(),
                  gva_f = numeric(),
                  gvaAct_f = numeric(),
                  gvamargin_f = numeric(),
                  gva_FTE_f = numeric(),
                  ccw_f = numeric(),
                  ccwCr_f = numeric(),
                  wageg_f = numeric(),
                  wagen_f = numeric(),
                  wageg_FTE_f = numeric(),
                  wageg_h_f = numeric(),
                  gp_f = numeric(),
                  gpAct_f = numeric(),
                  gpmargin_f = numeric(),
                  ncf_f = numeric(),
                  np_f = numeric(),
                  npmargin_f = numeric(),
                  prof_f = numeric(),
                  npmargin_trend_f = numeric(),
                  ssTot_f = numeric(),
                  ps_f = numeric(),
                  psAct_f = numeric(),
                  sts_f = numeric(),
                  stsAct_f = numeric(),
                  BER_f = numeric(),
                  CR_BER_f = numeric(),
                  fuelEff_f = numeric(),
                  ratio_fvol_gva_f = numeric(),
                  ratio_gp_gva_f = numeric(),
                  ratio_GVL_K_f = numeric(),
                  ratio_gp_K_f = numeric(),
                  RoFTA_f = numeric(),
                  ROI_f = numeric(),
                  ratio_np_K_f = numeric(),
                  ratio_GVL_cnb_ue_f = numeric(),
                  YTOT_fm = numeric(),
                  reconcilSPP = character(),
                  quotaExp_f = numeric(),
                  allocEff_f_m = numeric(),
                  GoFish = numeric())
  ) #,
	# validity=val.iamOutput
)



# iamOutputRep ####
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Definition de l'objet Output global avec replicats et test de validite
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# val.iamOutputRep <- function(object){
#
# #? remplir avec les tests de validit? ? appliquer aux donn?es de sortie
#   #if ... stop(...)
#
# 	return(TRUE)
#
# }


#' Class "iamOutputRep"
#'
#'
#' @slot desc Copy of the desc slot from \code{\link[IAM]{iamInput-class}}
#' @slot arguments Copy of the arguments slot from \code{\link[IAM]{iamArgs-class}}
#' @slot specific Copy of the specific slot from \code{\link[IAM]{iamInput-class}}
#' @slot outputSp # TODO
#' @slot output # TODO
#'
#' @examples
#' showClass("iamOutputRep")
#'
#'
#' @author Mathieu Merzereaud
#' @name iamOutputRep-class
#' @export
setClass("iamOutputRep",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de d?finition "esp?ce"
		output = "list"          # indicateur sans niveau de d?finition "esp?ce"
	),
	prototype(
		desc = "My output with replicates",
		arguments = list(Recruitment = list(),
                     Replicates = list(),
                     Scenario = list(),
                     Gestion = list(),
                     Eco = list()),
		specific = list(Species = character(),
                    StaticSpp=character(),
    	              Fleet = character(),
                  	Metier = character(),
                  	MetierEco = character(),
                    Ages = list(),
                   	Cat = list(),
                    t_init = double(),
                    NbSteps = integer(),
                   	times = integer(),
                    Q=integer()),
		outputSp = list(F = list(),               #mortalit? par p?che non corrig?e (-> captures)
                    Fr = list(),              #mortalit? par p?che corrig?e (-> morts)
                    Fothi = list(),           #mortalit? "autres flottilles, autres m?tiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalit? totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et m?tiers mod?lis?s
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et m?tiers mod?lis?s
                    Ytot = list(),            #captures totales en poids
                    D = list(),               #rejets totaux en poids pour les flottilles et m?tiers mod?lis?s
                    Li = list(),              #d?barquements totaux aux ?ges en poids pour les flottilles et m?tiers bio mod?lis?s
                    GVL_f_m_e = list(),       #CA total par esp?ce dynamique, flottille et m?tier ?co
                    statY = list(),           #captures par flottille, m?tier pour les esp?ces sans dynamique
                    statL = list(),           #d?barquements par flottille, m?tier pour les esp?ces sans dynamique
                    statD = list(),           #rejets par flottille, m?tier pour les esp?ces sans dynamique
                    statGVL_f_m = list(),     #CA total par esp?ce statique , flottille et m?tier ?co
                    PQuot = list(),
		                TradedQ = list()),
    output = list(nbv_f = list(),                   #Nb de navires par flottille
                  nbds_f = list(),                  #Nb de jdm moyen par an et par flottille
                  nbv_f_m = list(),                 #Nb de navires par flottille-m?tier
                  nbds_f_m = list(),                #Nb de jdm moyen par an et par flottille-m?tier
                  GVLtot_f_m = list(),              #CA total par flottille et m?tier
                  GVLtot_f = list(),                #CA total par flottille
                  GVLav_f = list(),                 #CA moyen par navire d'une flottille
                  vcst_f = list(),                  #Co?ts variables par navire d'une flottille
                  vcst_f_m = list(),                #Co?ts variables par navire d'une flottille-m?tier
                  rtbs_f = list(),                  #Reste ? partager par navire d'une flottille
                  rtbsAct_f = list(),               #Reste ? partager actualis? par navire d'une flottille
                  cs_f = list(),                    #Surplus de l'?quipage par navire d'une flottille
                  csAct_f = list(),                 #Surplus de l'?quipage actualis? par navire d'une flottille
                  gva_f = list(),                   #Valeur Ajout?e Brute pour un navire d'une flottille
                  gvaAct_f = list(),                #Valeur Ajout?e Brute actualis?e pour un navire d'une flottille
                  ccwCr_f = list(),                 #Co?ts du personnel par marin
                  wagen_f = list(),                 #Salaire net par marin
                  gcf_f = list(),                   #Exc?dent brut d'exploitation par navire d'une flottille
                  gcfAct_f = list(),                #Exc?dent brut d'exploitation actualis? par navire d'une flottille
                  gp_f = list(),                    #R?sultat net d'exploitation par navire d'une flottille
                  ps_f = list(),                    #Surplus total du producteur d'une flottille
                  psAct_f = list(),                 #Surplus total du producteur actualis? d'une flottille
                  sts_f = list(),                   #Surplus de l??tat associ? ? une flottille
                  stsAct_f = list())                #Surplus de l??tat actualis? associ? ? une flottille
  )# ,
	# validity=val.iamOutputRep
)


