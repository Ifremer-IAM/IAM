
#==================================================
# Définition de l'objet Output global sans réplicats et test de validité 
#==================================================

val.iamOutput <- function(object){

#à remplir avec les tests de validité à appliquer aux données de sortie
  #if ... stop(...)
  
	return(TRUE)

}



setClass("iamOutput",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de définition "espèce"
		output = "list"          # indicateur sans niveau de définition "espèce"
	),
	prototype(
		desc = "My output",
		arguments = list(Recruitment = list(),
                     Replicates = list(), 
                     Scenario = list(), 
                     Gestion = list(), 
                     Eco = list()),
		specific = list(Species = character(),
    	              Fleet = character(),
                  	Metier = character(),
                  	MetierEco = character(),	
                    Ages = list(),
                   	Cat = list(),
                    t_init = double(),
                    NbSteps = integer(),
                   	times = integer()),
		outputSp = list(F = list(),               #mortalité par pêche non corrigée (-> captures)
                    Fr = list(),              #mortalité par pêche corrigée (-> morts)
                    Fothi = list(),           #mortalité "autres flottilles, autres métiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalité totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et métiers modélisés
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et métiers modélisés
                    Ytot = list(),            #captures totales en poids
                    D = list(),               #rejets totaux en poids pour les flottilles et métiers modélisés
                    Li = list(),              #débarquements totaux aux âges en poids pour les flottilles et métiers bio modélisés
                    Lc = list(),              #débarquements totaux en poids par catégorie pour les flottilles et métiers bio modélisés        
                    Lcm = list(),             #débarquements totaux en poids par catégorie pour les flottilles et métiers éco modélisés  
                    P = list(),               #prix moyen par espèce et catégorie
                    GVL_f_m_e = list()),      #CA total par espèce, flottille et métier éco
    output = list(nbv_f = numeric(),                   #Nb de navires par flottille
                  nbds_f = numeric(),                  #Nb de jdm moyen par an et par flottille
                  nbv_f_m = numeric(),                 #Nb de navires par flottille-métier
                  nbds_f_m = numeric(),                #Nb de jdm moyen par an et par flottille-métier
                  Lbio_f = numeric(),                  #Débarquements totaux par flottille
                  GVLtot_f_m = numeric(),              #CA total par flottille et métier
                  GVLav_f_m = numeric(),               #CA moyen par navire d'une flottille-métier
                  GVLtot_f = numeric(),                #CA total par flottille
                  GVLav_f = numeric(),                 #CA moyen par navire d'une flottille
                  GVLoths_f = numeric(),               #CA "autres espèces" par flottille
                  NGVLav_f_m = numeric(),              #CA net par navire d'une flottille-métier 
                  NGVLav_f = numeric(),                #CA net par navire d'une flottille
                  vcst_f_m = numeric(),                #Coûts variables par navire d'une flottille-métier
                  vcst_f = numeric(),                  #Coûts variables par navire d'une flottille
                  rtbs_f_m = numeric(),                #Reste à partager par navire d'une flottille-métier
                  rtbs_f = numeric(),                  #Reste à partager par navire d'une flottille
                  rtbsAct_f = numeric(),               #Reste à partager actualisé par navire d'une flottille
                  cshrT_f_m = numeric(),               #Part équipage par navire d'une flottille-métier
                  cshrT_f = numeric(),                 #Part équipage par navire d'une flottille
                  sshr_f_m = numeric(),                #Part armement par navire d'une flottille-métier
                  sshr_f = numeric(),                  #Part armement par navire d'une flottille
                  ncshr_f = numeric(),                 #Part équipage nette par navire d'une flottille
                  ocl_f = numeric(),                   #Coût d'opportunité du travail par navire d'une flottille
                  cs_f = numeric(),                    #Surplus de l'équipage par navire d'une flottille
                  csAct_f = numeric(),                 #Surplus de l'équipage actualisé par navire d'une flottille
                  csTot_f = numeric(),                 #Surplus total de l'équipage d'une flottille
                  gva_f = numeric(),                   #Valeur Ajoutée Brute pour un navire d'une flottille
                  gvaAct_f = numeric(),                #Valeur Ajoutée Brute actualisée pour un navire d'une flottille
                  ccw_f = numeric(),                   #Coûts de personnel
                  ccwCr_f = numeric(),                 #Coûts du personnel par marin
                  wageg_f = numeric(),                 #Salaire brut par marin
                  wagen_f = numeric(),                 #Salaire net par marin
                  gcf_f = numeric(),                   #Excédent brut d'exploitation par navire d'une flottille
                  gcfAct_f = numeric(),                #Excédent brut d'exploitation actualisé par navire d'une flottille
                  ngcf_f = numeric(),                  #Excédent net d'exploitation par navire d'une flottille
                  gp_f = numeric(),                    #Résultat net d'exploitation par navire d'une flottille
                  ssTot_f = numeric(),                 #Surplus total de l'armateur d'une flottille 
                  ps_f = numeric(),                    #Surplus total du producteur d'une flottille
                  psAct_f = numeric(),                 #Surplus total du producteur actualisé d'une flottille
                  sts_f = numeric(),                   #Surplus de l’État associé à une flottille 
                  stsAct_f = numeric(),                #Surplus de l’État actualisé associé à une flottille 
                  ber_f = numeric(),                   #GVLf tel que gpf = 0
                  ratio_gva_GVL_f = numeric(),         #Ratio valeur ajoutée brute / CA d'un navire moyen d'une flottille 
                  ratio_gcf_GVL_f = numeric(),         #Ratio excédent brut d'exploitation / CA d'un navire moyen d'une flottille
                  ratio_fc_GVL_f = numeric(),          #Ratio coût de carburant / CA d'un navire moyen d'une flottille
                  ratio_oilc_GVL_f = numeric(),        #Ratio coût d'huile  / CA d'un navire moyen d'une flottille
                  ratio_bc_GVL_f = numeric(),          #Ratio coût d'appâts  / CA d'un navire moyen d'une flottille
                  ratio_foc_GVL_f = numeric(),         #Ratio coût de vivres  / CA d'un navire moyen d'une flottille
                  ratio_icec_GVL_f = numeric(),        #Ratio coût de glace  / CA d'un navire moyen d'une flottille
                  ratio_gc_GVL_f = numeric(),          #Ratio coût engins / CA d'un navire moyen d'une flottille
                  ratio_vc_GVL_f = numeric(),          #Ratio coût gréement / CA d'un navire moyen d'une flottille
                  ratio_rep_GVL_f = numeric(),         #Ratio coût entretien-réparation / CA d'un navire moyen d'une flottille
                  ratio_mngc_GVL_f = numeric(),        #Ratio cotisations centre de gestion / CA d'un navire moyen d'une flottille
                  ratio_licc_GVL_f = numeric(),        #Ratio coût licence / CA d'un navire moyen d'une flottille
                  ratio_fvol_GVL_f = numeric(),        #Ratio volume de carburant  / CA d'un navire moyen d'une flottille
                  ratio_fvol_Lbio_f = numeric(),       #Ratio volume de carburant  / débarquements totaux  d'une flottille
                  ratio_fvol_gva_f = numeric(),        #Ratio volume de carburant  / valeur ajoutée brute d'un navire moyen d'une flottille
                  ratio_gcf_gva_f = numeric(),         #Ratio excédent brut d'exploitation  / valeur ajoutée brute d'un navire moyen d'une flottille 
                  ratio_K_cnb_f = numeric(),           #Ratio valeur d'assurance  / effectif moyen d'un navire moyen d'une flottille
                  ratio_GVL_K_f = numeric(),           #Ratio CA / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_gcf_K_f = numeric(),           #Ratio excédent brut d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_ngcf_K_f = numeric(),          #Ratio excédent net d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_gp_K_f = numeric(),            #Ratio résultat net d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_GVL_cnb_ue_f = numeric())      #Ratio CA / effectif moyen  / unité d'effort  d'un navire moyen d'une flottille 
  ),
	validity=val.iamOutput
)




#==================================================
# Définition de l'objet Output global avec réplicats et test de validité 
#==================================================

val.iamOutputRep <- function(object){

#à remplir avec les tests de validité à appliquer aux données de sortie
  #if ... stop(...)
  
	return(TRUE)

}



setClass("iamOutputRep",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de définition "espèce"
		output = "list"          # indicateur sans niveau de définition "espèce"
	),
	prototype(
		desc = "My output with replicates",
		arguments = list(Recruitment = list(),
                     Replicates = list(), 
                     Scenario = list(), 
                     Gestion = list(), 
                     Eco = list()),
		specific = list(Species = character(),
    	              Fleet = character(),
                  	Metier = character(),
                  	MetierEco = character(),	
                    Ages = list(),
                   	Cat = list(),
                    t_init = double(),
                    NbSteps = integer(),
                   	times = integer()),
		outputSp = list(F = list(),               #mortalité par pêche non corrigée (-> captures)
                    Fr = list(),              #mortalité par pêche corrigée (-> morts)
                    Fothi = list(),           #mortalité "autres flottilles, autres métiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalité totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et métiers modélisés
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et métiers modélisés
                    Ytot = list(),            #captures totales en poids
                    D = list(),               #rejets totaux en poids pour les flottilles et métiers modélisés
                    Li = list(),              #débarquements totaux aux âges en poids pour les flottilles et métiers bio modélisés
                    GVL_f_m_e = list()),      #CA total par espèce, flottille et métier éco
    output = list(nbv_f = list(),                   #Nb de navires par flottille
                  nbds_f = list(),                  #Nb de jdm moyen par an et par flottille
                  nbv_f_m = list(),                 #Nb de navires par flottille-métier
                  nbds_f_m = list(),                #Nb de jdm moyen par an et par flottille-métier
                  GVLtot_f_m = list(),              #CA total par flottille et métier
                  GVLtot_f = list(),                #CA total par flottille
                  GVLav_f = list(),                 #CA moyen par navire d'une flottille
                  vcst_f = list(),                  #Coûts variables par navire d'une flottille
                  vcst_f_m = list(),                #Coûts variables par navire d'une flottille-métier
                  rtbs_f = list(),                  #Reste à partager par navire d'une flottille
                  rtbsAct_f = list(),               #Reste à partager actualisé par navire d'une flottille
                  cs_f = list(),                    #Surplus de l'équipage par navire d'une flottille
                  csAct_f = list(),                 #Surplus de l'équipage actualisé par navire d'une flottille
                  gva_f = list(),                   #Valeur Ajoutée Brute pour un navire d'une flottille
                  gvaAct_f = list(),                #Valeur Ajoutée Brute actualisée pour un navire d'une flottille
                  ccwCr_f = list(),                 #Coûts du personnel par marin
                  wagen_f = list(),                 #Salaire net par marin
                  gcf_f = list(),                   #Excédent brut d'exploitation par navire d'une flottille
                  gcfAct_f = list(),                #Excédent brut d'exploitation actualisé par navire d'une flottille
                  gp_f = list(),                    #Résultat net d'exploitation par navire d'une flottille
                  ps_f = list(),                    #Surplus total du producteur d'une flottille
                  psAct_f = list(),                 #Surplus total du producteur actualisé d'une flottille
                  sts_f = list(),                   #Surplus de l’État associé à une flottille 
                  stsAct_f = list())                #Surplus de l’État actualisé associé à une flottille 
  ),
	validity=val.iamOutputRep
)


