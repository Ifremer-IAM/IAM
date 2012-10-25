
#==================================================
# D�finition de l'objet Output global sans r�plicats et test de validit� 
#==================================================

val.iamOutput <- function(object){

#� remplir avec les tests de validit� � appliquer aux donn�es de sortie
  #if ... stop(...)
  
	return(TRUE)

}



setClass("iamOutput",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de d�finition "esp�ce"
		output = "list"          # indicateur sans niveau de d�finition "esp�ce"
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
		outputSp = list(F = list(),               #mortalit� par p�che non corrig�e (-> captures)
                    Fr = list(),              #mortalit� par p�che corrig�e (-> morts)
                    Fothi = list(),           #mortalit� "autres flottilles, autres m�tiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalit� totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et m�tiers mod�lis�s
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et m�tiers mod�lis�s
                    Ytot = list(),            #captures totales en poids
                    D = list(),               #rejets totaux en poids pour les flottilles et m�tiers mod�lis�s
                    Li = list(),              #d�barquements totaux aux �ges en poids pour les flottilles et m�tiers bio mod�lis�s
                    Lc = list(),              #d�barquements totaux en poids par cat�gorie pour les flottilles et m�tiers bio mod�lis�s        
                    Lcm = list(),             #d�barquements totaux en poids par cat�gorie pour les flottilles et m�tiers �co mod�lis�s  
                    P = list(),               #prix moyen par esp�ce et cat�gorie
                    GVL_f_m_e = list()),      #CA total par esp�ce, flottille et m�tier �co
    output = list(nbv_f = numeric(),                   #Nb de navires par flottille
                  nbds_f = numeric(),                  #Nb de jdm moyen par an et par flottille
                  nbv_f_m = numeric(),                 #Nb de navires par flottille-m�tier
                  nbds_f_m = numeric(),                #Nb de jdm moyen par an et par flottille-m�tier
                  Lbio_f = numeric(),                  #D�barquements totaux par flottille
                  GVLtot_f_m = numeric(),              #CA total par flottille et m�tier
                  GVLav_f_m = numeric(),               #CA moyen par navire d'une flottille-m�tier
                  GVLtot_f = numeric(),                #CA total par flottille
                  GVLav_f = numeric(),                 #CA moyen par navire d'une flottille
                  GVLoths_f = numeric(),               #CA "autres esp�ces" par flottille
                  NGVLav_f_m = numeric(),              #CA net par navire d'une flottille-m�tier 
                  NGVLav_f = numeric(),                #CA net par navire d'une flottille
                  vcst_f_m = numeric(),                #Co�ts variables par navire d'une flottille-m�tier
                  vcst_f = numeric(),                  #Co�ts variables par navire d'une flottille
                  rtbs_f_m = numeric(),                #Reste � partager par navire d'une flottille-m�tier
                  rtbs_f = numeric(),                  #Reste � partager par navire d'une flottille
                  rtbsAct_f = numeric(),               #Reste � partager actualis� par navire d'une flottille
                  cshrT_f_m = numeric(),               #Part �quipage par navire d'une flottille-m�tier
                  cshrT_f = numeric(),                 #Part �quipage par navire d'une flottille
                  sshr_f_m = numeric(),                #Part armement par navire d'une flottille-m�tier
                  sshr_f = numeric(),                  #Part armement par navire d'une flottille
                  ncshr_f = numeric(),                 #Part �quipage nette par navire d'une flottille
                  ocl_f = numeric(),                   #Co�t d'opportunit� du travail par navire d'une flottille
                  cs_f = numeric(),                    #Surplus de l'�quipage par navire d'une flottille
                  csAct_f = numeric(),                 #Surplus de l'�quipage actualis� par navire d'une flottille
                  csTot_f = numeric(),                 #Surplus total de l'�quipage d'une flottille
                  gva_f = numeric(),                   #Valeur Ajout�e Brute pour un navire d'une flottille
                  gvaAct_f = numeric(),                #Valeur Ajout�e Brute actualis�e pour un navire d'une flottille
                  ccw_f = numeric(),                   #Co�ts de personnel
                  ccwCr_f = numeric(),                 #Co�ts du personnel par marin
                  wageg_f = numeric(),                 #Salaire brut par marin
                  wagen_f = numeric(),                 #Salaire net par marin
                  gcf_f = numeric(),                   #Exc�dent brut d'exploitation par navire d'une flottille
                  gcfAct_f = numeric(),                #Exc�dent brut d'exploitation actualis� par navire d'une flottille
                  ngcf_f = numeric(),                  #Exc�dent net d'exploitation par navire d'une flottille
                  gp_f = numeric(),                    #R�sultat net d'exploitation par navire d'une flottille
                  ssTot_f = numeric(),                 #Surplus total de l'armateur d'une flottille 
                  ps_f = numeric(),                    #Surplus total du producteur d'une flottille
                  psAct_f = numeric(),                 #Surplus total du producteur actualis� d'une flottille
                  sts_f = numeric(),                   #Surplus de l��tat associ� � une flottille 
                  stsAct_f = numeric(),                #Surplus de l��tat actualis� associ� � une flottille 
                  ber_f = numeric(),                   #GVLf tel que gpf = 0
                  ratio_gva_GVL_f = numeric(),         #Ratio valeur ajout�e brute / CA d'un navire moyen d'une flottille 
                  ratio_gcf_GVL_f = numeric(),         #Ratio exc�dent brut d'exploitation / CA d'un navire moyen d'une flottille
                  ratio_fc_GVL_f = numeric(),          #Ratio co�t de carburant / CA d'un navire moyen d'une flottille
                  ratio_oilc_GVL_f = numeric(),        #Ratio co�t d'huile  / CA d'un navire moyen d'une flottille
                  ratio_bc_GVL_f = numeric(),          #Ratio co�t d'app�ts  / CA d'un navire moyen d'une flottille
                  ratio_foc_GVL_f = numeric(),         #Ratio co�t de vivres  / CA d'un navire moyen d'une flottille
                  ratio_icec_GVL_f = numeric(),        #Ratio co�t de glace  / CA d'un navire moyen d'une flottille
                  ratio_gc_GVL_f = numeric(),          #Ratio co�t engins / CA d'un navire moyen d'une flottille
                  ratio_vc_GVL_f = numeric(),          #Ratio co�t gr�ement / CA d'un navire moyen d'une flottille
                  ratio_rep_GVL_f = numeric(),         #Ratio co�t entretien-r�paration / CA d'un navire moyen d'une flottille
                  ratio_mngc_GVL_f = numeric(),        #Ratio cotisations centre de gestion / CA d'un navire moyen d'une flottille
                  ratio_licc_GVL_f = numeric(),        #Ratio co�t licence / CA d'un navire moyen d'une flottille
                  ratio_fvol_GVL_f = numeric(),        #Ratio volume de carburant  / CA d'un navire moyen d'une flottille
                  ratio_fvol_Lbio_f = numeric(),       #Ratio volume de carburant  / d�barquements totaux  d'une flottille
                  ratio_fvol_gva_f = numeric(),        #Ratio volume de carburant  / valeur ajout�e brute d'un navire moyen d'une flottille
                  ratio_gcf_gva_f = numeric(),         #Ratio exc�dent brut d'exploitation  / valeur ajout�e brute d'un navire moyen d'une flottille 
                  ratio_K_cnb_f = numeric(),           #Ratio valeur d'assurance  / effectif moyen d'un navire moyen d'une flottille
                  ratio_GVL_K_f = numeric(),           #Ratio CA / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_gcf_K_f = numeric(),           #Ratio exc�dent brut d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_ngcf_K_f = numeric(),          #Ratio exc�dent net d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_gp_K_f = numeric(),            #Ratio r�sultat net d'exploitation / valeur d'assurance d'un navire moyen d'une flottille
                  ratio_GVL_cnb_ue_f = numeric())      #Ratio CA / effectif moyen  / unit� d'effort  d'un navire moyen d'une flottille 
  ),
	validity=val.iamOutput
)




#==================================================
# D�finition de l'objet Output global avec r�plicats et test de validit� 
#==================================================

val.iamOutputRep <- function(object){

#� remplir avec les tests de validit� � appliquer aux donn�es de sortie
  #if ... stop(...)
  
	return(TRUE)

}



setClass("iamOutputRep",
	representation(
		desc = "character",
		arguments = "list",
		specific = "list",
		outputSp = "list",       # indicateur avec niveau de d�finition "esp�ce"
		output = "list"          # indicateur sans niveau de d�finition "esp�ce"
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
		outputSp = list(F = list(),               #mortalit� par p�che non corrig�e (-> captures)
                    Fr = list(),              #mortalit� par p�che corrig�e (-> morts)
                    Fothi = list(),           #mortalit� "autres flottilles, autres m�tiers"
                    Fbar = list(),            #indice Fbar
                    Z = list(),               #mortalit� totale (F+M)
                    N = list(),               #effectifs totaux en nombre
                    B = list(),               #biomasse
                    SSB = list(),             #biomasse reproductrice
                    C = list(),               #captures totales en nombre pour les flottilles et m�tiers mod�lis�s
                    Ctot = list(),            #captures totales en nombre
                    Y = list(),               #captures totales en poids pour les flottilles et m�tiers mod�lis�s
                    Ytot = list(),            #captures totales en poids
                    D = list(),               #rejets totaux en poids pour les flottilles et m�tiers mod�lis�s
                    Li = list(),              #d�barquements totaux aux �ges en poids pour les flottilles et m�tiers bio mod�lis�s
                    GVL_f_m_e = list()),      #CA total par esp�ce, flottille et m�tier �co
    output = list(nbv_f = list(),                   #Nb de navires par flottille
                  nbds_f = list(),                  #Nb de jdm moyen par an et par flottille
                  nbv_f_m = list(),                 #Nb de navires par flottille-m�tier
                  nbds_f_m = list(),                #Nb de jdm moyen par an et par flottille-m�tier
                  GVLtot_f_m = list(),              #CA total par flottille et m�tier
                  GVLtot_f = list(),                #CA total par flottille
                  GVLav_f = list(),                 #CA moyen par navire d'une flottille
                  vcst_f = list(),                  #Co�ts variables par navire d'une flottille
                  vcst_f_m = list(),                #Co�ts variables par navire d'une flottille-m�tier
                  rtbs_f = list(),                  #Reste � partager par navire d'une flottille
                  rtbsAct_f = list(),               #Reste � partager actualis� par navire d'une flottille
                  cs_f = list(),                    #Surplus de l'�quipage par navire d'une flottille
                  csAct_f = list(),                 #Surplus de l'�quipage actualis� par navire d'une flottille
                  gva_f = list(),                   #Valeur Ajout�e Brute pour un navire d'une flottille
                  gvaAct_f = list(),                #Valeur Ajout�e Brute actualis�e pour un navire d'une flottille
                  ccwCr_f = list(),                 #Co�ts du personnel par marin
                  wagen_f = list(),                 #Salaire net par marin
                  gcf_f = list(),                   #Exc�dent brut d'exploitation par navire d'une flottille
                  gcfAct_f = list(),                #Exc�dent brut d'exploitation actualis� par navire d'une flottille
                  gp_f = list(),                    #R�sultat net d'exploitation par navire d'une flottille
                  ps_f = list(),                    #Surplus total du producteur d'une flottille
                  psAct_f = list(),                 #Surplus total du producteur actualis� d'une flottille
                  sts_f = list(),                   #Surplus de l��tat associ� � une flottille 
                  stsAct_f = list())                #Surplus de l��tat actualis� associ� � une flottille 
  ),
	validity=val.iamOutputRep
)


