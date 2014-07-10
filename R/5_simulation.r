

#-------------------------------------------------------------------------------

# M�thode de simulation, avec en entr�e deux objets de classe 'iamArgs' et 'iamInput'

#-------------------------------------------------------------------------------


setGeneric("IAM.model", function(objArgs, objInput, ...){
	standardGeneric("IAM.model")
	}
)


setMethod("IAM.model", signature("iamArgs","iamInput"),function(objArgs, objInput, desc=as.character(NA), mOTH=1, TACbyF=as.double(NA), 
                  parBehav=list(active=as.integer(0),type=as.integer(3),FMT=NULL,MU=NULL,MUpos=as.integer(0),ALPHA=NULL),
                  parOptQuot=list(active=as.integer(0),pxQuIni=NA, pxQuMin=0, pxQuMax=NA, lambda=NA, ftol=0.0000001),
                  tacControl=list(tolVarTACinf=NA,tolVarTACsup=NA,corVarTACval=NA,corVarTACnby=2,Blim=NA,Bmax=NA,BlimTrigger=as.integer(0)),
                  stochPrice=list(), #liste d'�l�ments nomm�s par esp�ce consid�r�e, chaque �l�ment �tant une liste selon le sch�ma :
                                #list(type=NA (ou 1 ou 2,...), distr=c("norm",NA,NA,NA) (ou "exp" ou...), parA=c(0,NA,NA,NA), parB=c(1,NA,NA,NA), parC=c(NA,NA,NA,NA)) 
                  ...){
	
	
	
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

if (is.null(parOptQuot$active)) parOptQuot$active <- as.integer(0)
if (is.null(parOptQuot$pxQuIni)) parOptQuot$pxQuIni <- NA
if (is.null(parOptQuot$pxQuMin)) parOptQuot$pxQuMin <- 0
if (is.null(parOptQuot$pxQuMax)) parOptQuot$pxQuMax <- NA
if (is.null(parOptQuot$lambda)) parOptQuot$lambda <- NA
if (is.null(parOptQuot$ftol)) parOptQuot$ftol <- 0.0000001

if (is.null(parBehav$active)) parBehav$active <- as.integer(0)
if (is.null(parBehav$type)) parBehav$type <- as.integer(3)
if (is.null(parBehav$MUpos)) parBehav$MUpos <- as.integer(0)

Rectyp <- unlist(lapply(objArgs@arguments$Recruitment,function(x) x$simuSTOCHactive * x$typeSIMUstoch))

mOth <- rep(mOTH,length=length(objArgs@specific$Species)) # ; mOth[match(objArgs@arguments$Gestion$espece,objArgs@specific$Species)] <- mOTH

nT <- objInput@specific$NbSteps

out <-  .Call("IAM", objInput@input, objInput@specific, objInput@stochastic, objInput@scenario[[scenar]],
                    RecType1=as.integer(Rectyp==1), RecType2=as.integer(Rectyp==2), RecType3=as.integer(Rectyp==3),
                    as.integer(objArgs@arguments$Scenario$active), as.integer(objArgs@arguments$Replicates$active),
                    as.integer(objArgs@arguments$Replicates$nbIter), as.integer(objArgs@arguments$Gestion$active),
                    as.double(objArgs@arguments$Gestion$mf),as.double(mOth),
                    as.double(c(objArgs@arguments$Gestion$inf,objArgs@arguments$Gestion$sup)),
                    as.double(objArgs@arguments$Gestion$tac),as.double(objArgs@arguments$Gestion$fbar),
                    as.integer(c(eTemp = match(objArgs@arguments$Gestion$espece,objArgs@specific$Species)-1,
                                 var = match(objArgs@arguments$Gestion$control,c("Nb jdm","Nb navires")),
                                 trgt = match(objArgs@arguments$Gestion$target,c("TAC","Fbar","TAC->Fbar")),
                                 delay = objArgs@arguments$Gestion$delay,
                                 upd = objArgs@arguments$Gestion$upd, level = objArgs@arguments$Gestion$level)),
                    as.integer(objArgs@arguments$Eco$type-1),
                    as.integer(c(adj = objArgs@arguments$Eco$adj,
                                 lev = objArgs@arguments$Eco$lev,
                                 ue_choice = objArgs@arguments$Eco$ue_choice,
                                 oths = objArgs@arguments$Eco$oths,
                                 othsFM = objArgs@arguments$Eco$othsFM,
                                 perscCalc = objArgs@arguments$Eco$perscCalc,
                                 report = objArgs@arguments$Eco$report)),
                    as.double(objArgs@arguments$Eco$dr), 
                    as.integer(unlist(lapply(objArgs@arguments$Recruitment,function(x) x$modSRactive))),
                    lapply(objArgs@arguments$Recruitment,function(x) as.double(c(rep(x$parAmodSR,length=nT),rep(x$parBmodSR,length=nT),
                            rep(x$parCmodSR,length=nT),rep(x$wnNOISEmodSR,length=nT),rep(x$noiseTypeSR,length=nT)))),#modif MM 27/08/2013 : permet de d�finir un jeu de param�tres SR par ann�e en vectorisant (� la main) chaque composante
                    lapply(objArgs@arguments$Recruitment,function(x) 
                                as.integer(match(x$typeMODsr,c("Mean","Hockey-Stick","Beverton-Holt","Ricker","Shepherd","Quadratic-HS","Smooth-HS")))),
                    as.double(objArgs@arguments$Gestion$mfm),
                    as.double(TACbyF),                    #devra �tre int�gr�e dans objArgs, dimension = nbF+1
                    parBehav,
                    list(active=as.integer(parOptQuot$active),pxQuIni=as.double(parOptQuot$pxQuIni), pxQuMin=as.double(parOptQuot$pxQuMin), 
                          pxQuMax=as.double(parOptQuot$pxQuMax), lambda=as.double(parOptQuot$lambda), ftol=as.double(parOptQuot$ftol)),                           #fonctionne en conjugaison avec TACbyF
                    list(tolVarTACinf=as.double(tacControl$tolVarTACinf),tolVarTACsup=as.double(tacControl$tolVarTACsup),
                          corVarTACval=as.double(tacControl$corVarTACval),corVarTACnby=as.integer(tacControl$corVarTACnby),
                          Blim=as.double(tacControl$Blim),Bmax=as.double(tacControl$Bmax),BlimTrigger=as.integer(tacControl$BlimTrigger)),
                    newStochPrice,       #liste d'�l�ments esp�ce (pas forc�ment toutes pr�sentes, liste vide aussi possible) 
                                         #de format d�crit par la ligne de code de construction de 'newStochPrice'
                    as.character(objArgs@arguments$Replicates$SELECTvar) 
              )
              

if (objArgs@arguments$Replicates$active==1) {     #objet de classe 'iamOutputRep'

  if (is.na(desc)) desc <- "My iamOutputRep object"
  
  return(new("iamOutputRep", desc=desc, arguments=objArgs@arguments, specific=objArgs@specific, 
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
                    GVL_f_m_e = out$GVL_fme),      
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

    return(new("iamOutput", desc=desc, arguments=objArgs@arguments, specific=objArgs@specific,
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
                    Lcm = out$Lcm,             
                    P = out$P,             
                    GVL_f_m_e = out$E$GVL_f_m_e,
                    PQuot = out$PQuot),    
                output = list(
                  typeGest = out$typeGest,
                  nbv_f = out$Eff$nbv_f,              
                  nbds_f = out$Eff$nbds_f,                 
                  nbv_f_m = out$Eff$nbv_f_m,          
                  nbds_f_m = out$Eff$nbds_f_m,               
                  Lbio_f = out$E$Lbio_f,             
                  GVLtot_f_m = out$E$GVLtot_f_m,            
                  GVLav_f_m = out$E$GVLav_f_m,          
                  GVLtot_f = out$E$GVLtot_f,                
                  GVLav_f = out$E$GVLav_f,          
                  GVLoths_f = out$GVLoths_f,         
                  NGVLav_f_m = out$E$NGVLav_f_m,           
                  NGVLav_f = out$E$NGVLav_f,            
                  vcst_f_m = out$E$vcst_f_m,         
                  vcst_f = out$E$vcst_f,             
                  rtbs_f_m = out$E$rtbs_f_m,      
                  rtbs_f = out$E$rtbs_f,             
                  rtbsAct_f = out$E$rtbsAct_f,          
                  cshrT_f_m = out$E$cshrT_f_m,            
                  cshrT_f = out$E$cshrT_f,             
                  sshr_f_m = out$E$sshr_f_m,           
                  sshr_f = out$E$sshr_f,              
                  ncshr_f = out$E$ncshr_f,                
                  ocl_f = out$E$ocl_f,             
                  cs_f = out$E$cs_f,               
                  csAct_f = out$E$csAct_f,          
                  csTot_f = out$E$csTot_f,
                  gva_f = out$E$gva_f,                
                  gvaAct_f = out$E$gvaAct_f,           
                  ccw_f = out$E$ccw_f,        
                  ccwCr_f = out$E$ccwCr_f,            
                  wageg_f = out$E$wageg_f,             
                  wagen_f = out$E$wagen_f,           
                  gcf_f = out$E$gcf_f,             
                  gcfAct_f = out$E$gcfAct_f,           
                  ngcf_f = out$E$ngcf_f,             
                  gp_f = out$E$gp_f,                   
                  ssTot_f = out$E$ssTot_f,             
                  ps_f = out$E$ps_f,              
                  psAct_f = out$E$psAct_f,               
                  sts_f = out$E$sts_f,           
                  stsAct_f = out$E$stsAct_f,              
                  ber_f = out$E$ber_f,        
                  ratio_gva_GVL_f = out$E$ratio_gva_GVL_f,        
                  ratio_gcf_GVL_f = out$E$ratio_gcf_GVL_f,        
                  ratio_fc_GVL_f = out$E$ratio_fc_GVL_f,      
                  ratio_oilc_GVL_f = out$E$ratio_oilc_GVL_f,    
                  ratio_bc_GVL_f = out$E$ratio_bc_GVL_f,      
                  ratio_foc_GVL_f = out$E$ratio_foc_GVL_f,        
                  ratio_icec_GVL_f = out$E$ratio_icec_GVL_f,      
                  ratio_gc_GVL_f = out$E$ratio_gc_GVL_f,        
                  ratio_vc_GVL_f = out$E$ratio_vc_GVL_f,       
                  ratio_rep_GVL_f = out$E$ratio_rep_GVL_f,       
                  ratio_mngc_GVL_f = out$E$ratio_mngc_GVL_f,     
                  ratio_licc_GVL_f = out$E$ratio_licc_GVL_f,      
                  ratio_fvol_GVL_f = out$E$ratio_fvol_GVL_f,     
                  ratio_fvol_Lbio_f = out$E$ratio_fvol_Lbio_f,     
                  ratio_fvol_gva_f = out$E$ratio_fvol_gva_f,      
                  ratio_gcf_gva_f = out$E$ratio_gcf_gva_f,       
                  ratio_K_cnb_f = out$E$ratio_K_cnb_f,         
                  ratio_GVL_K_f = out$E$ratio_GVL_K_f,      
                  ratio_gcf_K_f = out$E$ratio_gcf_K_f,        
                  ratio_ngcf_K_f = out$E$ratio_ngcf_K_f,      
                  ratio_gp_K_f = out$E$ratio_gp_K_f,           
                  ratio_GVL_cnb_ue_f = out$E$ratio_GVL_cnb_ue_f)   
  ))
}              
                                 
})



#---------------
#Examples
#---------------


#out <- IAM.input("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Inputs_SIAD_SEL_2.xls",t_init=2010,nbStep=21)
#
#arg <- IAM.args(out)
#
#mod <- IAM.model(arg,out)
#
