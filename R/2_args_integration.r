

GUIiam <- function(INPUT,classInp="args",desc=as.character(NA)) {     #or classInp="param"


require(tcltk)
require(tcltk2)
e1 <- new.env()

ALLVarRep <- c("B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth","mu_nbds","mu_nbv","N","Eff",
                 "GVL_fme","GVLtot_fm","GVLav_f","vcst_fm","vcst_f","rtbs_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
                 "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f","GVLtot_f","wagen_f","L_efmit","D_efmit",
                 "Fr_fmi","C_efmit")

if (classInp=="args") {

  input <- INPUT@arguments
  
  spp <- names(input$Rec)
  years <- as.numeric(colnames(input$Gest$tac_fbar))  
  fleets <- colnames(input$Gest$mf)
  metiers <- colnames(input$Gest$mfm)
  ALLscenario <- input$Scen$ALLscenario
  SELECTscen <- input$Scen$SELECTscen
  SELECTvar <- input$Rep$SELECTvar         #character avec intitulés des variables
  tacfbar <- rbind(TAC=input$Gest$tac,Fbar=input$Gest$fbar) 
  level <- input$Gest$level
  matF <- matrix(input$Gest$mf,nrow=1) ; colnames(matF) <- names(input$Gest$mf)
  matFM <- input$Gest$mfm

} else {

  input <- INPUT

  spp <- input@specific$Species 
  RecIni <- unlist(lapply(input@input,function(x) as.numeric(x$N_it0[1]))[spp])
  RecIniQ <- unlist(lapply(input@input,function(x) if (length(x$N_it0)>=4) as.numeric(sum(x$N_it0[1:4])) else NA)[spp])
  RecIni[input@specific$trim%in%1] <- RecIniQ[input@specific$trim%in%1] 
  years <- input@specific$times 
  fleets <- input@specific$Fleet
  metier <- input@specific$MetierEco
  ALLscenario <- names(input@scenario) 
  SELECTscen <- 1
  SELECTvar <- c("B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth","mu_nbds","mu_nbv","N","Eff",
                 "GVL_fme","GVLtot_fm","GVLav_f","rtbs_f","vcst_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
                 "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f","GVLtot_f","wagen_f")
  tacfbar <- matrix(0,nrow=2,ncol=length(years),dimnames=list(c("TAC","Fbar"),years))
  level <- 0
  matF <- matrix(1,nrow=1,ncol=length(fleets),dimnames=list(NULL,fleets))
  matFM <- input@input$Fleet$nbv_f_m ; matFM[!is.na(matFM)] <- 1

}


BASE <- tktoplevel(height=500,width=1000)     #interface entière
tkwm.title(BASE,"IAM")

base1 <- tkframe(BASE)   
base2 <- tkframe(BASE)
base2_1 <- tkframe(base2)
#frame Economic  ---------------------------------------------------------------

  #variables

if (classInp=="args") {  
  
  EcoDisable <- tclVar(as.character(input$Eco$active))
  Type <- tclVar(as.character(input$Eco$type))
  Adj <- tclVar(as.character(input$Eco$adj))
  Lev <- tclVar(as.character(input$Eco$lev))
  Ue_choice <- tclVar(as.character(input$Eco$ue_choice))
  Oths <- tclVar(as.character(input$Eco$oths))
  OthsFM <- tclVar(as.character(input$Eco$othsFM))
  PerscCalc <- tclVar(as.character(input$Eco$perscCalc))
  Report <- tclVar(as.character(input$Eco$report))
  Dr <- tclVar(as.character(input$Eco$dr))

} else {

  EcoDisable <- tclVar("0")
  Type <- tclVar("1")
  Adj <- tclVar("1")
  Lev <- tclVar("1")
  Ue_choice <- tclVar("1")
  Oths <- tclVar("0")
  OthsFM <- tclVar("0")
  PerscCalc <- tclVar("0")
  Report <- tclVar("0")
  Dr <- tclVar("0")

} 

  #frame Eco

set.eco.state<-function(){
    on.off<-tclvalue(EcoDisable)
    if(on.off=="1") {tkconfigure(butType1, state="normal") ;
                     tkconfigure(butType2, state="normal");
                     tkconfigure(butAdj1, state="normal");
                     tkconfigure(butAdj2, state="normal");
                     tkconfigure(butLev1, state="normal");
                     tkconfigure(butLev2, state="normal");
                     tkconfigure(butUe1, state="normal");
                     tkconfigure(butUe2, state="normal");
                     tkconfigure(butOths1, state="normal");
                     tkconfigure(butOths2, state="normal");
                     tkconfigure(butOthsFM1, state="normal");
                     tkconfigure(butOthsFM2, state="normal");
                     tkconfigure(butPerscCalc1, state="normal");
                     tkconfigure(butPerscCalc2, state="normal");
                     tkconfigure(butPerscCalc3, state="normal");
                     tkconfigure(butReport1, state="normal");
                     tkconfigure(butReport2, state="normal");
                     tkconfigure(spinboxDR, state="normal")}

    if(on.off=="0") {tkconfigure(butType1, state="disabled") ;
                     tkconfigure(butType2, state="disabled");
                     tkconfigure(butAdj1, state="disabled");
                     tkconfigure(butAdj2, state="disabled");
                     tkconfigure(butLev1, state="disabled");
                     tkconfigure(butLev2, state="disabled");
                     tkconfigure(butUe1, state="disabled");
                     tkconfigure(butUe2, state="disabled");
                     tkconfigure(butOths1, state="disabled");
                     tkconfigure(butOths2, state="disabled");
                     tkconfigure(butOthsFM1, state="disabled");
                     tkconfigure(butOthsFM2, state="disabled");
                     tkconfigure(butPerscCalc1, state="disabled");
                     tkconfigure(butPerscCalc2, state="disabled");
                     tkconfigure(butPerscCalc3, state="disabled");
                     tkconfigure(butReport1, state="disabled");
                     tkconfigure(butReport2, state="disabled");
                     tkconfigure(spinboxDR, state="disabled")}
}
butEco <- ttkcheckbutton(base2_1,text="Economic",variable=EcoDisable,command=set.eco.state)
frmEco <- ttklabelframe(base2_1,labelwidget=butEco,borderwidth=2,height=500,width=300)

  #Type
frmType <- tkframe(frmEco)
butType1 <- ttkradiobutton(frmType,text="Complete",value="1",variable=Type,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butType2 <- ttkradiobutton(frmType,text="DCF",value="2",variable=Type,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butType1,ttklabel(frmType,text=""),butType2,side="left")

  #adj
frmAdj <- tkframe(frmEco)
butAdj1 <- ttkradiobutton(frmAdj,text="1",value="1",variable=Adj,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butAdj2 <- ttkradiobutton(frmAdj,text="2",value="2",variable=Adj,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butAdj1,ttklabel(frmAdj,text="      "),butAdj2,side="left")

  #lev
frmLev <- tkframe(frmEco)
butLev1 <- ttkradiobutton(frmLev,text="1",value="1",variable=Lev,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butLev2 <- ttkradiobutton(frmLev,text="2",value="2",variable=Lev,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butLev1,ttklabel(frmLev,text="      "),butLev2,side="left")

  #ue_choice
frmUe <- tkframe(frmEco)
butUe1 <- ttkradiobutton(frmUe,text="1",value="1",variable=Ue_choice,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butUe2 <- ttkradiobutton(frmUe,text="2",value="2",variable=Ue_choice,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butUe1,ttklabel(frmUe,text="      "),butUe2,side="left")

  #oths
frmOths <- tkframe(frmEco)
butOths1 <- ttkradiobutton(frmOths,text="0",value="0",variable=Oths,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butOths2 <- ttkradiobutton(frmOths,text="1",value="1",variable=Oths,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butOths1,ttklabel(frmOths,text="      "),butOths2,side="left")

  #othsFM
frmOthsFM <- tkframe(frmEco)
butOthsFM1 <- ttkradiobutton(frmOthsFM,text="0",value="0",variable=OthsFM,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butOthsFM2 <- ttkradiobutton(frmOthsFM,text="1",value="1",variable=OthsFM,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butOthsFM1,ttklabel(frmOthsFM,text="      "),butOthsFM2,side="left")

  #perscCalc
frmPerscCalc <- tkframe(frmEco)
butPerscCalc1 <- ttkradiobutton(frmPerscCalc,text="0",value="0",variable=PerscCalc,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butPerscCalc2 <- ttkradiobutton(frmPerscCalc,text="1",value="1",variable=PerscCalc,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butPerscCalc3 <- ttkradiobutton(frmPerscCalc,text="2",value="2",variable=PerscCalc,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))

tkpack(butPerscCalc1,ttklabel(frmPerscCalc,text="  "),butPerscCalc2,ttklabel(frmPerscCalc,text="  "),butPerscCalc3,side="left")

  #report
frmReport <- tkframe(frmEco)
butReport1 <- ttkradiobutton(frmReport,text="0",value="0",variable=Report,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
butReport2 <- ttkradiobutton(frmReport,text="1",value="1",variable=Report,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(butReport1,ttklabel(frmReport,text="      "),butReport2,side="left")

  #taux d'actualisation
frmDR <- tkframe(frmEco)
spinboxDR <- tk2spinbox(frmDR,from=-2,to=2,increment=0.001,width=9,
              state=ifelse(tclvalue(EcoDisable)=="0","disabled","normal"))
tkpack(spinboxDR,side="left")
tkconfigure(spinboxDR, textvariable = Dr)

tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="Type"),frmType)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="adj"),frmAdj)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="lev"),frmLev)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="ue_choice"),frmUe)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="oths"),frmOths)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="othsFM"),frmOthsFM)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="perscCalc"),frmPerscCalc)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="report"),frmReport)
tkgrid(tk2label(frmEco,text="    "))
tkgrid(tk2label(frmEco,text="Taux d'actualisation"),frmDR)
tkgrid(tk2label(frmEco,text="    "))

#tkpack(frmEco)

#-------------------------------------------------------------------------------





#frame Gestion -----------------------------------------------------------------

if (classInp=="args") {  
 
  #variables
  GestDisable <- tclVar(as.character(input$Gest$active))
  Controle <- tclVar(as.character(input$Gest$control))
  Target <- tclVar(as.character(input$Gest$target))
  Espece <- tclVar(as.character(input$Gest$espece))
  Level <- tclVar(as.character(input$Gest$level))
  Delay <- tclVar(as.character(input$Gest$delay))
  Update <- tclVar(as.character(input$Gest$upd))
  Super <- tclVar(as.character(input$Gest$sup))
  Infer <- tclVar(as.character(input$Gest$inf))
  assign("TACFbar",tacfbar,envir=e1) 
  assign("mF",matF,envir=e1)
  assign("mFM",matFM,envir=e1) 

} else {

    #variables
  GestDisable <- tclVar("0")
  Controle <- tclVar("Nb navires")
  Target <- tclVar("TAC")
  Espece <- tclVar(spp[1])
  Level <- tclVar("0")
  Delay <- tclVar("2")
  Update <- tclVar("1")
  Super <- tclVar("0")
  Infer <- tclVar("0")
  assign("TACFbar",tacfbar,envir=e1) 
  assign("mF",matF,envir=e1)
  assign("mFM",matFM,envir=e1) 

}

  #frame Gestion

set.ges.state<-function(){
    on.off<-tclvalue(GestDisable)
    if(on.off=="1") {tkconfigure(comboControl, state="normal") ;
                     tkconfigure(comboTarget, state="normal");
                     tkconfigure(comboSpecies, state="normal");
                     tkconfigure(butLvl0, state="normal");
                     tkconfigure(butLvl1, state="normal");
                     tkconfigure(spinboxDelay, state="normal");
                     tkconfigure(butUpd1, state="normal");
                     tkconfigure(butUpd2, state="normal");
                     tkconfigure(spinboxSup, state="normal");
                     tkconfigure(spinboxInf, state="normal");
                     tkconfigure(TFbut, state="normal");
                     tkconfigure(MFbut, state="normal")}

    if(on.off=="0") {tkconfigure(comboControl, state="disabled") ;
                     tkconfigure(comboTarget, state="disabled");
                     tkconfigure(comboSpecies, state="disabled");
                     tkconfigure(butLvl0, state="disabled");
                     tkconfigure(butLvl1, state="disabled");
                     tkconfigure(spinboxDelay, state="disabled");
                     tkconfigure(butUpd1, state="disabled");
                     tkconfigure(butUpd2, state="disabled");
                     tkconfigure(spinboxSup, state="disabled");
                     tkconfigure(spinboxInf, state="disabled");
                     tkconfigure(TFbut, state="disabled");
                     tkconfigure(MFbut, state="disabled")}
}

butGes <- ttkcheckbutton(base2_1,text="Gestion",variable=GestDisable,command=set.ges.state)
frmGes <- ttklabelframe(base2_1,labelwidget=butGes,borderwidth=2,height=500,width=300)

  #Contrôle
comboControl <- tk2combobox(frmGes,state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkconfigure(comboControl, textvariable = Controle, values=c("Nb navires","Nb jdm"))


  #Cible
comboTarget <- tk2combobox(frmGes,state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkconfigure(comboTarget, textvariable = Target, values=c("TAC","Fbar","TAC->Fbar"))

  #espèce
comboSpecies <- tk2combobox(frmGes,state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkconfigure(comboSpecies, textvariable = Espece, values=spp)

  #level
frmLvl <- tkframe(frmGes)
butLvl0 <- ttkradiobutton(frmLvl,text="f",value="0",variable=Level,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
butLvl1 <- ttkradiobutton(frmLvl,text="fm",value="1",variable=Level,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkpack(butLvl0,ttklabel(frmLvl,text="      "),butLvl1,side="left")

  #delay
spinboxDelay <- tk2spinbox(frmGes,from=1,to=15,increment=1,width=9,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkconfigure(spinboxDelay, textvariable = Delay)

  #update
frmUpd <- tkframe(frmGes)
butUpd1 <- ttkradiobutton(frmUpd,text="Oui",value="1",variable=Update,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
butUpd2 <- ttkradiobutton(frmUpd,text="Non",value="2",variable=Update,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))
tkpack(butUpd1,ttklabel(frmUpd,text="      "),butUpd2,side="left")

  #borne sup
spinboxSup <- tk2spinbox(frmGes,from=-10,to=10,increment=0.1,width=9,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"),
              command=function() {if (as.numeric(tclvalue(Super))<as.numeric(tclvalue(Infer)))
                                      tclvalue(Infer) <- tclvalue(Super)})
tkconfigure(spinboxSup, textvariable = Super)

  #borne inf
spinboxInf <- tk2spinbox(frmGes,from=-10,to=10,increment=0.1,width=9,
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"),
              command=function() {if (as.numeric(tclvalue(Infer))>as.numeric(tclvalue(Super)))
                                      tclvalue(Super) <- tclvalue(Infer)})
tkconfigure(spinboxInf, textvariable = Infer)

  #TAC/Fbar
TFbut <- tk2button(frmGes, text="TAC/Fbar", width=10, command=function() {temp <- get("TACFbar",envir=e1)
                                                                          res <- fix(temp)
                                                                          assign("TACFbar",res,envir=e1)
                                                                          tkfocus(BASE)},
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))

  #mF
MFbut <- tk2button(frmGes, text="mF/mFM", width=10, command=function() {
                                                                    if (tclvalue(Level)=="0") {
                                                                      temp <- get("mF",envir=e1)
                                                                      res <- fix(temp)
                                                                      assign("mF",res,envir=e1)
                                                                      tkfocus(BASE)
                                                                    } else {
                                                                      temp <- get("mFM",envir=e1)
                                                                      res <- fix(temp)
                                                                      assign("mFM",res,envir=e1)
                                                                      tkfocus(BASE)                                                      
                                                                    }},
              state=ifelse(tclvalue(GestDisable)=="0","disabled","normal"))


tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Contrôle"),comboControl)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Cible"),comboTarget)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Espèce"),comboSpecies)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Level"),frmLvl)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Delay"),spinboxDelay)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Update"),frmUpd)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Borne supérieure"),spinboxSup)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="Borne inférieure"),spinboxInf)
tkgrid(tk2label(frmGes,text="  "))
tkgrid(tk2label(frmGes,text="  "))
tkgrid(TFbut,MFbut)
tkgrid(tk2label(frmGes,text="  "))

#tkpack(frmEco,frmGes,side="left",fill="both")

#-------------------------------------------------------------------------------




#frame Recrutement -------------------------------------------------------------



  #variables
  

RecDisable <- tclVar("0")

    #par espèce
ll <- lapply(spp,function(x) NA)
ModSRDisable <- SimuStochDisable <- TypeModSR <- ParAModSR <- ParBModSR <- ParCModSR <- WNoiseModSR <- ParNoiseSR <- TypeSimuStoch <- ll
names(ModSRDisable) <- names(SimuStochDisable) <- names(TypeModSR) <- names(ParAModSR) <- names(ParBModSR) <-
  names(ParCModSR) <- names(WNoiseModSR) <- names(ParNoiseSR) <- names(TypeSimuStoch) <- spp

if (classInp=="args") {

  for (i in spp) {
    ModSRDisable[[i]] <- tclVar(as.character(input$Rec[[i]]$modSRactive))
    SimuStochDisable[[i]] <- tclVar(as.character(input$Rec[[i]]$simuSTOCHactive))
    TypeModSR[[i]] <- tclVar(as.character(input$Rec[[i]]$typeMODsr))
    ParAModSR[[i]] <- tclVar(as.character(input$Rec[[i]]$parAmodSR))
    ParBModSR[[i]] <- tclVar(as.character(input$Rec[[i]]$parBmodSR))
    ParCModSR[[i]] <- tclVar(as.character(input$Rec[[i]]$parCmodSR))
    WNoiseModSR[[i]] <- tclVar(as.character(input$Rec[[i]]$wnNOISEmodSR))
    ParNoiseSR[[i]] <- tclVar(as.character(input$Rec[[i]]$noiseTypeSR))
    TypeSimuStoch[[i]] <- tclVar(as.character(input$Rec[[i]]$typeSIMUstoch))
  }


} else {

  for (i in spp) {
    ModSRDisable[[i]] <- tclVar("1")
    SimuStochDisable[[i]] <- tclVar("0")
    TypeModSR[[i]] <- tclVar("Mean")
    ParAModSR[[i]] <- tclVar(RecIni[i])
    ParBModSR[[i]] <- tclVar("0")
    ParCModSR[[i]] <- tclVar("0")
    WNoiseModSR[[i]] <- tclVar("0")
    ParNoiseSR[[i]] <- tclVar("1")
    TypeSimuStoch[[i]] <- tclVar("1")
  }

}

  #frame Recrutement

#set.rec.state<-function(){
#    on.off<-tclvalue(RecDisable)
#    if(on.off=="1") {}
#
#    if(on.off=="0") {}
#}

set.modSR.state<-function(index){
    on.off<-tclvalue(ModSRDisable[[index]])
    if(on.off=="1") {eval(parse('',text=paste("tkconfigure(comboTypeRec",index,", state=\"normal\");",     
                     "tkconfigure(spinboxParA",index,", state=\"normal\");",
                     "tkconfigure(spinboxParB",index,", state=\"normal\");",
                     "tkconfigure(spinboxParC",index,", state=\"normal\");",
                     "tkconfigure(spinboxWN",index,", state=\"normal\");",
                     "tkconfigure(butWN",index,", state=\"normal\");",
                     "tkconfigure(butWLN",index,", state=\"normal\");",
                     "tkconfigure(butType1_",index,", state=\"disabled\");",
                     "tkconfigure(butType2_",index,", state=\"disabled\");",
                     "tkconfigure(butType3_",index,", state=\"disabled\");",
                     "tclvalue(SimuStochDisable[[",index,"]]) <- \"0\"",sep="")))
                     }

    if(on.off=="0") {eval(parse('',text=paste("tkconfigure(comboTypeRec",index,", state=\"disabled\");",     
                     "tkconfigure(spinboxParA",index,", state=\"disabled\");",
                     "tkconfigure(spinboxParB",index,", state=\"disabled\");",
                     "tkconfigure(spinboxParC",index,", state=\"disabled\");",
                     "tkconfigure(spinboxWN",index,", state=\"disabled\");",
                     "tkconfigure(butWN",index,", state=\"disabled\");",
                     "tkconfigure(butWLN",index,", state=\"disabled\");",
                     "tkconfigure(butType1_",index,", state=\"normal\");",
                     "tkconfigure(butType2_",index,", state=\"normal\");",
                     "tkconfigure(butType3_",index,", state=\"normal\");",
                     "tclvalue(SimuStochDisable[[",index,"]]) <- \"1\"",sep="")))
                     }
}

set.simu.state<-function(index){
    on.off<-tclvalue(SimuStochDisable[[index]])
    if(on.off=="1") {eval(parse('',text=paste("tkconfigure(comboTypeRec",index,", state=\"disabled\");",     
                     "tkconfigure(spinboxParA",index,", state=\"disabled\");",
                     "tkconfigure(spinboxParB",index,", state=\"disabled\");",
                     "tkconfigure(spinboxParC",index,", state=\"disabled\");",
                     "tkconfigure(spinboxWN",index,", state=\"disabled\");",
                     "tkconfigure(butWN",index,", state=\"disabled\");",
                     "tkconfigure(butWLN",index,", state=\"disabled\");",
                     "tkconfigure(butType1_",index,", state=\"normal\");",
                     "tkconfigure(butType2_",index,", state=\"normal\");",
                     "tkconfigure(butType3_",index,", state=\"normal\");",
                     "tclvalue(ModSRDisable[[",index,"]]) <- \"0\"",sep="")))
                     }

    if(on.off=="0") {eval(parse('',text=paste("tkconfigure(comboTypeRec",index,", state=\"normal\");",     
                     "tkconfigure(spinboxParA",index,", state=\"normal\");",
                     "tkconfigure(spinboxParB",index,", state=\"normal\");",
                     "tkconfigure(spinboxParC",index,", state=\"normal\");",
                     "tkconfigure(spinboxWN",index,", state=\"normal\");",
                     "tkconfigure(butWN",index,", state=\"normal\");",
                     "tkconfigure(butWLN",index,", state=\"normal\");",
                     "tkconfigure(butType1_",index,", state=\"disabled\");",
                     "tkconfigure(butType2_",index,", state=\"disabled\");",
                     "tkconfigure(butType3_",index,", state=\"disabled\");",
                     "tclvalue(ModSRDisable[[",index,"]]) <- \"1\"",sep="")))
                     }
}


titRec <- tk2label(base1,text="Recrutement")
frmRec <- ttklabelframe(base1,labelwidget=titRec,borderwidth=2)


nb <- tk2notebook(frmRec, tabs = spp)     #à modifier en frmRec
tkpack(tk2label(frmRec,text="    "),nb, tk2label(frmRec,text="    "),fill = "both", expand = 1,side="top")

for (i in 1:length(spp)){

  eval(parse('',text=paste("tbook",i," <- tk2notetab(nb, spp[",i,"])",sep="")))
  eval(parse('',text=paste("butModSR",i," <- ttkcheckbutton(tbook",i,
            ",text=\"Model SR\",variable=ModSRDisable[[",i,"]],command=function() set.modSR.state(",i,"))",sep="")))
  eval(parse('',text=paste("frmModSR",i," <- ttklabelframe(tbook",i,",labelwidget=butModSR",i,",borderwidth=2)",sep="")))

    #Type
  eval(parse('',text=paste("comboTypeRec",i," <- tk2combobox(tbook",i,
            ",width=12,state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("tkconfigure(comboTypeRec",i,", textvariable = TypeModSR[[",i,
            "]], values=c(\"Mean\",\"Hockey-Stick\",\"Beverton-Holt\",\"Ricker\",\"Shepherd\",\"Quadratic-HS\"))",sep="")))

    #Paramètre a
  eval(parse('',text=paste("spinboxParA",i," <- tk2spinbox(tbook",i,",from=-1000000,to=1000000,increment=0.0001,width=12,",
                    "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
                #tk2entry(tbook1,width=12,state=ifelse(tclvalue(ModSRDisable[[1]])=="0","disabled","normal"))
  eval(parse('',text=paste("tkconfigure(spinboxParA",i,", textvariable = ParAModSR[[",i,"]])",sep="")))

    #Paramètre b
  eval(parse('',text=paste("spinboxParB",i," <- tk2spinbox(tbook",i,",from=-1000000,to=1000000,increment=0.0001,width=12,",
                    "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("tkconfigure(spinboxParB",i,", textvariable = ParBModSR[[",i,"]])",sep="")))

    #Paramètre c
  eval(parse('',text=paste("spinboxParC",i," <- tk2spinbox(tbook",i,",from=-1000000,to=1000000,increment=0.0001,width=12,",
                    "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("tkconfigure(spinboxParC",i,", textvariable = ParCModSR[[",i,"]])",sep="")))

    #Paramètre Wnoise
  eval(parse('',text=paste("spinboxWN",i," <- tk2spinbox(tbook",i,",from=-1000000,to=1000000,increment=0.0001,width=12,",
                    "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("tkconfigure(spinboxWN",i,", textvariable = WNoiseModSR[[",i,"]])",sep="")))


  eval(parse('',text=paste("butWN",i," <- ttkradiobutton(tbook",i,",text=\"Normal\",value=\"1\",variable=ParNoiseSR[[",i,"]],",
              "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("butWLN",i," <- ttkradiobutton(tbook",i,",text=\"LogNormal\",value=\"2\",variable=ParNoiseSR[[",i,"]],",
              "state=ifelse(tclvalue(ModSRDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))

  eval(parse('',text=paste("butSimuSt",i," <- ttkcheckbutton(tbook",i,",text=\"Simu Stoch\",variable=SimuStochDisable[[",i,
                    "]],command=function() set.simu.state(",i,"))",sep="")))
  eval(parse('',text=paste("frmSimuSt",i," <- ttklabelframe(tbook",i,",labelwidget=butSimuSt",i,",borderwidth=2)",sep="")))
    #Type
  eval(parse('',text=paste("butType1_",i," <- ttkradiobutton(frmSimuSt",i,",text=\"1\",value=\"1\",variable=TypeSimuStoch[[",i,
                    "]],state=ifelse(tclvalue(SimuStochDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("butType2_",i," <- ttkradiobutton(frmSimuSt",i,",text=\"2\",value=\"2\",variable=TypeSimuStoch[[",i,
                    "]],state=ifelse(tclvalue(SimuStochDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))
  eval(parse('',text=paste("butType3_",i," <- ttkradiobutton(frmSimuSt",i,",text=\"3\",value=\"3\",variable=TypeSimuStoch[[",i,
                    "]],state=ifelse(tclvalue(SimuStochDisable[[",i,"]])==\"0\",\"disabled\",\"normal\"))",sep="")))

  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\"Type\"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(butType1_",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(butType2_",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(butType3_",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmSimuSt",i,",text=\" \"))",sep="")))

  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Type\"),comboTypeRec",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Paramètre a\"),spinboxParA",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Paramètre b\"),spinboxParB",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Paramètre c\"),spinboxParC",i,")",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Bruit :  \"),butWN",i,",tk2label(frmModSR",i,",text=\"  \"),butWLN",i,")",sep="")))  
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\" \"))",sep="")))
  eval(parse('',text=paste("tkgrid(tk2label(frmModSR",i,",text=\"Ecart-type\"),spinboxWN",i,")",sep="")))
  
  eval(parse('',text=paste("tkgrid(frmModSR",i,",frmSimuSt",i,")",sep="")))

}
  #tkpack(frmRec)
  
#-------------------------------------------------------------------------------



#frame Bootstrap ---------------------------------------------------------------


  #variables

if (classInp=="args") {

  BootDisable <- tclVar(as.character(input$Rep$active))
  NbIt <- tclVar(as.character(input$Rep$nbIter))

} else {

  BootDisable <- tclVar("0")
  NbIt <- tclVar("500")

} 
  #frame Boot

set.boot.state<-function(){
    on.off<-tclvalue(BootDisable)
    if(on.off=="1") {tkconfigure(spinboxBoot, state="normal") ; tkconfigure(menuVar,state="normal")}

    if(on.off=="0") {tkconfigure(spinboxBoot, state="disabled") ; tkconfigure(menuVar,state="disabled")}
}

butBoot <- ttkcheckbutton(base1,text="Bootstrap",variable=BootDisable,command=set.boot.state)
frmBoot <- ttklabelframe(base1,labelwidget=butBoot,borderwidth=2)

  #NbIterations & liste variables
  
frmNbit <- tkframe(frmBoot)
spinboxBoot <- tk2spinbox(frmNbit,from=1,to=5000,increment=1,width=5,
              state=ifelse(tclvalue(BootDisable)=="0","disabled","normal"))
tkpack(spinboxBoot,side="left")
tkconfigure(spinboxBoot, textvariable = NbIt)

scr2 <- tkscrollbar(frmBoot,repeatinterval=5,command=function(...)tkyview(menuVar,...))  
menuVar <- tk2listbox(frmBoot,height=6,width=35,selectmode="multiple",background="white",yscrollcommand=function(...)tkset(scr2,...))

for (j in (1:length(ALLVarRep))) tkinsert(menuVar,"end",ALLVarRep[j])
indicVar <- match(SELECTvar,ALLVarRep)
if (any(is.na(indicVar))) stop("wrong variable name for replicates!!")
for (k in (1:length(indicVar))) tkselection.set(menuVar,indicVar[k]-1)

if (tclvalue(BootDisable)=="0") {tkconfigure(menuVar,state="disabled")} else {tkconfigure(menuVar,state="normal")}

tkgrid(tk2label(frmBoot,text="    "))
tkgrid(tk2label(frmBoot,text="       Nombre d'itérations                    "),frmNbit)
tkgrid(tk2label(frmBoot,text="    "))
tkgrid(tk2label(frmBoot,text="Variables de sortie              "),menuVar,scr2,tk2label(frmBoot,text="  "))
tkgrid(tk2label(frmBoot,text="    "))



#-------------------------------------------------------------------------------



#frame Scenario ---------------------------------------------------------------


  #variables
  
if (classInp=="args") {
 
  ScenDisable <- tclVar(as.character(input$Scen$active))

} else {

  ScenDisable <- tclVar("0")

}
  #frame Scenar

set.scen.state<-function(){
    on.off <- tclvalue(ScenDisable)
    if(on.off=="1") tkconfigure(menuScen, state="normal") 
    if(on.off=="0") tkconfigure(menuScen, state="disabled")
}

butScen <- ttkcheckbutton(base1,text="Scenario",variable=ScenDisable,command=set.scen.state)
frmScen <- ttklabelframe(base1,labelwidget=butScen,borderwidth=2)

  #Scenarii
scr <- tkscrollbar(frmScen,repeatinterval=5,command=function(...)tkyview(menuScen,...))  
menuScen <- tk2listbox(frmScen,height=6,width=35,selectmode="single",background="white",yscrollcommand=function(...)tkset(scr,...))


for (j in (1:length(ALLscenario))) tkinsert(menuScen,"end",ALLscenario[j])
for (k in (1:length(SELECTscen))) tkselection.set(menuScen,SELECTscen[k]-1)

if (tclvalue(ScenDisable)=="0") {tkconfigure(menuScen,state="disabled")} else {tkconfigure(menuScen,state="normal")}

tkgrid(tk2label(frmScen,text="    "))
tkgrid(tk2label(frmScen,text="  "),menuScen,scr,tk2label(frmScen,text="  "))
tkgrid(tk2label(frmScen,text="    "))

#tkpack(frmScen)

#accès aux sélections : tkcurselection(menuScen) 
# --> SELECTscen <- as.numeric(strsplit(tclvalue(tkcurselection(menuScen))," ")[[1]])+1

#-------------------------------------------------------------------------------



#1ère colonne (frame base1) : Recrutement/Bootstrap/Scénario


tkpack(frmRec,frmBoot,frmScen,tk2label(base1,text="    "),side="top",fill="both")     #base1
tkpack(frmGes,frmEco,side="left",fill="both")             #base2_1




#frame Buttons   ---------------------------------------------------------------

RETURN <- function() {

listEco <- list(active = as.integer(tclvalue(EcoDisable)),
               type = as.integer(tclvalue(Type)),
               adj = as.integer(tclvalue(Adj)),
               lev = as.integer(tclvalue(Lev)),
               ue_choice = as.integer(tclvalue(Ue_choice)),
               oths = as.integer(tclvalue(Oths)),
               othsFM = as.integer(tclvalue(OthsFM)),
               perscCalc = as.integer(tclvalue(PerscCalc)),
               report = as.integer(tclvalue(Report)), 
               dr =  as.double(tclvalue(Dr)))

listGestion <- list(active = as.integer(tclvalue(GestDisable)),
                   control = tclvalue(Controle),
                   target = tclvalue(Target),
                   espece = tclvalue(Espece),
                   delay = as.integer(tclvalue(Delay)),
                   level = as.integer(tclvalue(Level)),
                   upd = as.integer(tclvalue(Update)),
                   sup = as.double(tclvalue(Super)),
                   inf = as.double(tclvalue(Infer)),
                   tac = get("TACFbar",envir=e1)["TAC",],
                   fbar = get("TACFbar",envir=e1)["Fbar",],
                   mf = get("mF",envir=e1)[1,],
                   mfm = get("mFM",envir=e1))
                   
attributes(listGestion$tac)$DimCst <- as.integer(c(0,0,0,length(listGestion$tac)))
attributes(listGestion$fbar)$DimCst <- as.integer(c(0,0,0,length(listGestion$fbar)))
attributes(listGestion$mf)$DimCst <- as.integer(c(length(listGestion$mf),0,0,0))
attributes(listGestion$mfm)$DimCst <- as.integer(c(dim(listGestion$mfm),0,0))
                   
listRec <- lapply(spp,function(x) NA) ; names(listRec) <- spp

for (i in 1:length(spp)) {
    
    listRec[[i]] <- list(modSRactive = as.integer(tclvalue(ModSRDisable[[i]])),
                         typeMODsr = tclvalue(TypeModSR[[i]]),
                         parAmodSR = as.double(tclvalue(ParAModSR[[i]])),
                         parBmodSR = as.double(tclvalue(ParBModSR[[i]])),
                         parCmodSR = as.double(tclvalue(ParCModSR[[i]])),
                         wnNOISEmodSR = as.double(tclvalue(WNoiseModSR[[i]])),
                         noiseTypeSR = as.double(tclvalue(ParNoiseSR[[i]])),
                         simuSTOCHactive = as.integer(tclvalue(SimuStochDisable[[i]])),
                         typeSIMUstoch = as.integer(tclvalue(TypeSimuStoch[[i]])))
}


listBoot <- list(active = as.integer(tclvalue(BootDisable)),
                 nbIter = as.integer(tclvalue(NbIt)),
                 SELECTvar = ALLVarRep[as.numeric(strsplit(tclvalue(tkcurselection(menuVar))," ")[[1]])+1])


listScen <- list(active = as.integer(tclvalue(ScenDisable)),
                 ALLscenario = ALLscenario,
                 SELECTscen = as.numeric(strsplit(tclvalue(tkcurselection(menuScen))," ")[[1]])+1)


tkdestroy(BASE)

param <- list(Recruitment=listRec,Replicates=listBoot,Scenario=listScen,Gestion=listGestion,Eco=listEco)
 
if (classInp=="args") {
  if (!is.na(desc)) {newDesc <- desc} else {newDesc <- INPUT@desc}
  OUT <- new("iamArgs", desc=newDesc, arguments=param, specific=INPUT@specific)
} else {
  OUT <- new("iamArgs", desc=desc, arguments=param, specific=INPUT@specific)
}


assign("LL",OUT,envir=e1)
}




But <- tkframe(base2)

OK.but <- tkbutton(But,text="            OK            ",command=RETURN)
Cancel.but <- tkbutton(But,text="         Cancel         ",command=function() {tkdestroy(BASE) ; assign("LL",INPUT,envir=e1)})

tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text="         "))
tkgrid(tk2label(But,text=paste(rep(" ",40),collapse="")),OK.but,
    tk2label(But,text=paste(rep(" ",10),collapse="")),Cancel.but,
    tk2label(But,text=paste(rep(" ",10),collapse="")))
tkgrid(tk2label(But,text="         "))

tkpack(base2_1,But,side="top",fill="both")     #base2
tkpack(base1,base2,side="left",fill="both")     #BASE


#-------------------------------------------------------------------------------

tkwait.window(BASE)
invisible(get("LL",envir=e1))      
}



                                        

#modelIAM <- function(outGUI) {
#
#out <- outGUI
#
#varFoth <- as.double(rep(0,length(out$input@specific$Species)))                                        # à intégrer dans l'interface
#varFoth[match(out$param$Gest$espece,out$input@specific$Species)] <- as.double(1)                             #
#
#
#system.time(outPGS <- .Call("essai", out$input@input, out$input@specific, out$input@stochastic, 
#              eval(parse('',text=paste(c("inp@scenario",out$param$Scen$ALLscenario[out$param$Scen$SELECTscen]),collapse="$"))), #attention, ne marche qu'avec un seul choix de scenario
#              RecType1=unlist(lapply(out$param$Rec,function(x) as.integer(x$modSRactive==0 & x$typeSIMUstoch==1))), #le recrutement stochastique sera de type 1 O/N(tirage aléatoire dans l'historique indépendant par espèce
#              RecType2=unlist(lapply(out$param$Rec,function(x) as.integer(x$modSRactive==0 & x$typeSIMUstoch==2))), #le recrutement stochastique sera de type 2 O/N
#              RecType3=unlist(lapply(out$param$Rec,function(x) as.integer(x$modSRactive==0 & x$typeSIMUstoch==3))), #le recrutement stochastique sera de type 1 O/N
#              as.integer(out$param$Scen$active), #pas de scénario caractérisé dans l'onglet 'Scénario' du fichier de paramètres
#              as.integer(out$param$Boot$active), #on enclenche la procédure itérative
#              as.integer(out$param$Boot$nbIter),#nombre d'itérations
#              as.integer(out$param$Gest$active), #on enclenche le module de Gestion
#              as.double(out$param$Gest$mf), #le taux de variation appliqué à l'effort de chaque flottille pour atteindre le TAC reste constant
#              varFoth, #attention !!!! -> à modifier #on ne fait varier conjointement aux efforts par flottille que la mortalité par pêche de la sole générée par les autres flottilles
#              as.double(c(out$param$Gest$inf,out$param$Gest$sup)), #bornes/contraintes du multiplicateur
#              as.double(out$param$Gest$tac_fbar["TAC",]), #vecteur des valeurs de TAC à atteindre pour chaque année simulée (la variable "délai" déterminera la première année d'application)
#              as.double(out$param$Gest$tac_fbar["Fbar",]), #vecteur des valeurs de Fbar à atteindre pour chaque année simulée (la variable "délai" déterminera la première année d'application)
#              as.integer(c(                        #vecteur des paramètres de gestion
#                           eTemp = match(out$param$Gest$espece,out$input@specific$Species)-1,                          #eTemp = espèce concernée par la mesure de gestion
#                           var = match(out$param$Gest$control,c("Nb jdm","Nb navires")),      #var = variable de controle (1=nbds, 2=nbv)
#                           trgt = match(out$param$Gest$target,c("TAC","Fbar","TAC->Fbar")),   #trgt = variable cible (1=TAC, 2=Fbar)
#                           delay = out$param$Gest$delay,                                      #delay = délai de première application de la mesure de gestion
#                           upd = out$param$Gest$upd)),                  #à corriger                       #upd = expression des multiplicateurs de sorties (appliqués à une valeur initiale=1, ou à la valeur antérieure=2)                                                                              
#              as.integer(out$param$Eco$type-1),  #type de modèle éco appliqué (DCF=1, complet=0)
#              as.integer(c(adj = out$param$Eco$adj, lev = out$param$Eco$lev, ue_choice = out$param$Eco$ue_choice, oths = out$param$Eco$oths, 
#                            othsFM = out$param$Eco$othsFM, perscCalc = out$param$Eco$perscCalc, report = out$param$Eco$report)), #paramètres du module éco
#              as.double(out$param$Eco$dr), #taux d'actualisation à appliquer aux indicateurs éco
#              as.integer(unlist(lapply(out$param$Rec,function(x) as.integer(x$modSRactive==1)))), #emploi d'un modèle de recrutement (incompatible avec recrutement stochastique)
#              lapply(out$param$Rec,function(x) as.double(x[3:6])), #paramètres SR par espèce (par milliers)
#              lapply(out$param$Rec,function(x) as.integer(match(x$typeMODsr,c("Mean","Hockey-Stick","Beverton-Holt","Ricker","Shepherd"))))))  #type de relation SR par espèce
#
#
#invisible(outPGS)
#
#}
#



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Méthodes de constitution des objets 'iamArgs' par interface graphique

#Arguments...

setGeneric("IAM.args", function(object, specific, ...){
	standardGeneric("IAM.args")
	}
)

  # étape d'initialisation

setMethod("IAM.args", signature("iamInput","missing"),function(object, desc=as.character(NA), ...){
	
  return(GUIiam(object,classInp="param",desc=desc))

})

  # ou...
  
setMethod("IAM.args", signature("character","character"),function(object, specific, desc=as.character(NA), ...){
	
  if (substring(object,nchar(object)-3,nchar(object))!=".txt") stop("'object' must be a .txt file!!")
  if (substring(specific,nchar(specific)-3,nchar(specific))!=".txt") stop("'specific' must be a .txt file!!")
  
  specif <- .Call("Fun",normalizePath(specific),NULL)
  argum <- .Call("Fun",normalizePath(object),specif)
  
  return(new("iamArgs", desc=desc, arguments=argum, specific=specif))

})
  


  # étape de modification 

setMethod("IAM.args", signature("iamArgs","missing"),function(object, desc=as.character(NA), ...){
	
  return(GUIiam(object,classInp="args",desc=desc))

})



#---------------
#Examples
#---------------

##importation de l'input 
#out <- IAM.input("Z:/Projet/Projet SIAD/Param bio_eco/Modele/Inputs_SIAD_SEL_2.xls",t_init=2010,nbStep=21)
#
##initialisation des arguments
#argOut <- IAM.args(out)
#
##modification des arguments
#updArg <- IAM.args(argOut)
#
##importation des arguments
#impArg <- IAM.args("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expArgs/args.txt",
#                   "C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expArgs/specific.txt",desc="My args")
#
