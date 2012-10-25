#WARNING : assumption --> for each fleet f, sum_m vrbl_f_m = vrbl_f





#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
     
#  "iamOutputRep"

#  GLOBAL  #
#----------#

iamToFLR.StCtch.rep <- function(out.rep,inp,type="Stocks") {      #or type="Catches"

  require(FLCore)
  it <- out.rep@arguments$Replicates$nbIter
  specs <- out.rep@specific$Species
  
  flst <- list()
  
  for (spp in specs) {
    #must change age codification
    Age <- suppressWarnings(as.numeric(out.rep@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
    fls <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out.rep@specific$times,iter=1:it)))
  #catch
    if (!is.null(out.rep@outputSp$Ytot[[1]][[spp]])) 
          for (j in 1:it) catch(fls)[,,,,,j] <- apply(out.rep@outputSp$Ytot[[j]][[spp]],2,sum,na.rm=TRUE)
    catch(fls)@units <- "tonnes"
  #catch.n
    if (!is.null(out.rep@outputSp$Ctot[[1]][[spp]])) 
          for (j in 1:it) catch.n(fls)[,,,,,j] <- out.rep@outputSp$Ctot[[j]][[spp]]/1000
    catch.n(fls)@units <- "thousands"
  #catch.wt
    if (!is.null(out.rep@outputSp$Ytot[[1]][[spp]]) & !is.null(out.rep@outputSp$Ctot[[1]][[spp]]))
          for (j in 1:it) catch.wt(fls)[,,,,,j] <- out.rep@outputSp$Ytot[[j]][[spp]]/out.rep@outputSp$Ctot[[j]][[spp]]
    catch.wt(fls)@units <- "kg"
  #discards.n (d_i is needed)
    di <- inp@input[[spp]]$d_i
    if (attributes(di)$DimCst[2]>0) 
      print(paste(spp,"species! 'd_i' parameter is defined per metier! No discards & no landings values in FLStock object!")) 
    if (!is.null(out.rep@outputSp$Ctot[[1]][[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0) 
      for (j in 1:it) discards.n(fls)[,,,,,j] <- catch.n(fls)[,,,,,j] * rep(di,length=length(catch.n(fls)[,,,,,j]))
    discards.n(fls)@units <- "thousands"
  #discards (wD_i should also be used here) 
    wDi <- inp@input[[spp]]$wD_i
    if (!is.null(out.rep@outputSp$Ytot[[1]][[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0 & !all(is.na(wDi))) 
      for (j in 1:it) discards(fls)[,,,,,j] <- apply(discards.n(fls)[,,,,,j] * rep(wDi,length=length(discards.n(fls)[,,,,,j])),2,sum,na.rm=TRUE)  
    discards(fls)@units <- "tonnes"
  #discards.wt
    if (!all(is.na(wDi)))
      for (j in 1:it) discards.wt(fls)[,,,,,j] <- as.vector(wDi)
    discards.wt(fls)@units <- "kg"
  #landings.n
    if (!is.null(out.rep@outputSp$Ctot[[1]][[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0) 
      for (j in 1:it) landings.n(fls)[,,,,,j] <- catch.n(fls)[,,,,,j] - discards.n(fls)[,,,,,j]
    landings.n(fls)@units <- "thousands"
  #landings 
    if (!is.null(out.rep@outputSp$Ytot[[1]][[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0 & !all(is.na(wDi))) 
      for (j in 1:it) landings(fls)[,,,,,j] <- catch(fls)[,,,,,j] - discards(fls)[,,,,,j]  
    landings(fls)@units <- "tonnes"
  #landings.wt
    wLi <- inp@input[[spp]]$wL_i
    if (!all(is.na(wLi)))
      for (j in 1:it) landings.wt(fls)[,,,,,j] <- as.vector(wLi)
    landings.wt(fls)@units <- "kg"
  #stock
    if (!is.null(out.rep@outputSp$B[[1]][[spp]])) 
      for (j in 1:it) stock(fls)[,,,,,j] <- as.vector(out.rep@outputSp$B[[j]][[spp]])
    stock(fls)@units <- "tonnes"        
  #stock.n
    if (!is.null(out.rep@outputSp$N[[1]][[spp]])) 
      for (j in 1:it) stock.n(fls)[,,,,,j] <- as.vector(out.rep@outputSp$N[[j]][[spp]])/1000
    stock(fls)@units <- "thousands"        
  #stock.wt
    wSi <- inp@input[[spp]]$wStock_i
    if (!all(is.na(wSi)))
      for (j in 1:it) stock.wt(fls)[,,,,,j] <- as.vector(wSi)
    stock.wt(fls)@units <- "kg"
  #m
    mnat <- inp@input[[spp]]$M_i
    if (!all(is.na(mnat)))
      for (j in 1:it) m(fls)[,,,,,j] <- as.vector(mnat)
  #mat
    mat <- inp@input[[spp]]$mat_i
    if (!all(is.na(mat)))
      for (j in 1:it) mat(fls)[,,,,,j] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
    if (!is.null(out.rep@outputSp$Z[[1]][[spp]]) & !all(is.na(mnat))) 
      for (j in 1:it) harvest(fls)[,,,,,j] <- as.vector(out.rep@outputSp$Z[[j]][[spp]]) - as.vector(mnat)
    harvest(fls)@units <- "f"        
  #m.spwn
    m.spwn(fls)[] <- 0
  #harvest.spwn
    harvest.spwn(fls)[] <- 0
   
    flst[[spp]] <- fls
  }
  
  if (type=="Catches") {
    warning("No prices in output object! Complete, please!")
    return(FLCatches(lapply(flst,as,'FLCatch')))
  } else {
    return(FLStocks(flst)) 
  }
 
} 
 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#  Per Fleet  #
#-------------#

iamToFLR.Fl.rep <- function(out.rep,inp) {      

  require(FLCore)
  it <- out.rep@arguments$Replicates$nbIter
  specs <- out.rep@specific$Species
  fleet <- out.rep@specific$Fleet
  
  if (is.null(out.rep@output$vcst_f[[1]])) { #ecoDCF
    VCOST <- lapply(1:it,function(x) (out.rep@output$GVLav_f[[x]]*as.vector(1-inp@input$Fleet$lc_f/100) - out.rep@output$rtbs_f[[x]])/out.rep@output$nbds_f[[x]])
  } else {                                   #full Eco
    VCOST <- lapply(1:it,function(x) out.rep@output$vcst_f[[x]]/out.rep@output$nbds_f[[x]])
  }
  CREWSH <- lapply(1:it,function(x) out.rep@output$ccwCr_f[[x]]*as.vector(inp@input$Fleet$cnb_f))
  CPCTY <- out.rep@output$nbv_f
  FCOST <- lapply(1:it,function(x) out.rep@output$rtbs_f[[x]] - out.rep@output$ccwCr_f[[x]]*as.vector(inp@input$Fleet$cnb_f) - out.rep@output$gp_f[[x]])
  
  #Catches per fleet
  flfl <- list()
  
  for (fl in fleet) {
  
    flst.fm <- list()
  
    for (spp in specs) {
  #must change age codification
      Age <- suppressWarnings(as.numeric(out.rep@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
      fls1 <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out.rep@specific$times,iter=1:it)))
  #catch
      test1 <- attributes(out.rep@outputSp$Y[[1]][[spp]])$DimCst[2]==0
      if (!is.null(out.rep@outputSp$Y[[1]][[spp]])) 
        for (j in 1:it) catch(fls1)[,,,,,j] <- apply(out.rep@outputSp$Y[[j]][[spp]][fl,,,],ifelse(test1,2,3),sum,na.rm=TRUE)
      catch(fls1)@units <- "tonnes"
  #catch.n
      if (!is.null(out.rep@outputSp$C[[1]][[spp]])) 
        for (j in 1:it) catch.n(fls1)[,,,,,j] <- apply(out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000,if (test1) 1:2 else 2:3,sum,na.rm=TRUE)
      catch.n(fls1)@units <- "thousands"
  #catch.wt
      if (!is.null(out.rep@outputSp$Y[[1]][[spp]]) & !is.null(out.rep@outputSp$C[[1]][[spp]]))
        for (j in 1:it) catch.wt(fls1)[,,,,,j] <- apply(out.rep@outputSp$Y[[j]][[spp]][fl,,,],if (test1) 1:2 else 2:3,sum,na.rm=TRUE)/
                                                    apply(out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000,if (test1) 1:2 else 2:3,sum,na.rm=TRUE)
      catch.wt(fls1)@units <- "kg"
  #discards.n
      di <- inp@input[[spp]]$d_i
      try(if (attributes(di)$DimCst[2]>0 & attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]==0) print(paste(spp,"species! 'd_i' parameter should not be defined per metier! No discards & no landings values in FLFLeet object!")),silent=TRUE)   #n'arrive pas selon les conditions d'application
      if (!is.null(out.rep@outputSp$C[[1]][[spp]])) 
        if (!all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]==0) & attributes(di)$DimCst[4]==0) {
       
       #2 cases on C
          if (attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]>0) {
       
       #4 cases on di: 
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]==0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000) * as.vector(rep(di[fl,],each=attributes(out.rep@outputSp$C[[j]][[spp]])$DimCst[2])),2:3,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]>0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000) * as.vector(di[fl,,]),2:3,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000) * as.vector(rep(di[],each=attributes(out.rep@outputSp$C[[j]][[spp]])$DimCst[2])),2:3,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]>0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,,,]/1000) * as.vector(di[]),2:3,sum,na.rm=TRUE)
          }
       
          if (attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]==0) { #so, attributes(di)$DimCst[2]==0
       
       #2 cases on di: 
            if (attributes(di)$DimCst[1]>0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- (out.rep@outputSp$C[[j]][[spp]][fl,,]/1000) * as.vector(di[fl,])
            if (attributes(di)$DimCst[1]==0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- (out.rep@outputSp$C[[j]][[spp]][fl,,]/1000) * as.vector(di[])
          }
       
        }
    
      discards.n(fls1)@units <- "thousands"
  
  #discards  
      wDi <- inp@input[[spp]]$wD_i
      if (!is.null(out.rep@outputSp$D[[j]][[spp]]))
        if (!all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out.rep@outputSp$D[[1]][[spp]])$DimCst[2]==0))  
          for (j in 1:it) discards(fls1)[,,,,,j] <- apply(out.rep@outputSp$D[[j]][[spp]],if (test1) c(1,3) else c(1,4),sum,na.rm=TRUE)[fl,] 
      discards(fls1)@units <- "tonnes"
  #discards.wt
      if (!all(is.na(wDi)))
        for (j in 1:it) discards.wt(fls1)[,,,,,j] <- as.vector(wDi)
      discards.wt(fls1)@units <- "kg"
  #landings.n
      if (!is.null(out.rep@outputSp$C[[1]][[spp]])) 
        if (!all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]==0)) 
          for (j in 1:it) landings.n(fls1)[,,,,,j] <- catch.n(fls1)[,,,,,j] - discards.n(fls1)[,,,,,j]
      landings.n(fls1)@units <- "thousands"
  #landings 
      if (!is.null(out.rep@outputSp$C[[1]][[spp]]))
        if (!all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out.rep@outputSp$C[[1]][[spp]])$DimCst[2]==0)) 
          for (j in 1:it) landings(fls1)[,,,,,j] <- catch(fls1)[,,,,,j] - discards(fls1)[,,,,,j]  
      landings(fls1)@units <- "tonnes"
  #landings.wt
      wLi <- inp@input[[spp]]$wL_i
      if (!all(is.na(wLi)))
        for (j in 1:it) landings.wt(fls1)[,,,,,j] <- as.vector(wLi)
      landings.wt(fls1)@units <- "kg"
  #stock
      if (!is.null(out.rep@outputSp$B[[1]][[spp]])) 
        for (j in 1:it) stock(fls1)[,,,,,j] <- as.vector(out.rep@outputSp$B[[j]][[spp]])
      stock(fls1)@units <- "tonnes"        
  #stock.n
      if (!is.null(out.rep@outputSp$N[[1]][[spp]])) 
        for (j in 1:it) stock.n(fls1)[,,,,,j] <- as.vector(out.rep@outputSp$N[[j]][[spp]])/1000
      stock(fls1)@units <- "thousands"        
  #stock.wt
      wSi <- inp@input[[spp]]$wStock_i
      if (!all(is.na(wSi)))
        for (j in 1:it) stock.wt(fls1)[,,,,,j] <- as.vector(wSi)
      stock.wt(fls1)@units <- "kg"
  #m
      mnat <- inp@input[[spp]]$M_i
      if (!all(is.na(mnat)))
        for (j in 1:it) m(fls1)[,,,,,j] <- as.vector(mnat)
  #mat
      mat <- inp@input[[spp]]$mat_i
      if (!all(is.na(mat)))
        for (j in 1:it) mat(fls1)[,,,,,j] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
      if (!is.null(out.rep@outputSp$Fr[[1]][[spp]])) {
        if (attributes(out.rep@outputSp$Fr[[1]][[spp]])$DimCst[2]==0) 
          for (j in 1:it) harvest(fls1)[,,,,,j] <- out.rep@outputSp$Fr[[j]][[spp]][fl,,]
        if (attributes(out.rep@outputSp$Fr[[1]][[spp]])$DimCst[2]>0) 
          for (j in 1:it) harvest(fls1)[,,,,,j] <- apply(out.rep@outputSp$Fr[[j]][[spp]][fl,,,],2:3,sum,na.rm=TRUE)
      }
      harvest(fls1)@units <- "f"        
  #m.spwn
      m.spwn(fls1)[] <- 0
  #harvest.spwn
      harvest.spwn(fls1)[] <- 0
  
  
      flst.fm[[spp]] <- fls1
    }
  
  #catches for the fleet
    FLCS.fm <- FLCatches(lapply(flst.fm,as,'FLCatch'))
  
  #we must consider one 'global' métier level
    met.Glob <- FLMetier(FLCS.fm)
    met.Glob@gear <- "all"
    met.Glob@effshare[] <- 1
    met.Glob@vcost[] <- unlist(lapply(VCOST,function(x) x[fl,]))    
  
    effort.Glob <- FLQuant(unlist(lapply(out.rep@output$nbds_f,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it)) 
    fcost.Glob <- FLQuant(unlist(lapply(FCOST,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it)) 
    capacity.Glob <- FLQuant(unlist(lapply(CPCTY,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it))
    crewshare.Glob <- FLQuant(unlist(lapply(CREWSH,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it))
    
    flfl[[fl]] <- FLFleet(effort = effort.Glob, fcost = fcost.Glob, capacity = capacity.Glob, crewshare = crewshare.Glob, met.Glob)
  }
  
  warning("No prices in output object! Complete, please!")
  return(FLFleets(flfl))

}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#  Per Fleet & Metier (eco) #
#---------------------------#

iamToFLR.FlMet.rep <- function(out.rep,inp) {      

  require(FLCore)
  it <- out.rep@arguments$Replicates$nbIter
  specs <- out.rep@specific$Species
  fleet <- out.rep@specific$Fleet
  metEco <- out.rep@specific$MetierEco
  
  if (!is.null(out.rep@output$vcst_f_m[[1]]))            #not ecoDCF
    VCOSTm <- lapply(1:it,function(x) out.rep@output$vcst_f_m[[x]]/out.rep@output$nbds_f_m[[x]])
  
  #Catches per fleet & per metier
  flfl <- list()
  
  for (fl in fleet) {
  
    flmet <- list()
  
    for (met in metEco) {
  
      flst.fm <- list()
  
      for (spp in specs) {
  
  #matrix conversion 
        MATmm <- inp@input[[spp]]$mm
        nbMAT <- match(met,metEco)
        fls1 <- NULL
  
        if (any(nbMAT%in%MATmm[fl,])) {
  
  #must change age codification
          Age <- suppressWarnings(as.numeric(out.rep@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
          fls1 <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out.rep@specific$times,iter=1:it)))
  #catch
          test1 <- attributes(out.rep@outputSp$Y[[1]][[spp]])$DimCst[2]==0
          if (!is.null(out.rep@outputSp$Y[[1]][[spp]]) & !test1) 
            for (j in 1:it) catch(fls1)[,,,,,j] <- apply(out.rep@outputSp$Y[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],4,sum,na.rm=TRUE)
          catch(fls1)@units <- "tonnes"
  #catch.n
          if (!is.null(out.rep@outputSp$C[[1]][[spp]]) & !test1) 
            for (j in 1:it) catch.n(fls1)[,,,,,j] <- apply(out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000,3:4,sum,na.rm=TRUE)
          catch.n(fls1)@units <- "thousands"
  #catch.wt
          if (!is.null(out.rep@outputSp$Y[[1]][[spp]]) & !is.null(out.rep@outputSp$C[[1]][[spp]]) & !test1)
            for (j in 1:it) catch.wt(fls1)[,,,,,j] <- apply(out.rep@outputSp$Y[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/as.vector(out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,]/1000),3:4,mean,na.rm=TRUE)
          catch.wt(fls1)@units <- "kg"
  #discards.n       
          di <- inp@input[[spp]]$d_i
          if (!is.null(out.rep@outputSp$C[[1]][[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[4]==0 & !test1) {
       
       #4 cases on di: 
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]==0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(rep(di[fl,],each=sum(MATmm[fl,]%in%nbMAT))),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]>0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[fl,MATmm[fl,]%in%nbMAT,]),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[]),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]>0)   
              for (j in 1:it) discards.n(fls1)[,,,,,j] <- apply((out.rep@outputSp$C[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[MATmm[fl,]%in%nbMAT,]),3:4,sum,na.rm=TRUE) 
          }
    
          discards.n(fls1)@units <- "thousands"
  #discards  
          wDi <- inp@input[[spp]]$wD_i
          if (!is.null(out.rep@outputSp$D[[j]][[spp]]) & !all(is.na(di)) & !test1)  
            for (j in 1:it) discards(fls1)[,,,,,j] <- apply(out.rep@outputSp$D[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],4,sum,na.rm=TRUE) 
          discards(fls1)@units <- "tonnes"
  #discards.wt
          if (!all(is.na(wDi)))
            for (j in 1:it) discards.wt(fls1)[,,,,,j] <- as.vector(wDi)
          discards.wt(fls1)@units <- "kg"
  #landings.n
          if (!is.null(out.rep@outputSp$C[[1]][[spp]]) & !all(is.na(di)) & !test1) 
            for (j in 1:it) landings.n(fls1)[,,,,,j] <- catch.n(fls1)[,,,,,j] - discards.n(fls1)[,,,,,j]
          landings.n(fls1)@units <- "thousands"
  #landings 
          if (!is.null(out.rep@outputSp$C[[1]][[spp]]) & !all(is.na(di)) & !test1) 
            for (j in 1:it) landings(fls1)[,,,,,j] <- catch(fls1)[,,,,,j] - discards(fls1)[,,,,,j]  
          landings(fls1)@units <- "tonnes"
  #landings.wt
          wLi <- inp@input[[spp]]$wL_i
          if (!all(is.na(wLi)))
            for (j in 1:it) landings.wt(fls1)[,,,,,j] <- as.vector(wLi)
          landings.wt(fls1)@units <- "kg"
  #stock
          if (!is.null(out.rep@outputSp$B[[1]][[spp]])) 
            for (j in 1:it) stock(fls1)[,,,,,j] <- as.vector(out.rep@outputSp$B[[j]][[spp]])
          stock(fls1)@units <- "tonnes"        
  #stock.n
          if (!is.null(out.rep@outputSp$N[[1]][[spp]])) 
            for (j in 1:it) stock.n(fls1)[,,,,,j] <- as.vector(out.rep@outputSp$N[[j]][[spp]])/1000
          stock(fls1)@units <- "thousands"        
  #stock.wt
          wSi <- inp@input[[spp]]$wStock_i
          if (!all(is.na(wSi)))
            for (j in 1:it) stock.wt(fls1)[,,,,,j] <- as.vector(wSi)
          stock.wt(fls1)@units <- "kg"
  #m
          mnat <- inp@input[[spp]]$M_i
          if (!all(is.na(mnat)))
            for (j in 1:it) m(fls1)[,,,,,j] <- as.vector(mnat)
  #mat
          mat <- inp@input[[spp]]$mat_i
          if (!all(is.na(mat)))
            for (j in 1:it) mat(fls1)[,,,,,j] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
          if (!is.null(out.rep@outputSp$Fr[[1]][[spp]]) & !test1) 
            for (j in 1:it) harvest(fls1)[,,,,,j] <- apply(out.rep@outputSp$Fr[[j]][[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],3:4,sum,na.rm=TRUE)
    
          harvest(fls1)@units <- "f"        
  #m.spwn
          m.spwn(fls1)[] <- 0
  #harvest.spwn
          harvest.spwn(fls1)[] <- 0
        }
  
        if (!is.null(fls1)) flst.fm[[spp]] <- fls1
  
      }
  
      if (length(flst.fm)>0) {
    #catches for the fleet-metier
        FLCS.fm <- FLCatches(lapply(flst.fm,as,'FLCatch'))
    #'métier' level
        met.comp <- FLMetier(FLCS.fm)
      } else {
        met.comp <- FLMetier(FLQuant(NA,dimnames=list(year=out.rep@specific$times,iter=1:it)))
      }                                                                                              
  
      met.comp@gear <- met
      for (j in 1:it) met.comp@effshare[,,,,,j] <- out.rep@output$nbds_f_m[[j]][fl,met,]/as.vector(out.rep@output$nbds_f[[j]][fl,])
      if (!is.null(out.rep@output$vcst_f_m[[1]])) {
        for (j in 1:it) met.comp@vcost[,,,,,j] <- unlist(lapply(VCOSTm,function(x) x[fl,met,]))    
      } else {
        for (j in 1:it) met.comp@vcost[,,,,,j] <- NA
      }
      flmet[[met]] <- met.comp
    }
  
    CREWSH <- lapply(1:it,function(x) out.rep@output$ccwCr_f[[x]]*as.vector(inp@input$Fleet$cnb_f))
    CPCTY <- out.rep@output$nbv_f
    FCOST <- lapply(1:it,function(x) out.rep@output$rtbs_f[[x]] - out.rep@output$ccwCr_f[[x]]*as.vector(inp@input$Fleet$cnb_f) - out.rep@output$gp_f[[x]])
  
    effort.comp <- FLQuant(unlist(lapply(out.rep@output$nbds_f,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it)) 
    fcost.comp <- FLQuant(unlist(lapply(FCOST,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it)) 
    capacity.comp <- FLQuant(unlist(lapply(CPCTY,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it))
    crewshare.comp <- FLQuant(unlist(lapply(CREWSH,function(x) x[fl,])), 
                              dimnames=list(year=out.rep@specific$times,iter=1:it))
    
    flfl[[fl]] <- FLFleet(effort = effort.comp, fcost = fcost.comp, capacity = capacity.comp, crewshare = crewshare.comp, FLMetiers(flmet))
  }
  
  warning("No prices in output object! Complete, please!")
  return(FLFleets(flfl))

}

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
     
#  "iamOutput"

#  GLOBAL  #
#----------#

#  FLStock(s)

iamToFLR.StCtch <- function(out,inp,type="Stocks") {      #or type="Catches"

  require(FLCore)
  specs <- out@specific$Species
  
  flst <- list()
  
  for (spp in specs) {
    #must change age codification
    Age <- suppressWarnings(as.numeric(out@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
    fls <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out@specific$times,iter=1)))
  #catch
    if (!is.null(out@outputSp$Ytot[[spp]])) 
      catch(fls)[,,,,,] <- apply(out@outputSp$Ytot[[spp]],2,sum,na.rm=TRUE)
    catch(fls)@units <- "tonnes"
  #catch.n
    if (!is.null(out@outputSp$Ctot[[spp]])) 
      catch.n(fls)[,,,,,] <- out@outputSp$Ctot[[spp]]/1000
    catch.n(fls)@units <- "thousands"
  #catch.wt
    if (!is.null(out@outputSp$Ytot[[spp]]) & !is.null(out@outputSp$Ctot[[spp]]))
      catch.wt(fls)[,,,,,] <- out@outputSp$Ytot[[spp]]/out@outputSp$Ctot[[spp]]
    catch.wt(fls)@units <- "kg"
  #discards.n (d_i is needed)
    di <- inp@input[[spp]]$d_i
    if (attributes(di)$DimCst[2]>0) print(paste(spp,"species! 'd_i' parameter is defined per metier! No discards & no landings values in FLStock object!")) 
    if (!is.null(out@outputSp$Ctot[[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0) 
      discards.n(fls)[,,,,,] <- catch.n(fls)[,,,,,] * rep(di,length=length(catch.n(fls)[,,,,,]))
    discards.n(fls)@units <- "thousands"
  #discards (wD_i should also be used here) 
    wDi <- inp@input[[spp]]$wD_i
    if (!is.null(out@outputSp$Ytot[[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0 & !all(is.na(wDi))) 
      discards(fls)[,,,,,] <- apply(discards.n(fls)[,,,,,] * rep(wDi,length=length(discards.n(fls)[,,,,,])),2,sum,na.rm=TRUE)  
    discards(fls)@units <- "tonnes"
  #discards.wt
    if (!all(is.na(wDi)))
      discards.wt(fls)[,,,,,] <- as.vector(wDi)
    discards.wt(fls)@units <- "kg"
  #landings.n
    if (!is.null(out@outputSp$Ctot[[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0) 
      landings.n(fls)[,,,,,] <- catch.n(fls)[,,,,,] - discards.n(fls)[,,,,,]
    landings.n(fls)@units <- "thousands"
  #landings 
    if (!is.null(out@outputSp$Ytot[[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0 & !all(is.na(wDi))) 
      landings(fls)[,,,,,] <- catch(fls)[,,,,,] - discards(fls)[,,,,,]  
    landings(fls)@units <- "tonnes"
  #landings.wt
    wLi <- inp@input[[spp]]$wL_i
    if (!all(is.na(wLi)))
      landings.wt(fls)[,,,,,] <- as.vector(wLi)
    landings.wt(fls)@units <- "kg"
  #stock
    if (!is.null(out@outputSp$B[[spp]])) 
      stock(fls)[,,,,,] <- as.vector(out@outputSp$B[[spp]])
    stock(fls)@units <- "tonnes"        
  #stock.n
    if (!is.null(out@outputSp$N[[spp]])) 
      stock.n(fls)[,,,,,] <- as.vector(out@outputSp$N[[spp]])/1000
    stock(fls)@units <- "thousands"        
  #stock.wt
    wSi <- inp@input[[spp]]$wStock_i
    if (!all(is.na(wSi)))
      stock.wt(fls)[,,,,,] <- as.vector(wSi)
    stock.wt(fls)@units <- "kg"
  #m
    mnat <- inp@input[[spp]]$M_i
    if (!all(is.na(mnat)))
      m(fls)[,,,,,] <- as.vector(mnat)
  #mat
    mat <- inp@input[[spp]]$mat_i
    if (!all(is.na(mat)))
      mat(fls)[,,,,,] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
    if (!is.null(out@outputSp$Z[[spp]]) & !all(is.na(mnat))) 
      harvest(fls)[,,,,,] <- as.vector(out@outputSp$Z[[spp]]) - as.vector(mnat)
    harvest(fls)@units <- "f"        
  #m.spwn
    m.spwn(fls)[] <- 0
  #harvest.spwn
    harvest.spwn(fls)[] <- 0
  
    flst[[spp]] <- fls
  }
   
  if (type=="Catches") {
    warning("No prices in output object! Complete, please!")
    return(FLCatches(lapply(flst,as,'FLCatch')))
  } else {
    return(FLStocks(flst))
  }

} 
 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#  Per Fleet  #
#-------------#

iamToFLR.Fl <- function(out,inp) {      

  require(FLCore)
  it <- out@arguments$Replicates$nbIter
  specs <- out@specific$Species
  fleet <- out@specific$Fleet

  if (is.null(out@output$vcst_f)) { #ecoDCF
    VCOST <- (out@output$GVLav_f*as.vector(1-inp@input$Fleet$lc_f/100) - out@output$rtbs_f)/out@output$nbds_f
  } else {                                   #full Eco
    VCOST <- out@output$vcst_f/out@output$nbds_f
  }
  
  CREWSH <- out@output$ccwCr_f*as.vector(inp@input$Fleet$cnb_f)
  CPCTY <- out@output$nbv_f
  FCOST <- out@output$rtbs_f - out@output$ccwCr_f*as.vector(inp@input$Fleet$cnb_f) - out@output$gp_f
  
  #Catches per fleet
  flfl <- list()
  
  for (fl in fleet) {
  
    flst.fm <- list()
    
    for (spp in specs) {
  #must change age codification
      Age <- suppressWarnings(as.numeric(out@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
      fls1 <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out@specific$times,iter=1)))
  #catch
      test1 <- attributes(out@outputSp$Y[[spp]])$DimCst[2]==0
      if (!is.null(out@outputSp$Y[[spp]])) 
        catch(fls1)[,,,,,] <- apply(out@outputSp$Y[[spp]][fl,,,],ifelse(test1,2,3),sum,na.rm=TRUE)
      catch(fls1)@units <- "tonnes"
  #catch.n
      if (!is.null(out@outputSp$C[[spp]])) 
        catch.n(fls1)[,,,,,] <- apply(out@outputSp$C[[spp]][fl,,,]/1000,if (test1) 1:2 else 2:3,sum,na.rm=TRUE)
      catch.n(fls1)@units <- "thousands"
  #catch.wt
      if (!is.null(out@outputSp$Y[[spp]]) & !is.null(out@outputSp$C[[spp]]))
        catch.wt(fls1)[,,,,,] <- apply(out@outputSp$Y[[spp]][fl,,,],if (test1) 1:2 else 2:3,sum,na.rm=TRUE)/
                                   apply(out@outputSp$C[[spp]][fl,,,]/1000,if (test1) 1:2 else 2:3,sum,na.rm=TRUE)
      catch.wt(fls1)@units <- "kg"
  #discards.n
      di <- inp@input[[spp]]$d_i
      if (attributes(di)$DimCst[2]>0 & attributes(out@outputSp$C[[spp]])$DimCst[2]==0) print(paste(spp,"species! 'd_i' parameter should not be defined per metier! No discards & no landings values in FLFLeet object!"))   #n'arrive pas selon les conditions d'application
      if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out@outputSp$C[[spp]])$DimCst[2]==0) & attributes(di)$DimCst[4]==0) {
       
       #2 cases on C
        if (attributes(out@outputSp$C[[spp]])$DimCst[2]>0) {
       
       #4 cases on di: 
          if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]==0)   
            discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,,,]/1000) * as.vector(rep(di[fl,],each=attributes(out@outputSp$C[[spp]])$DimCst[2])),2:3,sum,na.rm=TRUE)
          if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]>0)   
            discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,,,]/1000) * as.vector(di[fl,,]),2:3,sum,na.rm=TRUE)
          if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0)   
            discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,,,]/1000) * as.vector(rep(di[],each=attributes(out@outputSp$C[[spp]])$DimCst[2])),2:3,sum,na.rm=TRUE)
          if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]>0)   
            discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,,,]/1000) * as.vector(di[]),2:3,sum,na.rm=TRUE)
       
        }
       
        if (attributes(out@outputSp$C[[spp]])$DimCst[2]==0) { #so, attributes(di)$DimCst[2]==0
       
       #2 cases on di: 
          if (attributes(di)$DimCst[1]>0)   
            discards.n(fls1)[,,,,,] <- (out@outputSp$C[[spp]][fl,,]/1000) * as.vector(di[fl,])
          if (attributes(di)$DimCst[1]==0)   
            discards.n(fls1)[,,,,,] <- (out@outputSp$C[[spp]][fl,,]/1000) * as.vector(di[])
        }
       
      }
    
      discards.n(fls1)@units <- "thousands"
  #discards  
      wDi <- inp@input[[spp]]$wD_i
      if (!is.null(out@outputSp$D[[spp]]) & !all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out@outputSp$C[[spp]])$DimCst[2]==0))  
        discards(fls1)[,,,,,] <- apply(out@outputSp$D[[spp]],if (test1) c(1,3) else c(1,4),sum,na.rm=TRUE)[fl,] 
      discards(fls1)@units <- "tonnes"
  #discards.wt
      if (!all(is.na(wDi)))
        discards.wt(fls1)[,,,,,] <- as.vector(wDi)
      discards.wt(fls1)@units <- "kg"
  #landings.n
      if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out@outputSp$C[[spp]])$DimCst[2]==0)) 
        landings.n(fls1)[,,,,,] <- catch.n(fls1)[,,,,,] - discards.n(fls1)[,,,,,]
      landings.n(fls1)@units <- "thousands"
  #landings 
      if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & !(attributes(di)$DimCst[2]>0 & attributes(out@outputSp$C[[spp]])$DimCst[2]==0)) 
        landings(fls1)[,,,,,] <- catch(fls1)[,,,,,] - discards(fls1)[,,,,,]  
      landings(fls1)@units <- "tonnes"
  #landings.wt
      wLi <- inp@input[[spp]]$wL_i
      if (!all(is.na(wLi)))
        landings.wt(fls1)[,,,,,] <- as.vector(wLi)
      landings.wt(fls1)@units <- "kg"
  #stock
      if (!is.null(out@outputSp$B[[spp]])) 
        stock(fls1)[,,,,,] <- as.vector(out@outputSp$B[[spp]])
      stock(fls1)@units <- "tonnes"        
  #stock.n
      if (!is.null(out@outputSp$N[[spp]])) 
        stock.n(fls1)[,,,,,] <- as.vector(out@outputSp$N[[spp]])/1000
      stock(fls1)@units <- "thousands"        
  #stock.wt
      wSi <- inp@input[[spp]]$wStock_i
      if (!all(is.na(wSi)))
        stock.wt(fls1)[,,,,,] <- as.vector(wSi)
      stock.wt(fls1)@units <- "kg"
  #m
      mnat <- inp@input[[spp]]$M_i
      if (!all(is.na(mnat)))
        m(fls1)[,,,,,] <- as.vector(mnat)
  #mat
      mat <- inp@input[[spp]]$mat_i
      if (!all(is.na(mat)))
        mat(fls1)[,,,,,] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
      if (!is.null(out@outputSp$Fr[[spp]])) {
        if (attributes(out@outputSp$Fr[[spp]])$DimCst[2]==0) 
          harvest(fls1)[,,,,,] <- out@outputSp$Fr[[spp]][fl,,]
        if (attributes(out@outputSp$Fr[[spp]])$DimCst[2]>0) 
          harvest(fls1)[,,,,,] <- apply(out@outputSp$Fr[[spp]][fl,,,],2:3,sum,na.rm=TRUE)
      }
      harvest(fls1)@units <- "f"        
  #m.spwn
      m.spwn(fls1)[] <- 0
  #harvest.spwn
      harvest.spwn(fls1)[] <- 0
  
      flst.fm[[spp]] <- fls1
    }
  
  #catches for the fleet
    FLCS.fm <- FLCatches(lapply(flst.fm,as,'FLCatch'))
  
  #we must consider one 'global' métier level
    met.Glob <- FLMetier(FLCS.fm)
    met.Glob@gear <- "all"
    met.Glob@effshare[] <- 1
    met.Glob@vcost[] <- VCOST[fl,]    
  
    effort.Glob <- FLQuant(out@output$nbds_f[fl,], 
                              dimnames=list(year=out@specific$times,iter=1)) 
    fcost.Glob <- FLQuant(FCOST[fl,], 
                              dimnames=list(year=out@specific$times,iter=1)) 
    capacity.Glob <- FLQuant(CPCTY[fl,], 
                              dimnames=list(year=out@specific$times,iter=1))
    crewshare.Glob <- FLQuant(CREWSH[fl,], 
                              dimnames=list(year=out@specific$times,iter=1))
    
    flfl[[fl]] <- FLFleet(effort = effort.Glob, fcost = fcost.Glob, capacity = capacity.Glob, crewshare = crewshare.Glob, met.Glob)
  }
  
  warning("No prices in output object! Complete, please!")
  return(FLFleets(flfl))
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#  Per Fleet & Metier (eco) #
#---------------------------#

iamToFLR.FlMet <- function(out,inp) {      

  require(FLCore)
  it <- out@arguments$Replicates$nbIter
  specs <- out@specific$Species
  fleet <- out@specific$Fleet
  metEco <- out@specific$MetierEco

  if (!is.null(out@output$vcst_f_m))            #not ecoDCF
    VCOSTm <- out@output$vcst_f_m/out@output$nbds_f_m

  CREWSH <- out@output$ccwCr_f*as.vector(inp@input$Fleet$cnb_f)
  CPCTY <- out@output$nbv_f
  FCOST <- out@output$rtbs_f - out@output$ccwCr_f*as.vector(inp@input$Fleet$cnb_f) - out@output$gp_f

  #Catches per fleet & per metier
  flfl <- list()
  
  for (fl in fleet) {
  
    flmet <- list()
  
    for (met in metEco) {
  
      flst.fm <- list()
  
      for (spp in specs) {
  
  #matrix conversion 
        MATmm <- inp@input[[spp]]$mm
        nbMAT <- match(met,metEco)
        fls1 <- NULL
  
        if (any(nbMAT%in%MATmm[fl,])) {
  
  #must change age codification
          Age <- suppressWarnings(as.numeric(out@specific$Ages[[spp]])) ; nas <- sum(is.na(Age)) ; if (nas>0) Age[is.na(Age)] <- max(Age,na.rm=TRUE)+(1:nas)
          fls1 <- FLStock(FLQuant(NA,dimnames=list(age=Age,year=out@specific$times,iter=1)))
  #catch
          test1 <- attributes(out@outputSp$Y[[spp]])$DimCst[2]==0
          if (!is.null(out@outputSp$Y[[spp]]) & !test1) 
            catch(fls1)[,,,,,] <- apply(out@outputSp$Y[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],4,sum,na.rm=TRUE)
          catch(fls1)@units <- "tonnes"
  #catch.n
          if (!is.null(out@outputSp$C[[spp]]) & !test1) 
            catch.n(fls1)[,,,,,] <- apply(out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000,3:4,sum,na.rm=TRUE)
          catch.n(fls1)@units <- "thousands"
  #catch.wt
          if (!is.null(out@outputSp$Y[[spp]]) & !is.null(out@outputSp$C[[spp]]) & !test1)
            catch.wt(fls1)[,,,,,] <- apply(out@outputSp$Y[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/as.vector(out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,]/1000),3:4,mean,na.rm=TRUE)
          catch.wt(fls1)@units <- "kg"
  #discards.n       
          di <- inp@input[[spp]]$d_i
          if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & attributes(di)$DimCst[4]==0 & !test1) {
       
       #4 cases on di: 
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]==0)   
              discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(rep(di[fl,],each=sum(MATmm[fl,]%in%nbMAT))),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]>0 & attributes(di)$DimCst[2]>0)   
              discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[fl,MATmm[fl,]%in%nbMAT,]),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]==0)   
              discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[]),3:4,sum,na.rm=TRUE)
            if (attributes(di)$DimCst[1]==0 & attributes(di)$DimCst[2]>0)   
              discards.n(fls1)[,,,,,] <- apply((out@outputSp$C[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE]/1000) * as.vector(di[MATmm[fl,]%in%nbMAT,]),3:4,sum,na.rm=TRUE) 
          }
    
          discards.n(fls1)@units <- "thousands"
  #discards  
          wDi <- inp@input[[spp]]$wD_i
          if (!is.null(out@outputSp$D[[spp]]) & !all(is.na(di)) & !test1)  
            discards(fls1)[,,,,,] <- apply(out@outputSp$D[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],4,sum,na.rm=TRUE) 
          discards(fls1)@units <- "tonnes"
  #discards.wt
          if (!all(is.na(wDi)))
            discards.wt(fls1)[,,,,,] <- as.vector(wDi)
          discards.wt(fls1)@units <- "kg"
  #landings.n
          if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & !test1) 
            landings.n(fls1)[,,,,,] <- catch.n(fls1)[,,,,,] - discards.n(fls1)[,,,,,]
          landings.n(fls1)@units <- "thousands"
  #landings 
          if (!is.null(out@outputSp$C[[spp]]) & !all(is.na(di)) & !test1) 
            landings(fls1)[,,,,,] <- catch(fls1)[,,,,,] - discards(fls1)[,,,,,]  
          landings(fls1)@units <- "tonnes"
  #landings.wt
          wLi <- inp@input[[spp]]$wL_i
          if (!all(is.na(wLi)))
            landings.wt(fls1)[,,,,,] <- as.vector(wLi)
          landings.wt(fls1)@units <- "kg"
  #stock
          if (!is.null(out@outputSp$B[[spp]])) 
            stock(fls1)[,,,,,] <- as.vector(out@outputSp$B[[spp]])
          stock(fls1)@units <- "tonnes"        
  #stock.n
          if (!is.null(out@outputSp$N[[spp]])) 
            stock.n(fls1)[,,,,,] <- as.vector(out@outputSp$N[[spp]])/1000
          stock(fls1)@units <- "thousands"        
  #stock.wt
          wSi <- inp@input[[spp]]$wStock_i
          if (!all(is.na(wSi)))
            stock.wt(fls1)[,,,,,] <- as.vector(wSi)
          stock.wt(fls1)@units <- "kg"
  #m
          mnat <- inp@input[[spp]]$M_i
          if (!all(is.na(mnat)))
            m(fls1)[,,,,,] <- as.vector(mnat)
  #mat
          mat <- inp@input[[spp]]$mat_i
          if (!all(is.na(mat)))
            mat(fls1)[,,,,,] <- as.vector(mat)
  #harvest (Fr, not F ???) --> considering survivals 
          if (!is.null(out@outputSp$Fr[[spp]]) & !test1) 
            harvest(fls1)[,,,,,] <- apply(out@outputSp$Fr[[spp]][fl,MATmm[fl,]%in%nbMAT,,,drop=FALSE],3:4,sum,na.rm=TRUE)
    
          harvest(fls1)@units <- "f"        
  #m.spwn
          m.spwn(fls1)[] <- 0
  #harvest.spwn
          harvest.spwn(fls1)[] <- 0
        }
  
        if (!is.null(fls1)) flst.fm[[spp]] <- fls1
  
      }
  
      if (length(flst.fm)>0) {
    #catches for the fleet-metier
        FLCS.fm <- FLCatches(lapply(flst.fm,as,'FLCatch'))
    #'métier' level
        met.comp <- FLMetier(FLCS.fm)
      } else {
        met.comp <- FLMetier(FLQuant(NA,dimnames=list(year=out@specific$times,iter=1)))
      }                                                                                              
  
      met.comp@gear <- met
      met.comp@effshare[,,,,,] <- out@output$nbds_f_m[fl,met,]/as.vector(out@output$nbds_f[fl,])
      if (!is.null(out@output$vcst_f_m)) {
        met.comp@vcost[,,,,,] <- VCOSTm[fl,met,]    
      } else {
        met.comp@vcost[,,,,,] <- NA
      }
      flmet[[met]] <- met.comp
    }
  
    effort.comp <- FLQuant(out@output$nbds_f[fl,], 
                              dimnames=list(year=out@specific$times,iter=1)) 
    fcost.comp <- FLQuant(FCOST[fl,], 
                              dimnames=list(year=out@specific$times,iter=1)) 
    capacity.comp <- FLQuant(CPCTY[fl,], 
                              dimnames=list(year=out@specific$times,iter=1))
    crewshare.comp <- FLQuant(CREWSH[fl,], 
                              dimnames=list(year=out@specific$times,iter=1))
    
    flfl[[fl]] <- FLFleet(effort = effort.comp, fcost = fcost.comp, capacity = capacity.comp, crewshare = crewshare.comp, FLMetiers(flmet))
  }
  
  warning("No prices in output object! Complete, please!")
  return(FLFleets(flfl))
  
}

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------


#date()
#
#aa <- iamToFLR.StCtch.rep(outRep,inp,"Stocks")
#bb <- iamToFLR.StCtch.rep(outRep,inp,"Catches")    
#cc <- iamToFLR.Fl.rep(outRep,inp)    
#dd <- iamToFLR.FlMet.rep(outRep,inp)     
#
#date()
#     
#aa2 <- iamToFLR.StCtch(out,inp,"Stocks")
#bb2 <- iamToFLR.StCtch(out,inp,"Catches")    
#cc2 <- iamToFLR.Fl(out,inp)    
#dd2 <- iamToFLR.FlMet(out,inp)     
#
#date()



setGeneric("iamToFLR", function(object,input,type="Stock",...){      #type="Stock","Catch","Fleet","FleetMetier"

  standardGeneric("iamToFLR")

})



setMethod("iamToFLR", signature("iamOutputRep","iamInput"), function(object,input,type="Stock",...){

if (type=="Stock") return(iamToFLR.StCtch.rep(object,input,"Stocks"))
if (type=="Catch") return(iamToFLR.StCtch.rep(object,input,"Catches"))
if (type=="Fleet") return(iamToFLR.Fl.rep(object,input))
if (type=="FleetMetier") return(iamToFLR.FlMet.rep(object,input))
if (!type%in%c("Stock","Catch","Fleet","FleetMetier")) return(NULL)

})         
  
  
  
                                                                     
setMethod("iamToFLR", signature("iamOutput","iamInput"), function(object,input,type="Stock",...){

if (type=="Stock") return(iamToFLR.StCtch(object,input,"Stocks"))
if (type=="Catch") return(iamToFLR.StCtch(object,input,"Catches"))
if (type=="Fleet") return(iamToFLR.Fl(object,input))
if (type=="FleetMetier") return(iamToFLR.FlMet(object,input))
if (!type%in%c("Stock","Catch","Fleet","FleetMetier")) return(NULL)

})         



