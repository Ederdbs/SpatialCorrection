################################################################################
#### Simulation genotype, systematic checks, and analysis function 
#### Eder Silva and Alencar Xavier - September 2020
################################################################################
x <-  c('ggplot2','sp','plyr','reshape','INLA','dplyr','doParallel','GGally',
        'lubridate','parallel','verification','AlphaSimR','bWGR','Matrix',
        'foreach')
xl <- lapply(x, require, character.only = TRUE,quietly=TRUE)
### Covariance  Matern function from Elias  Krainski 2019-----------------------
book.rMatern <- function(n, coords, sigma, range, kappa = sqrt(8*nu)/range,
                         variance = sigma^2, nu=1){
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m) - lgamma(nu))* besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),matrix(rnorm(nrow(coords)*n), ncol=n))))
}
### Vanraden matrix genomic relationship ---------------------------------------
GmatrixM <- function(M){
  Z <- apply(as.matrix(M),2,function(x){x-mean(x)})
  G <- tcrossprod(Z)
  diag(G) <- diag(G)*1.01
  G <- G/mean(diag(G))
  return(G)
}
### simulate diagonal ----------------------------------------------------------
repla_diag <- function(r = 50,p = 50,f = 0.1){
  data <- expand.grid(y=1:r,x=1:p)
  data$type <- 0
  delta <- round(1/f,0)
  starP <- rep(c(seq(1,delta,by=2),seq(2,delta,by=2)),times= round((p/delta)+1))
  for(i in 1:p){
    data[data$x == i & data$y %in% seq(starP[i],r,by=delta),'type'] <- 1
  }
  return(data)
}
### simulate genotype ----------------------------------------------------------
simulate_gen <- function(r=NULL,
                         p=NULL,
                         f=NULL,
                         Nchecks=NULL,
                         sigma2e =NULL,
                         sigma2u = NULL,
                         range =NULL,
                         nu = NULL,
                         mgen=NULL,
                         sigma2ge=NULL,
                         unitD=F,
                         Sp=NULL,
                         Sr=NULL,
                         NF1=NULL,
                         NF2=NULL,
                         Nsnp=NULL,
                         Nrep=NULL){
  data <- repla_diag(r=r,p=p,f=f)
  data$type <- ifelse(data$type == 1,'DC','G')
  data$gen <- NULL
  data$gen[data$type=='DC'] <- paste('Check',sample(1:Nchecks,sum(data$type=='DC'),replace = TRUE),sep='')
  data$gen[data$type=='G'] <- paste('Gen',1:(round(sum(data$type=='G')/Nrep,0)),sep='')
  data$rep <- NA
  for(r in unique(data$gen)){
    data$rep[data$gen == r] <- 1:sum(data$gen == r)
  }
  ### Genotype simulated -------------------------------------------------------
  Ne <- 106  # Tsuda 2015 BMC genomics https://github.com/alenxav/Lectures/blob/master/Manuscripts/09_PGR_2018.pdf
  segSites  <-  as.integer(round(46000/20,0))  # Genome sequence of the palaeopolyploid soybean round(Nsnp/20*1.1,0)
  founderPop <- AlphaSimR::runMacs2(nInd = 20, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2612967/pdf/136.pdf
                                    nChr = 20,
                                    segSites = segSites,
                                    Ne = Ne,
                                    bp = 5.5e+07,#/(round(1115*10^6/20,0)), #Choi et al, 2007 TAG
                                    genLen = round(2550.3/20,0)/100, #Choi et al, 2007 TAG
                                    mutRate = 2.5e-08,#5.17*10^-3/1000000,  #Lavin, M., Herendeen, P. S. & Wojciechowski, M. F. Evolutionary rates
                                    histNe = c(500, 1500, 6000, 12000, 1e+05),  #T. E. CarterJr.T. HymowitzR. L. Nelson & Tracing soybean domestication history: From nucleotide to genome 2012
                                    histGen = c(100, 1000, 10000, 1e+05, 1e+06), #T. E. CarterJr.T. HymowitzR. L. Nelson & Tracing soybean domestication history: From nucleotide to genome 2012
                                    inbred = TRUE)
  SP <- SimParam$new(founderPop)
  SP$addSnpChip(nSnpPerChr=round(Nsnp/20,0))
  SP$addTraitA(nQtlPerChr=segSites*0.5, 
               mean = mgen, 
               var =  sigma2ge)
  SP$setVarE(h2=0.99999)
  pop <- newPop(founderPop, simParam=SP)
  pop <- randCross(pop, nCrosses=NF1, nProgeny = 1, simParam=SP)
  pop <- makeDH(pop, nDH=1, simParam=SP)
  pop <- selectCross(pop,nInd = NF1,nCrosses = NF1,use='pheno',simParam = SP)
  pop <- self(pop, nProgeny=NF2, simParam=SP)
  pop <- self(pop, nProgeny=1, simParam=SP)
  pop <- self(pop, nProgeny=1, simParam=SP)
  gen <- pullSnpGeno(pop, simParam = SP) 
  #plot(pop@gv,pop@pheno)
  genoS <- data.frame(gen=NA,
                      id=pop@id,
                      P1=pop@mother,
                      P2=pop@father,
                      mgen=mgen,
                      gv=pop@gv-mgen,
                      gen)
  row.names(genoS) <- NULL
  genoS <- genoS[sample(1:nrow(genoS),length(unique(data$gen))),]
  genoS$gen <- unique(data$gen)
  data <- merge(data,genoS,by.x='gen',by.y='gen',all.x=TRUE)
  data$xp <- data$x*Sp 
  data$yp <- data$y*Sr
  data$u_s <- book.rMatern(1,
                           data[,c('x','y')],
                           sigma=sqrt(sigma2u),
                           range=range,
                           kappa = sqrt(8*nu)/range,
                           variance = sigma2u,
                           nu=nu)
  data$u_e <- rnorm(nrow(data),0,sqrt(sigma2e))
  data$yield <- rowSums(data[,c('mgen','gv','u_s','u_e')])
  data <- cbind(data[,-grep('SNP',colnames(data))],data[,grep('SNP',colnames(data))])
  return(data)
}

insertPedM <- function (ped, founders = NULL){
  ped[, 1] <- as.character(ped[, 1])
  ped[, 2] <- as.character(ped[, 2])
  ped[, 3] <- as.character(ped[, 3])
  mmothers <- na.omit(ped[, 2][which(ped[, 2] %in% ped[, 1] == FALSE)])
  mfathers <- na.omit(ped[, 3][which(ped[, 3] %in% ped[, 1] == FALSE)])
  if (is.null(founders) == FALSE) {
    founders <- na.omit(founders[which(founders %in% ped[, 1] == FALSE)])
  }
  mparents <- unique(c(mmothers, mfathers, founders))
  nped <- ped[rep(1, length(mparents)), ]
  nped[, 1] <- mparents
  nped[, 2] <- NA
  nped[, 3] <- NA
  nped <- rbind(nped, ped)
  colnames(nped) <- colnames(ped)
  return(nped)
}

renum <- function(ped){## from packager ggroups
  colnames(ped) = c("ID", "SIRE", "DAM")
  for (i in 1:3) ped[, i] = paste0("x", ped[, i])
  ped[ped == "x0"] = 0
  newped = ped
  xrf = data.frame(ID = c(), newID = c())
  curr.set = ped[ped$SIRE == 0 & ped$DAM == 0, ]$ID
  xrf = data.frame(ID = curr.set, newID = 1:length(curr.set))
  ped = ped[!ped$ID %in% curr.set, ]
  gen = 1
  while (nrow(ped) > 0) {
    curr.set = ped[!ped$SIRE %in% ped$ID & !ped$DAM %in% ped$ID, ]
    curr.set = curr.set[order(curr.set$DAM, curr.set$SIRE), ]$ID
    xrf = rbind(xrf, data.frame(ID = curr.set, newID = (xrf[nrow(xrf), 
    ]$newID + 1):(xrf[nrow(xrf), ]$newID + length(curr.set))))
    ped = ped[!ped$ID %in% curr.set, ]
    gen = gen + 1
  }
  newped[] = xrf$newID[match(unlist(newped), xrf$ID)]
  newped[is.na(newped)] = 0
  newped = newped[order(newped$ID), ]
  xrf$ID = substring(xrf$ID, 2)
  #message("Found ", gen, " generations")
  return(list(newped = newped, xrf = xrf))
}

tabAinv <- function(ped, inbr){## from packager ggroups
  colnames(ped) = c("ID", "SIRE", "DAM")
  if (any(inbr < 0 | inbr > 1)) 
    stop("Inbreeding values should be between 0 and 1.")
  if (nrow(ped) != length(inbr)) 
    stop("Number of individuals in the pedigree does not match with the number of inbreeding values.")
  inbr = data.frame(ID = ped$ID, F = inbr)
  ped = merge(ped, inbr, by.x = "SIRE", by.y = "ID", 
              all.x = TRUE)
  ped = merge(ped, inbr, by.x = "DAM", by.y = "ID", 
              all.x = TRUE)
  ped = ped[order(ped$ID), ]
  colnames(ped)[4:5] = c("fs", "fd")
  ped[is.na(ped)] = 0
  Ai = data.frame(ID1 = ped$ID[1], ID2 = ped$ID[1], ai = 1)
  for (i in 2:nrow(ped)) {
    iID = ped$ID[i]
    si = ped$SIRE[i]
    di = ped$DAM[i]
    Fs = ped$fs[i]
    Fd = ped$fd[i]
    if (si > 0 & di == 0) {
      b = (3 - Fs)/4
      tmp = data.frame(ID1 = c(si, si, iID), ID2 = c(si, 
                                                     iID, iID), ai = c(0.25, -0.5, 1)/b)
      Ai = rbind(Ai, tmp)
    }
    if (si == 0 & di > 0) {
      b = (3 - Fd)/4
      tmp = data.frame(ID1 = c(di, di, iID), ID2 = c(di, 
                                                     iID, iID), ai = c(0.25, -0.5, 1)/b)
      Ai = rbind(Ai, tmp)
    }
    if (si > 0 & di > 0) {
      b = (2 - Fs - Fd)/4
      tmp = data.frame(ID1 = c(si, si, si, di, di, iID), 
                       ID2 = c(si, di, iID, di, iID, iID), ai = c(0.25, 
                                                                  0.25, -0.5, 0.25, -0.5, 1)/b)
      Ai = rbind(Ai, tmp)
    }
    if (si == 0 & di == 0) {
      Ai = rbind(Ai, c(iID, iID, 1))
    }
    if ((i%%1000) == 0) 
      message(i, " of ", nrow(ped))
  }
  if (i > 1000) 
    message(i, " of ", nrow(ped))
  Ai = transform(Ai, ID1 = pmin(Ai$ID1, Ai$ID2), ID2 = pmax(Ai$ID1, 
                                                            Ai$ID2))
  Ai = aggregate(Ai$ai, by = list(Ai$ID1, Ai$ID2), FUN = sum)
  colnames(Ai) = c("ID1", "ID2", "ai")
  return(Ai)
}
### Simulation and analysis ----------------------------------------------------
Spatial_simu <- function(i,...){
  cat(paste('Simulated: ',sprintf('%05i',i),Sys.time(),'\n'))
  Ncores = 1
  mgen <- 0
  ICCG = h2_s[i]
  ICCS = Spatial_per[i]
  VarG0 = ICCG
  VarE0 = 1-ICCG
  VarS0 = VarE0 * (1-ICCS)/ICCS
  VarY0 = VarG0 + VarS0 + VarE0
  sigma2ge <- (VarG0 / VarY0)*100
  sigma2e <- (VarE0 / VarY0)*100
  sigma2u <- (VarS0 / VarY0)*100
  range <- 5
  nu <- 1
  dat <- simulate_gen(r=r,
                      p=p,
                      f=f,
                      Nchecks = Nchecks,
                      sigma2e = sigma2e,
                      sigma2u = sigma2u,
                      range = range,
                      nu = nu,
                      mgen=mgen,
                      sigma2ge=sigma2ge,
                      unitD=F,
                      Sp=1,
                      Sr=5,
                      NF1=NF1,
                      NF2=NF2,
                      Nsnp=Nsnp,
                      Nrep=Nrep)
  dat <- dat[order(dat$y,dat$x),]
  Ma <- unique(dat[,grep('^gen$|SNP',colnames(dat))])
  GEN_label <- Ma$gen  
  
  ### Genotype matrix
  M <- Matrix::Matrix(as.matrix(Ma[,-1]-1),dimnames = list(GEN_label))## change to -1,0,1
  genoG <- GmatrixM(M)
  genoG_Invu <- solve(genoG)
  Cmatrix <- Matrix::Matrix(genoG_Invu,dimnames = list(GEN_label,GEN_label)) ### Genome relationship
  #image(Cmatrix,useRaster=TRUE)
  
  ### Pedigree matrix
  pedl <- unique(dat[,c('gen','P1','P2')])
  pedl <- insertPedM(pedl)
  row.names(pedl)  <- 1:nrow(pedl)
  pedl[is.na(pedl)] <- 0
  pedN <- renum(pedl)
  xx <- tabAinv(pedN$newped,rep(0,nrow(pedN$newped)))
  CmatrixP  <-  Matrix::sparseMatrix(i=xx[,1],j=xx[,2],x=xx[,3])
  dimnames(CmatrixP)[[1]] <- pedN$xrf$ID  
  dimnames(CmatrixP)[[2]] <- pedN$xrf$ID
  #image(CmatrixP)
  
  ### Independent matrix
  CmatrixI <- Matrix::Matrix(Diagonal(nrow(genoG)),dimnames = list(GEN_label,GEN_label)) # I independent 
  #image(CmatrixI)
  
  ### Spatial covariavel
  dat$block <- 1
  dat$Scov <- SPC(y = dat$yield,blk = dat$block,row = dat$x,col = dat$y, 2,2)
  
  ### Spatial Matrix covariable
  SpatialMatrix <- Matrix::Matrix(SPM(blk = dat$block,row = dat$x,col = dat$y,2,2))
  #image(SpatialMatrix)
  dat$id.Z <- 1:nrow(dat)
  Scalling_SM <- sum(apply(SpatialMatrix,2,var))
  
  ### Spatial kernel matrix
  D = as.matrix(stats::dist(dat[,c('x','y')])^2)
  K = exp(-D*0.25)
  diag(K) = diag(K)+0.01
  PrecK = Matrix::Matrix(solve(K))
  #image(PrecK,useRaster=TRUE)
  
  ### Storage results
  Time_started <- ymd_hms(Sys.time())
  results <- data.frame(LocationID = i,NPlots = length(na.omit(dat$yield)))
  results$mgen <- mgen
  results$sigma2_ge <- sigma2ge
  results$sigma2_e <- sigma2e
  results$H2_true <- sigma2ge/(sigma2ge+sigma2e)
  results$sigma2_u <- sigma2u
  results$kappa <- sqrt(8*nu)/range
  results$range <- range
  results$Spatial_per <- Spatial_per[i]
  lgprior <- list(prec = list(prior="loggamma", param = c(1,5e-05))) ## Defalt for generic0
  
  ##  No Spatial INLA  ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2ns <- yield ~ f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  m2ns <- inla(f2ns,data=dat,verbose=F,
               only.hyperparam = F,
               control.compute=list(dic=F,config=TRUE),
               control.predictor=list(compute = TRUE),
               num.threads = Ncores)
  #summary(m2ns)
  BLUPns <- m2ns$summary.random$gen[,1:3]
  colnames(BLUPns)[2:3] <- c('BLUPns','sdns')
  dat <- merge(dat,BLUPns[,1:3],by.x='gen',by.y='ID')
  results$Var_nsinla_erro <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2ns$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_nsinla_ge <- inla.emarginal(function(x) x,
                                          inla.tmarginal(function(x) 1/exp(x),
                                                         m2ns$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_nsinla <-  results$Var_nsinla_ge / (results$Var_nsinla_ge + results$Var_nsinla_erro)
  results$loglik_nsinla <-   m2ns$mlik[2]
  #print(results)
  
  ##  No Spatial INLA  GENOTYPE ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2nsG <-  yield ~ f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  m2nsG <- inla(f2nsG,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2nsG)
  BLUPnsG <- m2nsG$summary.random$gen[,1:3]
  colnames(BLUPnsG)[2:3] <- c('BLUPnsG','sdnsG')
  dat <- merge(dat,BLUPnsG[,1:3],by.x='gen',by.y='ID')
  results$Var_nsinlaG_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2nsG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_nsinlaG_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2nsG$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_nsinlaG <-  results$Var_nsinlaG_ge / (results$Var_nsinlaG_ge + results$Var_nsinlaG_erro)
  results$loglik_nsinlaG <-   m2nsG$mlik[2]
  #print(results)
  
  ##  No Spatial INLA  PEDIGREE ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2nsP <-  yield ~ f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  m2nsP <- inla(f2nsP,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2nsP)
  BLUPnsP <- m2nsP$summary.random$gen[,1:3]
  colnames(BLUPnsP)[2:3] <- c('BLUPnsP','sdnsP')
  dat <- merge(dat,BLUPnsP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_nsinlaP_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2nsP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_nsinlaP_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2nsP$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_nsinlaP <-  results$Var_nsinlaP_ge / (results$Var_nsinlaP_ge + results$Var_nsinlaP_erro)
  results$loglik_nsinlaP <-   m2nsP$mlik[2]
  #print(results)
  
  ##  Row + Col INla --------------------------------------------------------------------------
  ## function SPC packager NAM BWGR (models)
  dat <- dat[order(dat$y,dat$x),]
  f2rc <- yield ~ 
    f(x,model="iid",constr = T,hyper = lgprior) +
    f(y,model="iid",constr = T,hyper = lgprior) + 
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  
  m2rc <- inla(f2rc,data=dat,verbose=F,
               only.hyperparam = F,
               control.compute=list(dic=F,config=TRUE),
               control.predictor=list(compute = TRUE),
               num.thread= Ncores)
  #summary(m2rc)
  BLUPrc <- m2rc$summary.random$gen[,1:3]
  colnames(BLUPrc)[2:3] <- c('BLUPrc','sdrc')
  dat <- merge(dat,BLUPrc[,1:3],by.x='gen',by.y='ID')
  xef <- m2rc$summary.random$x[,1:3]
  colnames(xef) <- c('x','x_rc','xsd_rc')
  yef <- m2rc$summary.random$y[,1:3]
  colnames(yef) <- c('y','y_rc','ysd_rc')
  dat <- merge(dat,xef,by.x='x',by.y='x')
  dat <- merge(dat,yef,by.x='y',by.y='y')
  dat$u_spatial_rc <- rowSums(dat[,c('y_rc','x_rc')])
  dat$usd_spatial_rc <- rowSums(dat[,c('ysd_rc','xsd_rc')])
  results$Var_rcinla_erro <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2rc$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_rcinla_ge <- inla.emarginal(function(x) x,
                                          inla.tmarginal(function(x) 1/exp(x),
                                                         m2rc$internal.marginals.hyperpar[["Log precision for gen"]]))
  
  results$Var_rcinla_spatial <- inla.emarginal(function(x) x,
                                               inla.tmarginal(function(x) 1/exp(x),
                                                              m2rc$internal.marginals.hyperpar[["Log precision for x"]])) + 
    inla.emarginal(function(x) x,
                   inla.tmarginal(function(x) 1/exp(x),
                                  m2rc$internal.marginals.hyperpar[["Log precision for y"]]))
  
  results$H2_rcinla <-  results$Var_rcinla_ge / (results$Var_rcinla_ge + results$Var_rcinla_erro)
  results$loglik_rcinla <-   m2rc$mlik[2]
  #print(results)
  
  ##  Row + Col INla GENOTYPE --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2rcG <- yield ~ 
    f(x,model="iid",constr = T,hyper = lgprior) +
    f(y,model="iid",constr = T,hyper = lgprior) + 
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2rcG <- inla(f2rcG,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.thread= Ncores)
  #summary(m2rc)
  BLUPrcG <- m2rcG$summary.random$gen[,1:3]
  colnames(BLUPrcG)[2:3] <- c('BLUPrcG','sdrcG')
  dat <- merge(dat,BLUPrcG[,1:3],by.x='gen',by.y='ID')
  
  xef <- m2rcG$summary.random$x[,1:3]
  colnames(xef) <- c('x','x_rcG','xsd_rcG')
  yef <- m2rcG$summary.random$y[,1:3]
  colnames(yef) <- c('y','y_rcG','ysd_rcG')
  dat <- merge(dat,xef,by.x='x',by.y='x')
  dat <- merge(dat,yef,by.x='y',by.y='y')
  dat$u_spatial_rcG <- rowSums(dat[,c('y_rcG','x_rcG')])
  dat$usd_spatial_rcG <- rowSums(dat[,c('ysd_rcG','xsd_rcG')])
  results$Var_rcinlaG_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2rcG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_rcinlaG_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2rcG$internal.marginals.hyperpar[["Log precision for gen"]]))
  
  results$Var_rcinlaG_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2rcG$internal.marginals.hyperpar[["Log precision for x"]])) + 
    inla.emarginal(function(x) x,
                   inla.tmarginal(function(x) 1/exp(x),
                                  m2rcG$internal.marginals.hyperpar[["Log precision for y"]]))
  
  results$H2_rcinlaG <-  results$Var_rcinlaG_ge / (results$Var_rcinlaG_ge + results$Var_rcinlaG_erro)
  results$loglik_rcinlaG <-   m2rcG$mlik[2]
  #print(results)
  
  ##  Row + Col INla PEDIGRE --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2rcP <- yield ~ 
    f(x,model="iid",constr = T,hyper = lgprior) +
    f(y,model="iid",constr = T,hyper = lgprior) + 
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2rcP <- inla(f2rcP,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.thread= Ncores)
  #summary(m2rcP)
  BLUPrcP <- m2rcP$summary.random$gen[,1:3]
  colnames(BLUPrcP)[2:3] <- c('BLUPrcP','sdrcP')
  dat <- merge(dat,BLUPrcP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  
  xef <- m2rcP$summary.random$x[,1:3]
  colnames(xef) <- c('x','x_rcP','xsd_rcP')
  yef <- m2rcP$summary.random$y[,1:3]
  colnames(yef) <- c('y','y_rcP','ysd_rcP')
  dat <- merge(dat,xef,by.x='x',by.y='x')
  dat <- merge(dat,yef,by.x='y',by.y='y')
  dat$u_spatial_rcP <- rowSums(dat[,c('y_rcP','x_rcP')])
  dat$usd_spatial_rcP <- rowSums(dat[,c('ysd_rcP','xsd_rcP')])
  results$Var_rcinlaP_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2rcP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_rcinlaP_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2rcP$internal.marginals.hyperpar[["Log precision for gen"]]))
  
  results$Var_rcinlaP_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2rcP$internal.marginals.hyperpar[["Log precision for x"]])) + 
    inla.emarginal(function(x) x,
                   inla.tmarginal(function(x) 1/exp(x),
                                  m2rcP$internal.marginals.hyperpar[["Log precision for y"]]))
  
  results$H2_rcinlaP <-  results$Var_rcinlaP_ge / (results$Var_rcinlaP_ge + results$Var_rcinlaP_erro)
  results$loglik_rcinlaP <-   m2rcP$mlik[2]
  #print(results)
  ## Spatial Covariate  ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2sc <- yield ~ 
    f(Scov,model="linear")+
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  
  m2sc <- inla(f2sc,data=dat,verbose=F,
               only.hyperparam = F,
               control.compute=list(dic=F,config=TRUE),
               control.predictor=list(compute = TRUE),
               num.threads = Ncores)
  #summary(m2sc)
  BLUPsc <- m2sc$summary.random$gen[,1:3]
  colnames(BLUPsc)[2:3] <- c('BLUPsc','sdsc')
  dat$u_spatial_sc <- dat$Scov
  dat$usd_spatial_sc <- m2sc$summary.fixed$sd[2] ### ver sd!!!!!
  dat <- merge(dat,BLUPsc[,1:3],by.x='gen',by.y='ID')
  results$Var_scinla_erro <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2sc$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_scinla_ge <- inla.emarginal(function(x) x,
                                          inla.tmarginal(function(x) 1/exp(x),
                                                         m2sc$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_scinla <-  results$Var_scinla_ge / (results$Var_scinla_ge + results$Var_scinla_erro)
  results$loglik_scinla <-  m2sc$mlik[2]
  #print(results)
  ## Spatial Covariate Genotype ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2scG <- yield ~ 
    f(Scov,model="linear")+
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2scG <- inla(f2scG,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2scG)
  BLUPscG <- m2scG$summary.random$gen[,1:3]
  colnames(BLUPscG)[2:3] <- c('BLUPscG','sdscG')
  dat$u_spatial_scG <- dat$Scov
  dat$usd_spatial_scG <- m2scG$summary.fixed$sd[2] ### ver sd!!!!!
  dat <- merge(dat,BLUPscG[,1:3],by.x='gen',by.y='ID')
  results$Var_scGinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2scG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_scGinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2scG$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_scGinla <-  results$Var_scGinla_ge / (results$Var_scGinla_ge + results$Var_scGinla_erro)
  results$loglik_scGinla <-  m2scG$mlik[2]
  #print(results)
  ## Spatial Covariate Pedigree ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2scP <- yield ~ 
    f(Scov,model="linear")+
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2scP <- inla(f2scP,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2scP)
  BLUPscP <- m2scP$summary.random$gen[,1:3]
  colnames(BLUPscP)[2:3] <- c('BLUPscP','sdscP')
  dat$u_spatial_scP <- dat$Scov
  dat$usd_spatial_scP <- m2scP$summary.fixed$sd[2] ### ver sd!!!!!
  dat <- merge(dat,BLUPscP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_scPinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2scP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_scPinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2scP$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_scPinla <-  results$Var_scPinla_ge / (results$Var_scPinla_ge + results$Var_scPinla_erro)
  results$loglik_scPinla <-  m2scP$mlik[2]
  #print(results)
  ## Spatial Matrix Covariate ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2sm  <- yield ~ 
    f(id.Z,Z=SpatialMatrix,model="z")+
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  m2sm <- inla(f2sm,data=dat,verbose=F,
               only.hyperparam = F,
               control.compute=list(dic=F,config=TRUE),
               control.predictor=list(compute = TRUE),
               num.threads = Ncores)
  #summary(m2sm)
  xef <- m2sm$summary.random$id.Z[-(1:nrow(dat)),1:3]
  colnames(xef) <- c('x','u_sm','usd_sm')
  dat$u_spatial_sm <- as.vector(unlist(SpatialMatrix %*% m2sm$summary.random$id.Z[-(1:nrow(dat)),2]))
  dat$usd_spatial_sm <- m2sm$summary.random$id.Z[-(1:nrow(dat)),3]
  BLUPsm <- m2sm$summary.random$gen[,1:3]
  colnames(BLUPsm)[2:3] <- c('BLUPsm','sdsm')
  dat <- merge(dat,BLUPsm[,1:3],by.x='gen',by.y='ID')
  results$Var_sminla_erro <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2sm$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_sminla_ge <- inla.emarginal(function(x) x,
                                          inla.tmarginal(function(x) 1/exp(x),
                                                         m2sm$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_sminla_spatial <- inla.emarginal(function(x) x,
                                               inla.tmarginal(function(x) 1/exp(x),
                                                              m2sm$internal.marginals.hyperpar[["Log precision for id.Z"]]))*Scalling_SM
  results$H2_sminla <-  results$Var_sminla_ge / (results$Var_sminla_ge + results$Var_sminla_erro)
  results$loglik_sminla <-  m2sm$mlik[2]
  #print(results)
  ## Spatial Matrix Covariate Genotype----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2smG  <- yield ~ 
    f(id.Z,Z=SpatialMatrix,model="z")+
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2smG <- inla(f2smG,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2smG)
  xef <- m2smG$summary.random$id.Z[-(1:nrow(dat)),1:3]
  colnames(xef) <- c('x','u_smG','usd_smG')
  dat$u_spatial_smG <- as.vector(unlist(SpatialMatrix %*% m2smG$summary.random$id.Z[-(1:nrow(dat)),2]))
  dat$usd_spatial_smG <- m2smG$summary.random$id.Z[-(1:nrow(dat)),3]
  BLUPsmG <- m2smG$summary.random$gen[,1:3]
  colnames(BLUPsmG)[2:3] <- c('BLUPsmG','sdsmG')
  dat <- merge(dat,BLUPsmG[,1:3],by.x='gen',by.y='ID')
  results$Var_smGinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2smG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_smGinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2smG$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_smGinla_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2smG$internal.marginals.hyperpar[["Log precision for id.Z"]]))*Scalling_SM
  results$H2_smGinla <-  results$Var_smGinla_ge / (results$Var_smGinla_ge + results$Var_smGinla_erro)
  results$loglik_smGinla <-  m2smG$mlik[2]
  #print(results)
  ## Spatial Matrix Covariate Pedigreee ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2smP  <- yield ~ 
    f(id.Z,Z=SpatialMatrix,model="z")+
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2smP <- inla(f2smP,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2smP)
  xef <- m2smP$summary.random$id.Z[-(1:nrow(dat)),1:3]
  colnames(xef) <- c('x','u_smP','usd_smP')
  dat$u_spatial_smP <- as.vector(unlist(SpatialMatrix %*% m2smP$summary.random$id.Z[-(1:nrow(dat)),2]))
  dat$usd_spatial_smP <- m2smP$summary.random$id.Z[-(1:nrow(dat)),3]
  BLUPsmP <- m2smP$summary.random$gen[,1:3]
  colnames(BLUPsmP)[2:3] <- c('BLUPsmP','sdsmP')
  dat <- merge(dat,BLUPsmP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_smPinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2smP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_smPinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2smP$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_smPinla_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2smP$internal.marginals.hyperpar[["Log precision for id.Z"]]))*Scalling_SM
  results$H2_smPinla <-  results$Var_smPinla_ge / (results$Var_smPinla_ge + results$Var_smPinla_erro)
  results$loglik_smPinla <-  m2smP$mlik[2]
  #print(results)
  ### Gaussian Kernel --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2gk  <- yield ~ f(id.Z,model="generic0",Cmatrix=PrecK,constr=T, hyper = lgprior)+
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  
  m2gk <- inla(f2gk,data=dat,verbose=F,
               only.hyperparam = F,
               control.compute=list(dic=F,config=TRUE),
               control.predictor=list(compute = TRUE),
               num.threads = Ncores)
  #summary(m2gk)
  #dat$u_spatial_gk <- as.vector(unlist(PrecK %*% m2gk$summary.random$id.Z[,2]))
  dat$u_spatial_gk <- m2gk$summary.random$id.Z[,2]
  dat$usd_spatial_gk <- m2gk$summary.random$id.Z[,3]
  BLUPgk <- m2gk$summary.random$gen[,1:3]
  colnames(BLUPgk)[2:3] <- c('BLUPgk','sdgk')
  dat <- merge(dat,BLUPgk[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_gkinla_erro <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2gk$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_gkinla_ge <- inla.emarginal(function(x) x,
                                          inla.tmarginal(function(x) 1/exp(x),
                                                         m2gk$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_gkinla_spatial <- inla.emarginal(function(x) x,
                                               inla.tmarginal(function(x) 1/exp(x),
                                                              m2gk$internal.marginals.hyperpar[["Log precision for id.Z"]]))
  results$H2_gkinla <-  results$Var_gkinla_ge / (results$Var_gkinla_ge + results$Var_gkinla_erro)
  results$loglik_gkinla <-  m2gk$mlik[2]
  #print(results)
  
  ### Gaussian Kernel Genotype --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2gkG  <- yield ~ 
    f(id.Z,model="generic0",Cmatrix=PrecK,constr=T, hyper = lgprior)+
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2gkG <- inla(f2gkG,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2gkG)
  #dat$u_spatial_gkG <- as.vector(unlist(PrecK %*% m2gkG$summary.random$id.Z[,2]))
  dat$u_spatial_gkG <- m2gkG$summary.random$id.Z[,2]
  dat$usd_spatial_gkG <- m2gkG$summary.random$id.Z[,3]
  BLUPgkG <- m2gkG$summary.random$gen[,1:3]
  colnames(BLUPgkG)[2:3] <- c('BLUPgkG','sdgkG')
  dat <- merge(dat,BLUPgkG[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_gkGinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2gkG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_gkGinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2gkG$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_gkGinla_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2gkG$internal.marginals.hyperpar[["Log precision for id.Z"]]))
  results$H2_gkGinla <-  results$Var_gkGinla_ge / (results$Var_gkGinla_ge + results$Var_gkGinla_erro)
  results$loglik_gkGinla <-  m2gkG$mlik[2]
  #print(results)
  
  ### Gaussian Kernel Pedigree --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  f2gkP  <- yield ~ 
    f(id.Z,model="generic0",Cmatrix=PrecK,constr=T, hyper = lgprior)+
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2gkP <- inla(f2gkP,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2gkP)
  #dat$u_spatial_gkP <- as.vector(unlist(PrecK %*% m2gkP$summary.random$id.Z[,2]))
  dat$u_spatial_gkP <- m2gkP$summary.random$id.Z[,2]
  dat$usd_spatial_gkP <- m2gkP$summary.random$id.Z[,3]
  BLUPgkP <- m2gkP$summary.random$gen[,1:3]
  colnames(BLUPgkP)[2:3] <- c('BLUPgkP','sdgkP')
  dat <- merge(dat,BLUPgkP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_gkPinla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2gkP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_gkPinla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2gkP$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_gkPinla_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2gkP$internal.marginals.hyperpar[["Log precision for id.Z"]]))
  results$H2_gkPinla <-  results$Var_gkPinla_ge / (results$Var_gkPinla_ge + results$Var_gkPinla_erro)
  results$loglik_gkPinla <-  m2gkP$mlik[2]
  #print(results)
  ## AR1 x AR1 INLA  ----------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  dat$x <- as.numeric(dat$x)
  dat$y <- as.numeric(dat$y)
  f2ar1 <- yield ~ 
    f(x,model="ar1",group = y,control.group = list(model="ar1"))+
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  
  m2ar1 <- inla(f2ar1,data=dat,verbose=F,
                only.hyperparam = F,
                control.compute=list(dic=F,config=TRUE),
                control.predictor=list(compute = TRUE),
                num.threads = Ncores)
  #summary(m2ar1)
  xef <- m2ar1$summary.random$x[,1:3]
  colnames(xef) <- c('x','u_ar1','usd_ar1')
  dat$u_spatial_ar1 <- xef$u_ar1
  dat$usd_spatial_ar1 <- xef$usd_ar1
  BLUPar1 <- m2ar1$summary.random$gen[,1:3]
  colnames(BLUPar1)[2:3] <- c('BLUPar1','sdar1')
  dat <- merge(dat,BLUPar1[,1:3],by.x='gen',by.y='ID')
  results$Var_ar1inla_erro <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2ar1$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_ar1inla_ge <- inla.emarginal(function(x) x,
                                           inla.tmarginal(function(x) 1/exp(x),
                                                          m2ar1$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_ar1inla_spatial <- inla.emarginal(function(x) x,
                                                inla.tmarginal(function(x) 1/exp(x),
                                                               m2ar1$internal.marginals.hyperpar[["Log precision for x"]]))
  results$H2_ar1inla <-  results$Var_ar1inla_ge / (results$Var_ar1inla_ge + results$Var_ar1inla_erro)
  results$loglik_ar1inla <-  m2ar1$mlik[2]
  #print(results)
  ## AR1 x AR1 INLA  GENOTYPE --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  dat$x <- as.numeric(dat$x)
  dat$y <- as.numeric(dat$y)
  f2ar1G <- yield ~ 
    f(x,model="ar1",group = y,control.group = list(model="ar1"))+
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2ar1G <- inla(f2ar1G,data=dat,verbose=F,
                 only.hyperparam = F,
                 control.compute=list(dic=F,config=TRUE),
                 control.predictor=list(compute = TRUE),
                 num.threads = Ncores)
  #summary(m2ar1G)
  xef <- m2ar1G$summary.random$x[,1:3]
  colnames(xef) <- c('x','u_ar1G','usd_ar1G')
  dat$u_spatial_ar1G <- xef$u_ar1G
  dat$usd_spatial_ar1G <- xef$usd_ar1G
  
  BLUPar1G <- m2ar1G$summary.random$gen[,1:3]
  colnames(BLUPar1G)[2:3] <- c('BLUPar1G','sdar1G')
  dat <- merge(dat,BLUPar1G[,1:3],by.x='gen',by.y='ID')
  results$Var_ar1inlaG_erro <- inla.emarginal(function(x) x,
                                              inla.tmarginal(function(x) 1/exp(x),
                                                             m2ar1G$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_ar1inlaG_ge <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2ar1G$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_ar1inlaG_spatial <- inla.emarginal(function(x) x,
                                                 inla.tmarginal(function(x) 1/exp(x),
                                                                m2ar1G$internal.marginals.hyperpar[["Log precision for x"]]))
  results$H2_ar1inlaG <-  results$Var_ar1inlaG_ge / (results$Var_ar1inlaG_ge + results$Var_ar1inlaG_erro)
  results$loglik_ar1inlaG <-  m2ar1G$mlik[2]
  #print(results)
  ## AR1 x AR1 INLA  PEDIGRE --------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  dat$x <- as.numeric(dat$x)
  dat$y <- as.numeric(dat$y)
  f2ar1P <- yield ~ 
    f(x,model="ar1",group = y,control.group = list(model="ar1"))+
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2ar1P <- inla(f2ar1P,data=dat,verbose=F,
                 only.hyperparam = F,
                 control.compute=list(dic=F,config=TRUE),
                 control.predictor=list(compute = TRUE),
                 num.threads = Ncores)
  #summary(m2ar1P)
  xef <- m2ar1P$summary.random$x[,1:3]
  colnames(xef) <- c('x','u_ar1P','usd_ar1P')
  dat$u_spatial_ar1P <- xef$u_ar1P
  dat$usd_spatial_ar1P <- xef$usd_ar1P
  
  BLUPar1P <- m2ar1P$summary.random$gen[,1:3]
  colnames(BLUPar1P)[2:3] <- c('BLUPar1P','sdar1P')
  dat <- merge(dat,BLUPar1P[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  results$Var_ar1inlaP_erro <- inla.emarginal(function(x) x,
                                              inla.tmarginal(function(x) 1/exp(x),
                                                             m2ar1P$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_ar1inlaP_ge <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2ar1P$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$Var_ar1inlaP_spatial <- inla.emarginal(function(x) x,
                                                 inla.tmarginal(function(x) 1/exp(x),
                                                                m2ar1P$internal.marginals.hyperpar[["Log precision for x"]]))
  results$H2_ar1inlaP <-  results$Var_ar1inlaP_ge / (results$Var_ar1inlaP_ge + results$Var_ar1inlaP_erro)
  results$loglik_ar1inlaP <-  m2ar1P$mlik[2]
  #print(results)
  
  ### MESH and STACK  ---------------------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  xrange <- range(dat$x,na.rm=T)
  yrange <- range(dat$y,na.rm=T)
  bordas <- data.frame(x=c(xrange[1],xrange[2],xrange[2],xrange[1],xrange[1]),
                       y=c(yrange[1],yrange[1],yrange[2],yrange[2],yrange[1]))
  mp <- max(max(dat$x),max(dat$y))*0.05
  mesh <- inla.mesh.2d(loc=dat[,c('x','y')],
                       loc.domain=bordas,
                       offset=mp*0.1,
                       max.edge=mp,
                       cutoff=mp)
  spde <- inla.spde2.matern(mesh=mesh,alpha=2)
  A <- inla.spde.make.A(mesh=mesh,loc=as.matrix(dat[,c('x','y')]))
  stk.e <- inla.stack(tag='est',data=list(yield=dat$yield),
                      A=list(A,1,1,1),
                      effects=list(field=1:spde$n.spde,
                                   gen=as.factor(dat$gen),
                                   xp=dat$xp,
                                   yp=dat$yp))
  #plot(mesh,col='gray10');lines(bordas,lwd=3);points(dat[,c('x','y')],col='red',cex=1,pch='+')
  ### Model INLA SPDE -----------------------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  fspde <- yield ~ 
    f(field,model=spde)+
    f(gen,model="generic0",Cmatrix=CmatrixI,constr=T,values=GEN_label, hyper = lgprior)
  
  m2spde <- inla(fspde,data=inla.stack.data(stk.e),verbose=F,
                 only.hyperparam = F,
                 control.predictor=list(A=inla.stack.A(stk.e),compute = TRUE),
                 control.compute=list(dic=F,config=TRUE),
                 num.threads = Ncores)
  #summary(m2spde)
  BLUPspde <- m2spde$summary.random$gen[,1:3]
  colnames(BLUPspde)[2:3] <- c('BLUPspde','sdspde')
  dat <- merge(dat,BLUPspde[,1:3],by.x='gen',by.y='ID')
  
  proj_grid <- inla.mesh.projector(mesh,loc = as.matrix(dat[,c('x','y')]))
  dat$u_spatial_spde <- inla.mesh.project(proj_grid,m2spde$summary.random$field$mean)
  dat$usd_spatial_spde <- inla.mesh.project(proj_grid,m2spde$summary.random$field$sd)
  #plot(dat$u_s,dat$u_spatial_spde)
  results$Var_spdeinla_erro <- inla.emarginal(function(x) x,
                                              inla.tmarginal(function(x) 1/exp(x),
                                                             m2spde$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_spdeinla_ge <- inla.emarginal(function(x) x,
                                            inla.tmarginal(function(x) 1/exp(x),
                                                           m2spde$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_spdeinla <-  results$Var_spdeinla_ge / (results$Var_spdeinla_ge + results$Var_spdeinla_erro)
  spde.result <- inla.spde2.result(m2spde, "field", spde, do.transform=TRUE)
  results$Var_spdeinla_spatial <- inla.emarginal(function(x) x, 
                                                 spde.result$marginals.variance.nominal[["variance.nominal.1"]])
  results$Range_spdeinla <- inla.emarginal(function(x) x, 
                                           spde.result$marginals.range.nominal[["range.nominal.1"]])
  results$kappa <- sqrt(8 * 1)/results$Range_spdeinla
  results$loglik_spdeinla <- m2spde$mlik[2]
  #print(results)
  ### Model INLA SPDE GENOTYPE -------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  fspdeG <- yield ~ 
    f(field,model=spde)+
    f(gen,model="generic0",Cmatrix=Cmatrix,constr=T,values=GEN_label, hyper = lgprior)
  
  m2spdeG <- inla(fspdeG,data=inla.stack.data(stk.e),verbose=F,
                  only.hyperparam = F,
                  control.predictor=list(A=inla.stack.A(stk.e),compute = TRUE),
                  control.compute=list(dic=F,config=TRUE),
                  num.threads = Ncores)
  #summary(m2spde)
  BLUPspdeG <- m2spdeG$summary.random$gen[,1:3]
  colnames(BLUPspdeG)[2:3] <- c('BLUPspdeG','sdspdeG')
  dat <- merge(dat,BLUPspdeG[,1:3],by.x='gen',by.y='ID')
  dat$u_spatial_spdeG <- inla.mesh.project(proj_grid,m2spdeG$summary.random$field$mean)
  dat$usd_spatial_spdeG <- inla.mesh.project(proj_grid,m2spdeG$summary.random$field$sd)
  #plot(dat$u_s,dat$u_spatial_spde)
  results$Var_spdeinlaG_erro <- inla.emarginal(function(x) x,
                                               inla.tmarginal(function(x) 1/exp(x),
                                                              m2spdeG$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_spdeinlaG_ge <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2spdeG$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_spdeinlaG <-  results$Var_spdeinlaG_ge / (results$Var_spdeinlaG_ge + results$Var_spdeinlaG_erro)
  
  spde.resultG <- inla.spde2.result(m2spdeG, "field", spde, do.transform=TRUE)
  results$Var_spdeinlaG_spatial <- inla.emarginal(function(x) x, 
                                                  spde.resultG$marginals.variance.nominal[["variance.nominal.1"]])
  results$Range_spdeinlaG <- inla.emarginal(function(x) x, 
                                            spde.resultG$marginals.range.nominal[["range.nominal.1"]])
  results$kappaG <- sqrt(8 * 1)/results$Range_spdeinlaG
  results$loglik_spdeinlaG <- m2spdeG$mlik[2]
  #print(results)
  ### Model INLA SPDE Pedigree  -------------------------------------------------
  dat <- dat[order(dat$y,dat$x),]
  fspdeP <- yield ~ 
    f(field,model=spde)+
    f(gen,model="generic0",Cmatrix=CmatrixP,values=pedN$xrf$ID,constr=T, hyper = lgprior)
  
  m2spdeP <- inla(fspdeP,data=inla.stack.data(stk.e),verbose=F,
                  only.hyperparam = F,
                  control.predictor=list(A=inla.stack.A(stk.e),compute = TRUE),
                  control.compute=list(dic=F,config=TRUE),
                  num.threads = Ncores)
  #summary(m2spdeP)
  BLUPspdeP <- m2spdeP$summary.random$gen[,1:3]
  colnames(BLUPspdeP)[2:3] <- c('BLUPspdeP','sdspdeP')
  dat <- merge(dat,BLUPspdeP[,1:3],by.x='gen',by.y='ID',all.x=TRUE)
  dat$u_spatial_spdeP <- inla.mesh.project(proj_grid,m2spdeP$summary.random$field$mean)
  dat$usd_spatial_spdeP <- inla.mesh.project(proj_grid,m2spdeP$summary.random$field$sd)
  #plot(dat$u_s,dat$u_spatial_spdeP)
  results$Var_spdeinlaP_erro <- inla.emarginal(function(x) x,
                                               inla.tmarginal(function(x) 1/exp(x),
                                                              m2spdeP$internal.marginals.hyperpar[["Log precision for the Gaussian observations"]]))
  results$Var_spdeinlaP_ge <- inla.emarginal(function(x) x,
                                             inla.tmarginal(function(x) 1/exp(x),
                                                            m2spdeP$internal.marginals.hyperpar[["Log precision for gen"]]))
  results$H2_spdeinlaP <-  results$Var_spdeinlaP_ge / (results$Var_spdeinlaP_ge + results$Var_spdeinlaP_erro)
  
  spde.resultP <- inla.spde2.result(m2spdeP, "field", spde, do.transform=TRUE)
  results$Var_spdeinlaP_spatial <- inla.emarginal(function(x) x, 
                                                  spde.resultP$marginals.variance.nominal[["variance.nominal.1"]])
  results$Range_spdeinlaP <- inla.emarginal(function(x) x, 
                                            spde.resultP$marginals.range.nominal[["range.nominal.1"]])
  results$kappaP <- sqrt(8 * 1)/results$Range_spdeinlaP
  results$loglik_spdeinlaP <- m2spdeP$mlik[2]
  #print(results)
  
  ### Rank ---------------------------------------------------------------------------
  dat$Rank <- 0
  dat <- dat[order(dat$gv,decreasing = T),]
  dat$Rank[1:Nsel] <- 1
  ### output values -----------------------------------------------------------------------------
  #Correlation ----------------------
  w = which(!duplicated(dat$gen))
  datS <- dat[w,grep("^gv|^BLUP",colnames(dat))]
  Corr <- cor(datS,method = "pearson")[1,-1]
  names(Corr) <- paste0('Corr_gv_',names(Corr))
  results <- data.frame(results,data.frame(t(Corr)))
  
  datS <- dat[,grep("^u_s|^u_spatial_",colnames(dat))]
  Corr <- cor(datS,method = "pearson")[1,-1]
  names(Corr) <- paste0('Corr_us_',gsub('u_spatial_','',names(datS[,-1])))
  results <- data.frame(results,data.frame(t(Corr)))
  
  ## CRPS----------------------------
  datS <- dat[,grep("^u_s|^u_spatial_|usd_spatial",colnames(dat))]
  models <- unique(gsub("u_s|u_spatial_|usd_spatial_","",colnames(datS)))
  models <-   models[nchar(models)>1]
  crps_r <- c()
  for(m in models){
    crps_r <- c(crps_r,crps(datS$u_s,dat[,c(paste0('u_spatial_',m),paste0('usd_spatial_',m))])$CRPS)
  }
  names(crps_r) <- paste0('CRPS_us_',models)
  results <- data.frame(results,data.frame(t(crps_r)))
  
  datS <- dat[,grep("^gv|^BLUP|^sd",colnames(dat))]
  models <- unique(gsub("gv|^BLUP|^sd","",colnames(datS)))
  models <-   models[nchar(models)>1]
  crps_r <- c()
  for(m in models){
    crps_r <- c(crps_r,crps(datS$gv,dat[,c(paste0('BLUP',m),paste0('sd',m))])$CRPS)
  }
  names(crps_r) <- paste0('CRPS_gv_',models)
  results <- data.frame(results,data.frame(t(crps_r)))
  
  datS <- dat[,grep("^BLUP|Rank",colnames(dat))]
  models <- unique(gsub("Rank|^BLUP","",colnames(datS)))
  models <-   models[nchar(models)>1]
  crps_r <- c()
  for(m in models){
    crps_r <- c(crps_r,sum(datS[order(datS[,paste0('BLUP',m)],decreasing = T),'Rank'][1:Nsel])/Nsel)
  }
  names(crps_r) <- paste0('Nsel_BLUP_',models)
  results <- data.frame(results,data.frame(t(crps_r)))
  
  # Bias ----------------------------
  datS <- results[,grep("^H2_",colnames(results))]
  difft <- as.numeric(datS) - as.numeric(datS['H2_true'])
  names(difft) <- paste0('Bias_',gsub('inla','',colnames(datS)))
  results <- data.frame(results,data.frame(t(difft[difft!=0])))
  
  datS <- results[,grep("_ge$|sigma2_ge",colnames(results))]
  difft <- as.numeric(datS) - as.numeric(datS['sigma2_ge'])
  names(difft) <- paste0('Bias_sigma2_ge_',gsub('Var_|_ge|inla','',colnames(datS)))
  results <- data.frame(results,data.frame(t(difft[difft!=0])))
  
  datS <- results[,grep("_erro$|sigma2_e",colnames(results))]
  difft <- as.numeric(datS) - as.numeric(datS['sigma2_e'])
  names(difft) <- paste0('Bias_sigma2_erro_',gsub('Var_|_erro|inla','',colnames(datS)))
  results <- data.frame(results,data.frame(t(difft[difft!=0])))
  
  datS <- results[,grep("_spatial$|sigma2_u",colnames(results))]
  difft <- as.numeric(datS) - as.numeric(datS['sigma2_u'])
  names(difft) <- paste0('Bias_sigma2_u_',gsub('Var_|_spatial|inla','',colnames(datS)))
  results <- data.frame(results,data.frame(t(difft[difft!=0])))
  
  datS <- results[,grep("^loglik|loglik_nsinla",colnames(results))]
  difft <- as.numeric(datS) / as.numeric(datS['loglik_nsinla'])
  names(difft) <- paste0('NOR',gsub('inla','',colnames(datS)))
  results <- data.frame(results,data.frame(t(difft[difft!=1])))
  results$Time <- as.numeric(as.duration(interval(Time_started,ymd_hms(Sys.time()))),"secs")
  #print(results)
  return(results)
}
########################################################################## end..