### Analysis Spatial variation -------------------------------------------------
rm(list=ls())
source('functions_simulation.R')
### Run parameters -------------------------------------------------------------
AOIl <- 1:100
Nsel <- 90 #selection intensity = 10%
r <- 33 ## Ranges 33
p <- 30 ## Plots 30
f <- 0.10
Nrep <- 1
Nchecks <- 3
NF1 <- 30
NF2 <- 30
mgen <- 0
Nsnp <- 20000
AOIlen <- length(AOIl)
Spatial_per <- rep(c(0.1,0.25,0.5,0.75,0.90),len=max(AOIl))
h2_s <- runif(length(AOIl),0.01,0.99)
#Spatial_simu(1,r=r,p=p,f=f,NF1=NF1,NF2=NF2,Nsnp=Nsnp)
### Running code ---------------------------------------------------------------
sysinf <- Sys.info()
os <- sysinf['sysname']
if(os == 'Windows'){
  cl <- makeCluster(6,revtunnel = TRUE)
}else{
  foreach::registerDoSEQ()
  hosts1 <- as.vector((unlist(strsplit(as.character(Sys.getenv("LSB_HOSTS"))," "))))
  host1 <- table(hosts1 )
  H <- names(host1)
  Hn <- as.numeric(host1)## one core free by node
  Hf <- c()
  for(ic in 1:length(H)){
    Hf <- c(Hf,rep(H[ic],l=Hn[ic]))
  }
  hosts1 <- Hf
  cl <-  parallel::makePSOCKcluster(names= hosts1,
                                    outfile = "debug.txt",
                                    master=nsl(Sys.info()['nodename']),
                                    revtunnel = TRUE,
                                    outfile='',
                                    useXDR = TRUE)
}
doParallel::registerDoParallel(cl=cl)
lP <- foreach(i = AOIl,
              .verbose=TRUE,
              #.export = ls(globalenv()),
              .inorder=FALSE,
              .packages = x) %dopar% {
                tryCatch(Spatial_simu(i,r=r,p=p,f=f,NF1=NF1,NF2=NF2,Nsnp=Nsnp),error=function(e){return(NA)})
                }
parallel::stopCluster(cl)
### Results --------------------------------------------------------------------
results <- ldply(lapply(lP[sapply(lP, function(x){sum(is.na(x))==0})],unlist))
colnames(results) <- gsub('inla','',colnames(results))
jobID <- Sys.getenv("LSB_JOBID")
write.csv(results,paste0(jobID,"_results.csv"),row.names = F)


pdf(paste0(jobID,'_Main_Results.pdf'),w=20,h=17)
limitRange <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...,size = 1) + 
    geom_abline(...) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1))
}

## Ge variance
Mag <- max(unlist(c(results[,grep('_ge$',colnames(results))])))
Mig <- min(unlist(c(results[,grep('_ge$',colnames(results))])))
Hg <- ggpairs(results[,grep('_ge$',colnames(results))],mapping=ggplot2::aes(colour = as.factor(results$Spatial_per)),
              upper = list(continuous = wrap("cor", method= "spearman")),
              lower=list(continuous=function(data, mapping, ...) { 
                ggplot(data = data,mapping = mapping,...) + 
                  geom_point(... ,size = 1) + 
                  geom_abline(...) + 
                  scale_y_continuous(limits = c(Mig,Mag)) + 
                  scale_x_continuous(limits = c(Mig,Mag))
              }))+
  labs(title='Genetic variance')
print(Hg)

## Error variance
MErro <- max(unlist(c(results[,c(grep("^sigma2_e",colnames(results)),grep("_erro$",colnames(results)))])))

He <- ggpairs(results[,c(grep("^sigma2_e",colnames(results)),grep("_erro$",colnames(results)))],mapping=ggplot2::aes(colour = as.factor(results$Spatial_per)),
              upper = list(continuous = wrap("cor", method= "spearman")),
              lower=list(continuous=function(data, mapping, ...) { 
                ggplot(data = data, mapping = mapping, ...) + 
                  geom_point(...,size = 1) + geom_abline(...) + 
                  scale_y_continuous(limits = c(0,MErro)) + scale_x_continuous(limits = c(0,MErro))
              }))+
  labs(title='Error')
print(He)
### heritability variance
Mah <- max(unlist(c(results[,grep('^H2_',colnames(results))])))
Mih <- min(unlist(c(results[,grep('^H2_',colnames(results))])))
Hh <- ggpairs(results[,grep('^H2_',colnames(results))],mapping=ggplot2::aes(colour = as.factor(results$Spatial_per)),
              upper = list(continuous = wrap("cor", method= "spearman")),
              lower=list(continuous=function(data, mapping, ...) { 
                ggplot(data = data, mapping = mapping, ...) + 
                  geom_point(...,size = 1) + geom_abline(...) + 
                  scale_y_continuous(limits = c(Mih,Mah)) + scale_x_continuous(limits = c(Mih,Mah))
              }))+
  labs(title='Herdability')
print(Hh)
### Spatial variance
Mas <- max(unlist(c(results[,c(grep('^sigma2_u',colnames(results)),grep('_spatial$',colnames(results)))])))
Mis <- min(unlist(c(results[,c(grep('^sigma2_u',colnames(results)),grep('_spatial$',colnames(results)))])))
Hs <- ggpairs(results[,c(grep('^sigma2_u',colnames(results)),grep('_spatial$',colnames(results)))],mapping=ggplot2::aes(colour = as.factor(results$Spatial_per)),
              upper = list(continuous = wrap("cor", method= "spearman")),
              lower=list(continuous=function(data, mapping, ...) { 
                ggplot(data = data, mapping = mapping, ...) + 
                  geom_point(...,size = 1) + geom_abline(...) + 
                  scale_y_continuous(limits = c(Mis,Mas)) + scale_x_continuous(limits = c(Mis,Mas))
              }))+
  labs(title='spatial variance')
print(Hs)
## likelihood 
Mal <- max(unlist(c(results[,grep('^loglik_',colnames(results))])))
Mil <- min(unlist(c(results[,grep('^loglik_',colnames(results))])))
Hl <- ggpairs(results[,grep('^loglik_',colnames(results))],mapping=ggplot2::aes(colour = as.factor(results$Spatial_per)),
              upper = list(continuous = wrap("cor", method= "spearman")),
              lower=list(continuous=function(data, mapping, ...) { 
                ggplot(data = data, mapping = mapping, ...) + 
                  geom_point(...,size = 1) + geom_abline(...) + 
                  scale_y_continuous(limits = c(Mil,Mal)) + scale_x_continuous(limits = c(Mil,Mal))
              }))+
  labs(title='log likehood')
print(Hl)

### GV correlation
resultsM <- melt(results[,c(grep("LocationID|Spatial_per|^Corr_gv_BLUP",colnames(results)))],id=c('LocationID','Spatial_per'))

BPc <- ggplot(resultsM,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='correlation true value and BLUP')
print(BPc)

### GV CRPS
resultsCRPS <- melt(results[,c(grep("LocationID|Spatial_per|^CRPS_gv_",colnames(results)))],
                    id=c('LocationID','Spatial_per'))
BPcrps <- ggplot(resultsCRPS,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='CRPS GV and predict')
print(BPcrps)

### us correlation
resultsCoru <- melt(results[,c(grep("LocationID|Spatial_per|^Corr_us_",colnames(results)))],
                    id=c('LocationID','Spatial_per'))
BPCoru <- ggplot(resultsCoru,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='Correlation true spatial effects and predict')
print(BPCoru)

## us CRPS
resultsCRPSu <- melt(results[,c(grep("LocationID|Spatial_per|^CRPS_us_",colnames(results)))],
                     id=c('LocationID','Spatial_per'))
BPcrpsu <- ggplot(resultsCRPSu,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='CRPS true spatial value')
print(BPcrpsu)
## Bias
resultsBiasge <- melt(results[,c(grep("LocationID|Spatial_per|^Bias_sigma2_ge_",colnames(results)))],
                      id=c('LocationID','Spatial_per'))
BPcrps <- ggplot(resultsBiasge,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='Bias Sigma2 GE')
print(BPcrps)

resultsBiase <- melt(results[,c(grep("LocationID|Spatial_per|^Bias_sigma2_erro_",colnames(results)))],
                     id=c('LocationID','Spatial_per'))
BPcrps <- ggplot(resultsBiase,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='Bias Sigma2 Erro')
print(BPcrps)

resultsBiasH2 <- melt(results[,c(grep("LocationID|Spatial_per|^Bias_H2_",colnames(results)))],
                      id=c('LocationID','Spatial_per'))
BPcrps <- ggplot(resultsBiasH2,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='Bias H2')
print(BPcrps)

resultsBiasU <- melt(results[,c(grep("LocationID|Spatial_per|^Bias_sigma2_u_",colnames(results)))],
                     id=c('LocationID','Spatial_per'))
BPbiasU <- ggplot(resultsBiasU,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='Bias sigma u')
print(BPbiasU)

## some selection 
resultsSel <- melt(results[,c(grep("LocationID|Spatial_per|^Nsel_BLUP",colnames(results)))],
                   id=c('LocationID','Spatial_per'))
BPsel <- ggplot(resultsSel,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title=paste('Selection top ',Nsel,' 100 GEs'))
print(BPsel)

## some selection 
resultsSel <- melt(results[,c(grep("LocationID|Spatial_per|^NORloglik_",colnames(results)))],
                   id=c('LocationID','Spatial_per'))
BPsel <- ggplot(resultsSel,aes(y=value,x=variable,color=variable))+
  geom_boxplot()+
  facet_wrap(~Spatial_per)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_grid(cols = vars(Spatial_per))+
  labs(title='likelihood models versus no spatial')
print(BPsel)
dev.off()
