library(rrBLUP)
library(WeightIt)
library(ggplot2)
library(reshape2)
source("s012.R")

# read in design matrices and pheno 

ZZ<-readRDS("input/d003.Rds")
psZZ<-readRDS("input/d004.Rds")
ZZ5<-readRDS("input/d005.Rds")

pheno<-read.csv("input/d002.csv", sep = ",", dec = ".", header = T)

#####################
# Additive tetraploid
#####################

# CV parameters 

# number of cross-validation runs 
NRUNS <- 1000

# traits 
traits<-colnames(pheno)[2:ncol(pheno)]

# number of markers tested 
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for (trait in traits) {
  
  assign(paste0("corr_A_", trait), data.frame(matrix(NA, nrow = NRUNS, ncol = 0)))
  print(paste0("trait = ", trait))
  
  # make yy
  yy<-data.frame(y = pheno[,trait])
  rownames(yy)<-pheno$genotype
  
  # remove genotypes in genotypic data that are not present in phenotypic data
  ZZfm<-ZZ[!rownames(ZZ) %in% setdiff(rownames(ZZ), rownames(yy)),]
  
  # get same order of genotypes on both 
  yy$name<-rownames(yy)
  yy<-yy[match(rownames(ZZfm), rownames(yy)),]
  yy<-data.frame(y=yy$y, row.names = rownames(yy))
  
  tops<-readRDS(paste0("input/", trait, "_topScores.Rds"))
  
  # for no markers added as fixed effects
  
  for(nm in nmarks){
    print(paste0("n.markers = ", nm))
    
    if(nm == 0){  
  
      crrs <- c()
      for (i in 1:NRUNS) {
        
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)
        
        yy.trn$X1 <- yy.trn$trait
        
        ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                         bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u  
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))       
        crrs <- append(crrs, rr)
        
      }
      temp<-get(paste0("corr_A_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_A_", trait), temp)
      
    } else { # for markers added as fixed effects 
      
      crrs <- c()
      XX<-as.matrix(ZZfm[,head(tops$marks, nm)])
      XX<-cbind(rep(1, nrow(XX)), XX)
      colnames(XX)<-c("intercept", head(tops$marks, nm))
      
      for (i in 1:NRUNS) {
        
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.val <- as.matrix(XX[vv,])
        XX.trn <- as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
        
        colnames(XX.trn)<-colnames(XX)
        colnames(XX.val)<-colnames(XX)
        
        yy.trn$X1 <- yy.trn$trait
        
        # check if fixed effects matrix has full rank, if not make into a full rank matrix
        if (qr(XX.trn)$rank < ncol(XX.trn)) {
          XX.trn <- make_full_rank(XX.trn)
          
          dim(XX.trn)
          nc <- c(colnames(XX.trn))
          dim(XX.val)
          XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
        } 
        
        tryCatch({
          
          ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                           bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
          
          y.hat<- as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u   
        
          # compare observed and predicted VL for validation set
          (rr <- cor(yy.val, y.hat))       
          crrs <- append(crrs, rr)
          
        }, error = function(e) {
          cat("ERROR :", conditionMessage(e), "\n")
        })
      }
      
      temp<-get(paste0("corr_A_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_A_", trait), temp)
    }
  }
  temp<-get(paste0("corr_A_", trait))
  colnames(temp)<-paste0("nmarks_", nmarks)
  assign(paste0("corr_A_", trait), temp)
}       
#############################################################################

for(i in 1:length(traits)){
  temp<-get(paste0("corr_A_", traits[i]))
  temp$tr<-rep(traits[i], nrow(temp))
  assign(paste0("corr_A_", traits[i]), temp)
}

for(i in 1:length(traits)){
  memp<-get(paste0("corr_A_", traits[i]))
  memp<-melt(memp)
  assign(paste0("mco_A_", traits[i]), memp)
}

corrsA<-data.frame()
for(i in 1:length(traits)){
  corrsA<-rbind(corrsA, get(paste0("mco_A_", traits[i])))
}

corrsA$tr<-factor(corrsA$tr, levels = c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "vigour", "plant_height"))

(p1<-ggplot(corrsA, aes(x = tr, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (A)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())

##################################################################
# pseudo diploid #################################################
##################################################################

# CV parameters 

# number of cross-validation runs 
NRUNS <- 1000

# traits 
traits<-colnames(pheno)[2:ncol(pheno)]

# number of markers tested 
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for (trait in traits) {
  
  assign(paste0("corr_ps_", trait), data.frame(matrix(NA, nrow = NRUNS, ncol = 0)))
  print(paste0("trait = ", trait))
  # make yy
  yy<-data.frame(y = pheno[,trait])
  rownames(yy)<-pheno$genotype
  
  # remove genotypes in genotypic data that are not present in phenotypic data
  ZZfm<-psZZ[!rownames(psZZ) %in% setdiff(rownames(psZZ), rownames(yy)),]
  
  # get same order of genotypes on both 
  yy$name<-rownames(yy)
  yy<-yy[match(rownames(ZZfm), rownames(yy)),]
  yy<-data.frame(y=yy$y, row.names = rownames(yy))
  
  tops<-readRDS(paste0("input/", trait, "_topScores.Rds"))
  
  # for no markers added as fixed effects
  
  for(nm in nmarks){
    print(paste0("n.markers = ", nm))
    
    if(nm == 0){  
      
      crrs <- c()
      
      for (i in 1:NRUNS) {
        
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)
        
        yy.trn$X1 <- yy.trn$trait
        
        ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                         bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u  
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))       
        crrs <- append(crrs, rr)
      }
      
      temp<-get(paste0("corr_ps_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_ps_", trait), temp)
      
    } else {
      
      crrs <- c()
      XX<-as.matrix(ZZfm[,head(tops$marks, nm)])
      XX<-cbind(rep(1, nrow(XX)), XX)
      colnames(XX)<-c("intercept", head(tops$marks, nm))
      
      for (i in 1:NRUNS) {
    
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        
        # set training and validation sets
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.val <- as.matrix(XX[vv,])
        XX.trn <- as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
        
        colnames(XX.trn)<-colnames(XX)
        colnames(XX.val)<-colnames(XX)
        
        yy.trn$X1 <- yy.trn$trait
        
        # check if fixed effects matrix has full rank, if not make into a full rank matrix
        if (qr(XX.trn)$rank < ncol(XX.trn)) {
          XX.trn <- make_full_rank(XX.trn)
          
          dim(XX.trn)
          nc <- c(colnames(XX.trn))
          dim(XX.val)
          XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
        } 
     
        tryCatch({
          
          ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                           bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
          
          y.hat<- as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u   
          
          # compare observed and predicted VL for validation set
          (rr <- cor(yy.val, y.hat))       
          crrs <- append(crrs, rr)
          
        }, error = function(e) {
          cat("ERROR :", conditionMessage(e), "\n")
        })
      }
      temp<-get(paste0("corr_ps_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_ps_", trait), temp)
    }
    
  }
  temp<-get(paste0("corr_ps_", trait))
  colnames(temp)<-paste0("nmarks_", nmarks)
  assign(paste0("corr_ps_", trait), temp)
}
###############################################################################

for(i in 1:length(traits)){
  temp<-get(paste0("corr_ps_", traits[i]))
  temp$tr<-rep(traits[i], nrow(temp))
  assign(paste0("corr_ps_", traits[i]), temp)
}

for(i in 1:length(traits)){
  memp<-get(paste0("corr_ps_", traits[i]))
  memp<-melt(memp)
  assign(paste0("mco_ps_", traits[i]), memp)
}

corrsps<-data.frame()
for(i in 1:length(traits)){
  corrsps<-rbind(corrsps, get(paste0("mco_ps_", traits[i])))
}

corrsps$tr<-factor(corrsps$tr, levels = c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "vigour", "plant_height"))

(p2<-ggplot(corrsps, aes(x = tr, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (ps)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())


##################################################################
# additive and non-additive ######################################
##################################################################

# CV parameters 

# number of cross-validation runs 
NRUNS <- 1000

# traits 
traits<-colnames(pheno)[2:ncol(pheno)]

# number of markers tested 
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for (trait in traits) {
  
  assign(paste0("corr_AD_", trait), data.frame(matrix(NA, nrow = NRUNS, ncol = 0)))
  print(paste0("trait = ", trait))
  
  # make yy
  yy<-data.frame(y = pheno[,trait])
  rownames(yy)<-pheno$genotype
  
  # remove genotypes in genotypic data that are not present in phenotypic data
  ZZfm<-ZZ5[!rownames(ZZ5) %in% setdiff(rownames(ZZ5), rownames(yy)),]
  
  # get same order of genotypes on both 
  yy$name<-rownames(yy)
  yy<-yy[match(rownames(ZZfm), rownames(yy)),]
  yy<-data.frame(y=yy$y, row.names = rownames(yy))
  
  tops<-readRDS(paste0("input/", trait, "_topScores.Rds"))
  
  fullXX<-ZZ[,head(tops$marks, 120)]
  fullXX5<-make.D.mat(fullXX)
  
  # for no markers added as fixed effects
  for(nm in nmarks){
    print(paste0("n.markers = ", nm))
    
    if(nm == 0){  
      
      crrs <- c()
      
      for (i in 1:NRUNS) {
        
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)
      
        yy.trn$X1 <- yy.trn$trait
        
        ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                         bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
        y.hat<- as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u  
      
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))       
        crrs <- append(crrs, rr)
      }
      
      temp<-get(paste0("corr_AD_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_AD_", trait), temp)
      
    } else {
      
      crrs <- c()
      XX<-fullXX5[,1:(nm*5)]
      XX<-cbind(intercept=rep(1, nrow(XX)), XX)
    
      for (i in 1:NRUNS) {
        
        #print(paste0("run = ", i))
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.val <- as.matrix(XX[vv,])
        XX.trn <- as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
        
        colnames(XX.trn)<-colnames(XX)
        colnames(XX.val)<-colnames(XX)
        
        yy.trn$X1 <- yy.trn$trait
        
        # check if fixed effects matrix has full rank, if not make into a full rank matrix
        if (qr(XX.trn)$rank < ncol(XX.trn)) {
          XX.trn <- make_full_rank(XX.trn)
          
          dim(XX.trn)
          nc <- c(colnames(XX.trn))
          dim(XX.val)
          XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
        } 

        ms1<-try(mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML"),
                 silent = TRUE
        )
        
        if (inherits(ms1, "try-error")) {
          cat("ERROR:", attr(ms1, "condition")$message, "\n")
          rr <- NA
          crrs <- append(crrs, rr)
        } else {
          y.hat <- XX.val %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u
          rr <- cor(yy.val, y.hat)
          crrs <- append(crrs, rr)
        }
      }
      temp<-get(paste0("corr_AD_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_AD_", trait), temp)
    }
  }
  temp<-get(paste0("corr_AD_", trait))
  colnames(temp)<-paste0("nmarks_", nmarks)
  assign(paste0("corr_AD_", trait), temp)
}        
###############################################################################

for(i in 1:length(traits)){
  temp<-get(paste0("corr_AD_", traits[i]))
  temp$tr<-rep(traits[i], nrow(temp))
  assign(paste0("corr_AD_", traits[i]), temp)
}

for(i in 1:length(traits)){
  memp<-get(paste0("corr_AD_", traits[i]))
  memp<-melt(memp)
  assign(paste0("mco_AD_", traits[i]), memp)
}

corrsAD<-data.frame()
for(i in 1:length(traits)){
  corrsAD<-rbind(corrsAD, get(paste0("mco_AD_", traits[i])))
}

corrsAD$tr<-factor(corrsAD$tr, levels = c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "vigour", "plant_height"))

(p3<-ggplot(corrsAD, aes(x = tr, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (AD)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())

######################################
# Additive tetraploid (random markers)
######################################

# CV parameters 

# number of cross-validation runs 
NRUNS <- 1000

# traits 
traits<-colnames(pheno)[2:ncol(pheno)]

# number of markers tested 
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for (trait in traits) {
  
  assign(paste0("corr_Ar_", trait), data.frame(matrix(NA, nrow = NRUNS, ncol = 0)))
  
  print(paste0("trait = ", trait))
  # make yy
  yy<-data.frame(y = pheno[,trait])
  rownames(yy)<-pheno$genotype
  
  # remove genotypes in genotypic data that are not present in phenotypic data
  ZZfm<-ZZ[!rownames(ZZ) %in% setdiff(rownames(ZZ), rownames(yy)),]
  
  # get same order of genotypes on both 
  yy$name<-rownames(yy)
  yy<-yy[match(rownames(ZZfm), rownames(yy)),]
  yy<-data.frame(y=yy$y, row.names = rownames(yy))
  
  #tops<-readRDS(paste0("input/", trait, "_maf_01_topScores_", map, ".Rds"))
  
  # for no markers added as fixed effects
  
  for(nm in nmarks){
    print(paste0("n.markers = ", nm))
    
    if(nm == 0){  
      
      crrs <- c()
      for (i in 1:NRUNS) {
        
        #print(paste0("run = ", i))
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        # set training and validation sets,
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)
        
        yy.trn$X1 <- yy.trn$trait
        
        ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                         bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u  
    
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))       
        crrs <- append(crrs, rr)
      }
      
      temp<-get(paste0("corr_Ar_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_Ar_", trait), temp)
      
    } else {
      
      crrs <- c()
      ranmar<-sample(colnames(ZZ), nm)
      XX<-as.matrix(ZZfm[,ranmar])
      colnames(XX)<-ranmar
      XX<-cbind(rep(1, nrow(XX)), XX)
      colnames(XX)[1]<-"intercept"

      for (i in 1:NRUNS) {
        
        setTxtProgressBar(pb, i) 
        flush.console()  
        set.seed(i)
        # set training and validation sets
        vv <- sample(rownames(yy), round(nrow(yy) / 5))
        yy.val <- data.frame(trait = yy[vv,])
        rownames(yy.val) <- vv
        
        yy.trn<-data.frame(trait = yy[setdiff(rownames(yy), rownames(yy.val)),])
        rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
        
        ZZfm.val <- ZZfm[vv,]
        ZZfm.trn <- ZZfm[setdiff(rownames(ZZfm), rownames(ZZfm.val)),]
        
        XX.val <- as.matrix(XX[vv,])
        XX.trn <- as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
        
        colnames(XX.trn)<-colnames(XX)
        colnames(XX.val)<-colnames(XX)
        
        yy.trn$X1 <- yy.trn$trait
        
        # check if fixed effects matrix has full rank, if not make into a full rank matrix
        if (qr(XX.trn)$rank < ncol(XX.trn)) {
          XX.trn <- make_full_rank(XX.trn)
          
          dim(XX.trn)
          nc <- c(colnames(XX.trn))
          dim(XX.val)
          XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
        } 
        tryCatch({
          ms1<-mixed.solve(yy.trn$X1, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML", 
                           bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
          
          y.hat<- as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u   
          
          # compare observed and predicted VL for validation set
          (rr <- cor(yy.val, y.hat))       
          crrs <- append(crrs, rr)
          
        }, error = function(e) {
          cat("ERROR :", conditionMessage(e), "\n")
        })
      }
      temp<-get(paste0("corr_Ar_", trait))
      temp<-cbind(temp, crrs)
      assign(paste0("corr_Ar_", trait), temp)
    }
  }
  
  temp<-get(paste0("corr_Ar_", trait))
  colnames(temp)<-paste0("nmarks_", nmarks)
  assign(paste0("corr_Ar_", trait), temp)
  
}

############################################################################

for(i in 1:length(traits)){
  temp<-get(paste0("corr_Ar_", traits[i]))
  temp$tr<-rep(traits[i], nrow(temp))
  assign(paste0("corr_Ar_", traits[i]), temp)
}

for(i in 1:length(traits)){
  memp<-get(paste0("corr_Ar_", traits[i]))
  memp<-melt(memp)
  assign(paste0("mco_Ar_", traits[i]), memp)
}

corrsAr<-data.frame()
for(i in 1:length(traits)){
  corrsAr<-rbind(corrsAr, get(paste0("mco_Ar_", traits[i])))
}

corrsAr$tr<-factor(corrsAr$tr, levels = c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "vigour", "plant_height"))

(p4<-ggplot(corrsAr, aes(x = tr, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (A random mks)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())
