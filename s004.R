library(rrBLUP)
library(ggplot2)
library(reshape2)
source("s012.R")

# read in geno and pheno data

geno<-read.csv("input/d001.csv", sep = ",", dec = ".", header = T)
pheno<-read.csv("input/d002.csv", sep = ",", dec = ".", header = T)

############################################
# format and filter geno data make ZZ matrix

geno0<-geno[,4:ncol(geno)]
rownames(geno0)<-geno$rn
maf <- rowMeans(geno0, na.rm=TRUE) / 4
maf <- ifelse(maf > 0.5, 1-maf, maf)
maf_threshold <- 0.1
keep_snps <- maf >= maf_threshold

geno1 <- geno0[keep_snps,]
ZZ<-t(geno1)

saveRDS(ZZ, "input/d003.Rds")

###############################################################
# make psZZ matrix, all heterozygotes are the same = 1

psZZ<-ZZ

psZZ[psZZ == 2] <- 1
psZZ[psZZ == 3] <- 1
psZZ[psZZ == 4] <- 2

saveRDS(psZZ, "input/d004.Rds")

###############################################################
# make ZZ5 matrix, 5 columns per marker, each allele one column

ZZ5<-make.D.mat(ZZ)

saveRDS(ZZ5, "input/d005.Rds")

###############
# CV parameters 

# number of cross-validation runs 
NRUNS <- 1000

# traits 
traits<-colnames(pheno)[2:ncol(pheno)]

# encodigs 
model<-c("pd", "A", "AD")

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

###############################################################################
# loop with three models ######################################################
###############################################################################

for (m in model){
  
  print(paste0("model = ", m))
  if(m == "pd"){
    ZZm <- psZZ
  } else if (m == "A"){
    ZZm <- ZZ
  } else if (m == "AD"){
    ZZm <- ZZ5
  }
  assign(paste0("co", m) , data.frame(matrix(NA, nrow = NRUNS, ncol = 0)))
  
  for (trait in traits) {
    
    print(paste0("trait = ", trait))
    # make yy
    yy<-data.frame(y = pheno[,trait])
    rownames(yy)<-pheno$genotype
    
    # remove genotypes in genotypic data that are not present in phenotypic data
    ZZfm<-ZZm[!rownames(ZZm) %in% setdiff(rownames(ZZm), rownames(yy)),]
    
    # get same order of genotypes on both 
    yy$name<-rownames(yy)
    yy<-yy[match(rownames(ZZfm), rownames(yy)),]
    yy<-data.frame(y=yy$y, row.names = rownames(yy))
    
    crrs <- c()
    # for no markers added as fixed effects
    
    for (i in 1:NRUNS) {
      

      setTxtProgressBar(pb, i) 
      flush.console()  
      set.seed(i)
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
    
    if(m == "pd"){
      copd <- cbind(copd, crrs)
    } else if (m == "A"){
      coA <- cbind(coA, crrs)
    } else if (m == "AD"){
      coAD <- cbind(coAD, crrs)
    }
    
  }
  if(m == "pd"){
    colnames(copd)<- traits
  } else if (m == "A"){
    colnames(coA)<- traits
  } else if (m == "AD"){
    colnames(coAD)<- traits
  }
  
} 

# format correlation results for plot

copd$model<-rep("pd", nrow(copd))
coA$model<-rep("A", nrow(coA))
coAD$model<-rep("AD", nrow(coAD))

mcopd<-melt(copd)
mcoA<-melt(coA)
mcoAD<-melt(coAD)

corrs<-rbind(mcopd, mcoA, mcoAD)
corrs$variable<-factor(corrs$variable, levels = c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "vigour", "plant_height"))

# plot results 

(p<-ggplot(corrs, aes(x = variable, y = value, fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(title = "Boxplot of correlations for each encoding",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())
