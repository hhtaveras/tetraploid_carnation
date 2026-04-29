library(GWASpoly)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(reshape2)
source("s010.R")
source("s011.R")

gfile<-"input/d001.csv"
pfile<-"input/d002.csv"

# GWASpoly
data <- read.GWASpoly(ploidy=4, pheno.file=pfile, geno.file=gfile,
                      format="numeric", n.traits=8, delim=",")

# GWAS loop for all triats
for(t in colnames(data@pheno)[2:ncol(data@pheno)]){
      
      # MAF filter and K matrix 
      params <- set.params(n.PC = 1, MAF = 0.1)    
      setKmat <- set.K(data,LOCO=F,n.core=4)
       
      
      # testing markers for significance
      data.scan <- GWASpoly(data=setKmat,
                                 models=c("additive", 
                                          "1-dom",
                                          "2-dom"),
                                 traits = t, params=params, n.core=4)
      
      
      
      forqq<-data.scan
      forqq@map$Chrom<-rep(1, length(forqq@map$Chrom))
      
      p1<-qqbo.plot(forqq,trait= t) + ggtitle(label=t)
      p11<-p1 + theme_bw(base_size = 11) + theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(), legend.position = "top",
        legend.title = element_blank(), 
        #plot.margin = unit(c(0.4,0,0,0), "cm"),
        legend.text = element_text(size = 9)
      )
      
      meff <- set.threshold(data.scan,method="M.eff",level=0.05, n.core = 4)
      
      p2 <- manbo.plot(meff,traits= t)
      p22<-p2 + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + 
          ylim(0, ceiling(max(meff@scores[[1]], na.rm = T))) + theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(), panel.grid = element_blank(),
            plot.margin = unit(c(1,0.5,0,0), "cm"), legend.position = "top",
            legend.text = element_text(size = 9)
          )
      
      p3<-plot_grid(p22, p11, rel_widths = c(1.8,1), rel_heights = c(1,1))
      
      ggsave( p3,
              file = paste0("figs/manqq", t, ".png"), width=16, height=6, units="in")
      
      # get markers by score
      df<-meff@scores[[1]]
      df$marks<-rownames(df)
      
      mm<-melt(df)
      mmo<-mm[order(mm$value, decreasing = T),]
      mmod<-mmo[!duplicated(mmo$marks),]
      
      saveRDS(mmod, paste0("input/", t, "_topScores.Rds"))
      
      # marker characterization
      
      tryCatch({
        
        qtl <- get.QTL(data=meff, traits= t, models= c("additive", "1-dom", "2-dom") ,bp.window=5e6)
        qtl
        
        fit.ans <- fit.QTL(data=meff,
                           trait=t,
                           qtl=qtl[,c("Marker","Model")],
                           fixed=data.frame(Effect="env",Type="factor"))
        
        if(length(unique(fit.ans$Marker)) == 1){
          
          marks<-as.matrix(data@geno[,colnames(data@geno) %in% unique(fit.ans$Marker)])
          colnames(marks)<-unique(fit.ans$Marker)
          
        } else {
          marks<-as.matrix(data@geno[,colnames(data@geno) %in% unique(fit.ans$Marker)])
        }
        
        # calculate MAF
        mall.freq<-function(x){
          
          mall<-sum(x, na.rm = T)
          nas<-sum(unlist(lapply(x, is.na)))
          totnall<-(length(x)-nas)*4
          maf <- mall/totnall
          return(maf)
        }
        
        (mafs<-apply(marks, 2, mall.freq))    
        
        mreps<-c()
        for(i in 1:length(fit.ans$Marker)){
          for(j in 1:length(mafs)){
            
            if(fit.ans$Marker[i]==names(mafs)[j])
              mreps<-append(mreps, mafs[j])
            
          }
        } 
        
        tab2<-data.frame(Marker=fit.ans$Marker,
                         Chr=fit.ans$Chrom,
                         Pos=fit.ans$Position,
                         MAF=mreps,
                         Score=qtl$Score,
                         Effect=qtl$Effect,
                         PVE=fit.ans$R2*100
                         
                         )
        rownames(tab2)<-NULL
        
        tab2o<-tab2[order(tab2$Chr, tab2$Pos),]
        
        saveRDS(tab2o, paste0("input/", t, "_sigMarks.Rds"))
        
      }, error = function(e) {
        cat("no significant markers \n")
        cat("ERROR :", conditionMessage(e), "\n")
      })
      
  }

# table with significant markers of all traits

traits<-c("earliness", "flower_filling", "flower_size", "branching", "habit", "plant_height", "vigour")

for(i in 1:length(traits)){
  assign(paste0("sigs_", traits[i]) , readRDS(paste0("input/", traits[i], "_sigMarks.Rds")))
}

for(i in 1:length(traits)){
  temp<-get(paste0("sigs_", traits[i]))
  temp$trait<-rep(traits[i], nrow(temp))
  assign(paste0("sigs_", traits[i]), temp)
}

tab1<-data.frame()
for(i in 1:length(traits)){
  temp<-get(paste0("sigs_", traits[i]))
  tab1<-rbind(tab1, temp)
}
tab1


