as.five.marks<-function(z){
  snps<-data.frame(matrix(NA, nrow = length(z), ncol = 0))
  for(i in 1:length(z)){
    if(z[i] == 0){
      snp<-c(1, 0, 0, 0, 0)
    } else if(z[i] == 1){
      snp<-c(0, 1, 0, 0, 0)
    } else if(z[i] == 2){
      snp<-c(0, 0, 1, 0, 0)
    } else if(z[i] == 3){
      snp<-c(0, 0, 0, 1, 0)
    } else if(z[i] == 4){
      snp<-c(0, 0, 0, 0, 1)
    } 
    snps<-rbind(snps, snp)
  } 
  colnames(snps)<-paste0(colnames(z)[1], ".", seq(1, 5))
  return(snps)
}

make.D.mat<-function(z){
  mark.list<-apply(z, 2, as.five.marks)
  mtx<-data.frame(matrix(NA, nrow = length(nrow(z)), ncol = 0))
  for(i in 1:length(mark.list)){
    colnames(mark.list[[i]])<-c(paste0(names(mark.list)[i], ".", seq(1, 5)))
  }
  names(mark.list)<-NULL
  mtx<-cbind(mtx, mark.list)
  rownames(mtx)<-rownames(z)
  return(mtx)
}
