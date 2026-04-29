

# read in pheno data 

pheno<-read.csv("input/d002.csv", sep = ",", dec = ".", header = T)


# histograms per trait 
for (i in 2:ncol(pheno)){
  
  hist(pheno[,i], main = colnames(pheno)[i], probability = T, border = "black")
  curve(dnorm(x, mean = mean(pheno[,i]), sd = sd(pheno[,i])),
        col = "red", lwd = 2, add = TRUE)
  
}

# correlations between traits   

cor_test_matrix <- function(df) {
  n <- ncol(df)
  p_mat <- matrix(NA, n, n)
  colnames(p_mat) <- colnames(df)
  rownames(p_mat) <- colnames(df)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        p_mat[i, j] <- 0
      } else {
        test <- cor.test(df[[i]], df[[j]], use = "pairwise.complete.obs")
        p_mat[i, j] <- test$p.value
      }
    }
  }
  p_mat
}

cor_mat<-cor(pheno[,2:9], use = "pairwise.complete.obs")
cor_mat

p_mat   <- cor_test_matrix(pheno[,2:9])
p_mat
