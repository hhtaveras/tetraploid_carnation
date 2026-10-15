library(GWASpoly)
library(parallel)
library(reshape2)

limit_threads <- function() {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
limit_threads()


gfile<-"input/d001.csv"
pfile<-"input/d002.csv"

data <- read.GWASpoly(ploidy=4, pheno.file=pfile, geno.file=gfile,
                      format="numeric", n.traits=8, delim=",")

data.full <- data                       # the read.GWASpoly object; never modified
traits    <- colnames(data.full@pheno)[2:9]   # check: the 8 trait columns
params    <- set.params(n.PC = 1, MAF = 0.1)
out.dir   <- "gwas_cv_runs"
dir.create(out.dir, showWarnings = FALSE)

run_one <- function(i) {
  limit_threads() 
  f <- file.path(out.dir, sprintf("run_%04d.rds", i))
  if (file.exists(f)) return(invisible(f))            # resume after a crash
  
  set.seed(i)
  ids <- data.full@pheno$genotype
  vv  <- sample(ids, round(length(ids) / 5))          # validation set
  
  d   <- data.full                                    # fresh copy every run
  trn <- !d@pheno$genotype %in% vv
  d@pheno <- d@pheno[trn, , drop = FALSE]
  d@geno  <- d@geno[d@pheno$genotype, , drop = FALSE]
  # if you use fixed effects: d@fixed <- d@fixed[trn, , drop = FALSE]
  
  invisible(capture.output(suppressMessages(
    dK <- set.K(d, LOCO = FALSE, n.core = 1)          # once per run
  )))
  
  top <- lapply(traits, function(t) {
    invisible(capture.output(suppressMessages(
      scan <- GWASpoly(data = dK, models = c("additive", "1-dom", "2-dom"),
                       traits = t, params = params, n.core = 1)
    )))
    df <- scan@scores[[1]]
    df$marks <- rownames(df)
    mm <- melt(df, id.vars = "marks")
    mm <- mm[order(mm$value, decreasing = TRUE), ]
    mm[!duplicated(mm$marks), ]                       # best model score per marker
  })
  names(top) <- traits
  
  saveRDS(list(seed = i, val = vv, top = top), f)     # val needed for the prediction step
  invisible(f)
}

res <- mclapply(1:1000, run_one, mc.cores = 6, mc.preschedule = FALSE)
