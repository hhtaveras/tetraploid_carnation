library(rrBLUP)
library(WeightIt)
library(ggplot2)
library(reshape2)

# read in design matrices and pheno

ZZ<-readRDS("input/d003.Rds")
psZZ<-readRDS("input/d004.Rds")
ZZ5<-readRDS("input/d005.Rds")

pheno<-read.csv("input/d002.csv", sep = ",", dec = ".", header = T)
rownames(pheno)<-pheno$genotype

# nested-CV GWAS results: one file per CV run, holding that run's validation
# set ($val) and the top marker scores computed using ONLY that run's 230
# training individuals ($top). Unlike the static input/<trait>_topScores.Rds
# used in s005.R (ranked once on all 287 genotypes), this means no validation
# phenotype ever influences which markers get used as fixed effects for that
# same fold - no selection leakage.
cv.dir<-"gwas_cv_runs"

# where the per-model correlation results get saved (one .Rds per model,
# each a named list of per-trait NRUNS x nmarks matrices)
out.dir<-"input"
dir.create(out.dir, showWarnings = FALSE)

# CV parameters

NRUNS <- 50

# traits
traits<-colnames(pheno)[2:ncol(pheno)]

# number of markers tested
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

#####################
# Additive tetraploid
#####################

for(trait in traits){
  assign(paste0("corr_A_", trait),
         matrix(NA, nrow = NRUNS, ncol = length(nmarks),
                dimnames = list(NULL, paste0("nmarks_", nmarks))))
}

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
# (loop over runs first so each run's .Rds file is read only once, then
# reused for every trait and every marker count)
for(i in 1:NRUNS){

  setTxtProgressBar(pb, i)
  flush.console()

  run<-readRDS(file.path(cv.dir, sprintf("run_%04d.rds", i)))
  vv<-run$val

  for(trait in traits){

    yy<-data.frame(y = pheno[,trait], row.names = pheno$genotype)
    ZZfm<-ZZ[rownames(yy),]

    trn<-setdiff(rownames(yy), vv)

    yy.val<-data.frame(trait = yy[vv,"y"], row.names = vv)
    yy.trn<-data.frame(trait = yy[trn,"y"], row.names = trn)

    ZZfm.val<-ZZfm[vv,]
    ZZfm.trn<-ZZfm[trn,]

    tops<-run$top[[trait]]
    # this fold's top-ranked markers can include a few SNPs that pass
    # GWASpoly's MAF filter on the 230 training individuals but fall just
    # below the 0.1 MAF threshold used (on all 287 genotypes) to build ZZ -
    # restrict to markers actually present in the design matrix first
    avail.marks<-tops$marks[tops$marks %in% colnames(ZZfm)]
    topmarks<-head(avail.marks, max(nmarks))

    for(nm_i in seq_along(nmarks)){

      nm<-nmarks[nm_i]

      if(nm == 0){

        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)

      } else { # for markers added as fixed effects

        XX<-as.matrix(ZZfm[,head(topmarks, nm)])
        colnames(XX)<-head(topmarks, nm)
        XX<-cbind(intercept = rep(1, nrow(XX)), XX)

        XX.trn<-XX[trn,]
        XX.val<-XX[vv,]

        # check if fixed effects matrix has full rank, if not make into a full rank matrix
        if(qr(XX.trn)$rank < ncol(XX.trn)){
          XX.trn<-make_full_rank(XX.trn)
          XX.val<-XX.val[,colnames(XX.trn), drop = FALSE]
        }
      }

      ms1<-tryCatch(
        mixed.solve(yy.trn$trait, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML",
                    bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE),
        error = function(e){ cat("ERROR :", conditionMessage(e), "\n"); NULL }
      )

      if(!is.null(ms1)){
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u
        rr<-cor(yy.val$trait, y.hat)
        temp<-get(paste0("corr_A_", trait))
        temp[i, nm_i]<-rr
        assign(paste0("corr_A_", trait), temp)
      }
    }
  }
}
#############################################################################

corr_A_list<-setNames(lapply(traits, function(tr) get(paste0("corr_A_", tr))), traits)
saveRDS(corr_A_list, file.path(out.dir, "corr_A.Rds"))

for(i in 1:length(traits)){
  temp<-as.data.frame(get(paste0("corr_A_", traits[i])))
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
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (A, nested-CV markers)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())

##################################################################
# pseudo diploid #################################################
##################################################################


# CV parameters

NRUNS <- 50
traits<-colnames(pheno)[2:ncol(pheno)]
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

for(trait in traits){
  assign(paste0("corr_ps_", trait),
         matrix(NA, nrow = NRUNS, ncol = length(nmarks),
                dimnames = list(NULL, paste0("nmarks_", nmarks))))
}

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for(i in 1:NRUNS){

  setTxtProgressBar(pb, i)
  flush.console()

  run<-readRDS(file.path(cv.dir, sprintf("run_%04d.rds", i)))
  vv<-run$val

  for(trait in traits){

    yy<-data.frame(y = pheno[,trait], row.names = pheno$genotype)
    ZZfm<-psZZ[rownames(yy),]

    trn<-setdiff(rownames(yy), vv)

    yy.val<-data.frame(trait = yy[vv,"y"], row.names = vv)
    yy.trn<-data.frame(trait = yy[trn,"y"], row.names = trn)

    ZZfm.val<-ZZfm[vv,]
    ZZfm.trn<-ZZfm[trn,]

    tops<-run$top[[trait]]
    # this fold's top-ranked markers can include a few SNPs that pass
    # GWASpoly's MAF filter on the 230 training individuals but fall just
    # below the 0.1 MAF threshold used (on all 287 genotypes) to build ZZ -
    # restrict to markers actually present in the design matrix first
    avail.marks<-tops$marks[tops$marks %in% colnames(ZZfm)]
    topmarks<-head(avail.marks, max(nmarks))

    for(nm_i in seq_along(nmarks)){

      nm<-nmarks[nm_i]

      if(nm == 0){

        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)

      } else {

        XX<-as.matrix(ZZfm[,head(topmarks, nm)])
        colnames(XX)<-head(topmarks, nm)
        XX<-cbind(intercept = rep(1, nrow(XX)), XX)

        XX.trn<-XX[trn,]
        XX.val<-XX[vv,]

        if(qr(XX.trn)$rank < ncol(XX.trn)){
          XX.trn<-make_full_rank(XX.trn)
          XX.val<-XX.val[,colnames(XX.trn), drop = FALSE]
        }
      }

      ms1<-tryCatch(
        mixed.solve(yy.trn$trait, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML",
                    bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE),
        error = function(e){ cat("ERROR :", conditionMessage(e), "\n"); NULL }
      )

      if(!is.null(ms1)){
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u
        rr<-cor(yy.val$trait, y.hat)
        temp<-get(paste0("corr_ps_", trait))
        temp[i, nm_i]<-rr
        assign(paste0("corr_ps_", trait), temp)
      }
    }
  }
}
###############################################################################

corr_ps_list<-setNames(lapply(traits, function(tr) get(paste0("corr_ps_", tr))), traits)
saveRDS(corr_ps_list, file.path(out.dir, "corr_ps.Rds"))

for(i in 1:length(traits)){
  temp<-as.data.frame(get(paste0("corr_ps_", traits[i])))
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
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (ps, nested-CV markers)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())


##################################################################
# additive and non-additive ######################################
##################################################################
# make.D.mat() one-hot-encodes each marker into 5 dosage columns (0-4 copies).
# Those 5 columns always sum to 1, so together with an intercept every single
# marker block is exactly collinear (the classic dummy-variable trap) -
# guaranteed rank deficiency on every nm > 0 fit, not just an occasional fold.
# We drop the dosage-0 column of every marker up front (4 full-rank columns
# per marker instead of 5), and reuse the already one-hot-encoded ZZ5
# (= make.D.mat(ZZ), see s004.R) by column name instead of recomputing
# make.D.mat() 1000x per trait.
#
# Even with 4 cols/marker instead of 5, the fixed-effect matrix still grows
# 4x faster than in the A/ps models. With ~230 training individuals per fold,
# nm = 100 or 120 markers would need 401-481 fixed-effect columns - more
# columns than training rows, which is not identifiable. Those two marker
# counts are dropped for this model only.

# CV parameters

NRUNS <- 50
traits<-colnames(pheno)[2:ncol(pheno)]
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

n.val.expected<-round(nrow(pheno) / 5)
n.trn.expected<-nrow(pheno) - n.val.expected
nmarks.AD<-nmarks[(nmarks * 4 + 1) < n.trn.expected]

if(length(nmarks.AD) < length(nmarks)){
  cat("AD model: dropping nmarks =", paste(setdiff(nmarks, nmarks.AD), collapse = ", "),
      "(4 fixed-effect columns/marker + intercept would be >=",
      n.trn.expected, "training individuals per fold) \n")
}

for(trait in traits){
  assign(paste0("corr_AD_", trait),
         matrix(NA, nrow = NRUNS, ncol = length(nmarks.AD),
                dimnames = list(NULL, paste0("nmarks_", nmarks.AD))))
}

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for(i in 1:NRUNS){

  setTxtProgressBar(pb, i)
  flush.console()

  run<-readRDS(file.path(cv.dir, sprintf("run_%04d.rds", i)))
  vv<-run$val

  for(trait in traits){

    yy<-data.frame(y = pheno[,trait], row.names = pheno$genotype)
    ZZfm<-ZZ5[rownames(yy),]

    trn<-setdiff(rownames(yy), vv)

    yy.val<-data.frame(trait = yy[vv,"y"], row.names = vv)
    yy.trn<-data.frame(trait = yy[trn,"y"], row.names = trn)

    ZZfm.val<-ZZfm[vv,]
    ZZfm.trn<-ZZfm[trn,]

    tops<-run$top[[trait]]
    # see note in the A block above: keep only markers that are actually in
    # the design matrices (ZZ / ZZ5) before taking the top N
    avail.marks<-tops$marks[tops$marks %in% colnames(ZZ)]
    topmarks<-head(avail.marks, max(nmarks.AD))

    # pull this run's top markers' 5 dosage-dummy columns out of the
    # already-computed ZZ5, then drop the dosage-0 (reference) column of
    # each marker so the design is full rank by construction
    cols5<-paste0(rep(topmarks, each = 5), ".", rep(1:5, times = length(topmarks)))
    fullXX5<-ZZ5[rownames(yy), cols5]
    fullXX5<-fullXX5[, -seq(1, ncol(fullXX5), by = 5), drop = FALSE]

    for(nm_i in seq_along(nmarks.AD)){

      nm<-nmarks.AD[nm_i]

      if(nm == 0){

        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)

      } else {

        XX<-as.matrix(fullXX5[,1:(nm * 4)])
        XX<-cbind(intercept = rep(1, nrow(XX)), XX)

        XX.trn<-XX[trn,]
        XX.val<-XX[vv,]

        # safety net only: a dosage class can still be absent by chance from
        # a given training fold even with the reference level dropped
        if(qr(XX.trn)$rank < ncol(XX.trn)){
          XX.trn<-make_full_rank(XX.trn)
          XX.val<-XX.val[,colnames(XX.trn), drop = FALSE]
        }
      }

      ms1<-tryCatch(
        mixed.solve(yy.trn$trait, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML",
                    bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE),
        error = function(e){ cat("ERROR :", conditionMessage(e), "\n"); NULL }
      )

      if(!is.null(ms1)){
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u
        rr<-cor(yy.val$trait, y.hat)
        temp<-get(paste0("corr_AD_", trait))
        temp[i, nm_i]<-rr
        assign(paste0("corr_AD_", trait), temp)
      }
    }
  }
}
###############################################################################

corr_AD_list<-setNames(lapply(traits, function(tr) get(paste0("corr_AD_", tr))), traits)
saveRDS(corr_AD_list, file.path(out.dir, "corr_AD.Rds"))

for(i in 1:length(traits)){
  temp<-as.data.frame(get(paste0("corr_AD_", traits[i])))
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
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (AD, nested-CV markers)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())

######################################
# Additive tetraploid (random markers)
######################################
# Negative control: markers are chosen without looking at phenotypes at all,
# so there is no selection leakage here regardless of nesting. The same
# random marker set is reused across all 1000 runs for a given trait/nm (as
# in s005.R), so run-to-run variation reflects only the CV split, not which
# markers were picked.

# CV parameters

NRUNS <- 1000
traits<-colnames(pheno)[2:ncol(pheno)]
nmarks<-c(0, 1, 2, 5, 10, 20, 50, 100, 120)

set.seed(1)
ranmar.list<-list()
for(trait in traits){
  ranmar.list[[trait]]<-list()
  for(nm in nmarks[nmarks > 0]){
    ranmar.list[[trait]][[paste0("nmarks_", nm)]]<-sample(colnames(ZZ), nm)
  }
}

for(trait in traits){
  assign(paste0("corr_Ar_", trait),
         matrix(NA, nrow = NRUNS, ncol = length(nmarks),
                dimnames = list(NULL, paste0("nmarks_", nmarks))))
}

pb <- txtProgressBar(min = 0, max = NRUNS, style = 3)

# start CV loop ###########################################################
for(i in 1:NRUNS){

  setTxtProgressBar(pb, i)
  flush.console()

  run<-readRDS(file.path(cv.dir, sprintf("run_%04d.rds", i)))
  vv<-run$val

  for(trait in traits){

    yy<-data.frame(y = pheno[,trait], row.names = pheno$genotype)
    ZZfm<-ZZ[rownames(yy),]

    trn<-setdiff(rownames(yy), vv)

    yy.val<-data.frame(trait = yy[vv,"y"], row.names = vv)
    yy.trn<-data.frame(trait = yy[trn,"y"], row.names = trn)

    ZZfm.val<-ZZfm[vv,]
    ZZfm.trn<-ZZfm[trn,]

    for(nm_i in seq_along(nmarks)){

      nm<-nmarks[nm_i]

      if(nm == 0){

        XX.trn<-matrix(rep(1, nrow(ZZfm.trn)), nrow = nrow(ZZfm.trn), byrow = T)
        XX.val<-matrix(rep(1, nrow(ZZfm.val)), nrow = nrow(ZZfm.val), byrow = T)

      } else {

        ranmar<-ranmar.list[[trait]][[paste0("nmarks_", nm)]]
        XX<-as.matrix(ZZfm[,ranmar])
        colnames(XX)<-ranmar
        XX<-cbind(intercept = rep(1, nrow(XX)), XX)

        XX.trn<-XX[trn,]
        XX.val<-XX[vv,]

        if(qr(XX.trn)$rank < ncol(XX.trn)){
          XX.trn<-make_full_rank(XX.trn)
          XX.val<-XX.val[,colnames(XX.trn), drop = FALSE]
        }
      }

      ms1<-tryCatch(
        mixed.solve(yy.trn$trait, Z=ZZfm.trn, K=NULL, X=XX.trn, method="REML",
                    bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE),
        error = function(e){ cat("ERROR :", conditionMessage(e), "\n"); NULL }
      )

      if(!is.null(ms1)){
        y.hat<-as.matrix(XX.val) %*% ms1$beta + as.matrix(ZZfm.val) %*% ms1$u
        rr<-cor(yy.val$trait, y.hat)
        temp<-get(paste0("corr_Ar_", trait))
        temp[i, nm_i]<-rr
        assign(paste0("corr_Ar_", trait), temp)
      }
    }
  }
}
############################################################################

corr_Ar_list<-setNames(lapply(traits, function(tr) get(paste0("corr_Ar_", tr))), traits)
saveRDS(corr_Ar_list, file.path(out.dir, "corr_Ar.Rds"))

for(i in 1:length(traits)){
  temp<-as.data.frame(get(paste0("corr_Ar_", traits[i])))
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
    labs(title = "Boxplot of correlations by N. markers fitted as fixed effects and trait (A random mks, nested-CV folds)",
         x = "Trait",
         y = expression(r(y, hat(y)))) +
    theme_minimal())
