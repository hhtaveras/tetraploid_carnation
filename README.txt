Scripts for the tetraploid carnation GWAS and genomic prediction manuscript.

s001.R
  Reads genotype data (input/d001.csv), computes pairwise modified Roger's
  distances (MRD) between cultivars, and produces the PCoA scatter plot and
  distance heatmap used in Figure 1.

s002.R
  Reads phenotype data (input/d002.csv), produces per-trait histograms and
  the pairwise trait correlation matrix with significance tests used in
  Figure 2.

s003.R
  Runs the genome-wide association study (GWASpoly) on all eight traits
  across three gene-action models (additive, 1-dom, 2-dom), applies the
  Meff significance threshold, and produces the Manhattan and QQ plots
  (Figure 3, Supplementary Fig. S1), the significant-marker table
  (Table 3), and per-trait topScores.Rds / sigMarks.Rds files. Sources
  s010.R and s011.R for the plotting functions.

s004.R
  Builds the three SNP marker encodings (pseudo-diploid, additive,
  additive+non-additive) and saves them as input/d003.Rds, input/d004.Rds,
  input/d005.Rds, which s016.R and s018.R depend on. Also runs the
  baseline (no GWAS markers) 1,000-run cross-validation genomic prediction
  used in Figure 4. Sources s012.R.

s010.R
  Helper functions for Manhattan plots (manbo.plot), sourced by s003.R.

s011.R
  Helper function for QQ plots (qqbo.plot), sourced by s003.R.

s012.R
  Helper function for one-hot encoding SNP dosage into five columns
  (as.five.marks), sourced by s004.R.

s015.R
  Generates the nested cross-validation GWAS runs, for each of 1,000
  80/20 train-validation splits, re-runs GWAS on the training individuals
  only and saves the result as gwas_cv_runs/run_NNNN.rds. This is what
  makes the downstream genomic prediction leakage-free.

s016.R
  Reads gwas_cv_runs/ (from s015.R) and the marker encodings (from
  s004.R), and runs the nested GWAS-assisted genomic prediction across
  all three encodings plus a random-marker control. Saves results as
  input/corr_A.Rds, input/corr_ps.Rds, input/corr_AD.Rds, input/corr_Ar.Rds.

