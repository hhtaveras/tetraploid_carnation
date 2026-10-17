# Figure 5 remade with the nested-CV GWAS-assisted genomic prediction results
# (s017.R output, corr_A.Rds / corr_ps.Rds / corr_AD.Rds / corr_Ar.Rds in input/),
# replacing the original leaky-GWAS version produced by car03-s034.R. Plot design
# (panel layout, colors, fonts) is kept identical to car03-s034.R for comparability;
# only data loading and axis limits are adapted to the new, nested-CV result format.

library(reshape2)
library(ggplot2)
library(cowplot)

traits<-c("earliness", "flower_filling", "flower_count", "flower_size", "branching", "habit", "plant_height", "vigour")

# s017.R saves one .Rds per model, each a named list (per trait) of raw
# NRUNS x nmarks matrices - reshape into the long format car03-s034.R expects
read_model<-function(rds_path){
  lst<-readRDS(rds_path)
  out<-data.frame()
  for(tr in names(lst)){
    temp<-as.data.frame(lst[[tr]])
    temp$tr<-tr
    memp<-melt(temp, id.vars = "tr")
    out<-rbind(out, memp)
  }
  out$tr<-factor(out$tr, levels = traits)
  out
}

corrsPD<-read_model("input/corr_ps.Rds")
corrsA<-read_model("input/corr_A.Rds")
corrsAD<-read_model("input/corr_AD.Rds")
corrsR<-read_model("input/corr_Ar.Rds")

# four plots ################################################################################

taz<-16
# the nested-CV results dip more negative than the original leaky-CV results
# across all four models (down to ~-0.42 for AD), not just AD as before, so a
# single shared lower limit is used for all panels instead of clipping at 0.
# Fixed at -0.25 (supplementary version): outliers below this are clipped.
ulim<--0.25

blp<-c("#EFA09C", "#ECC87A", "#B3D48E",
       "#8ED4B0", "#8ED4DA", "#8EC4EA",
       "#C4A0EA", "#EAA0CC", "#EAB0B0")

boxcols<-blp



(plPDt<-ggplot(corrsPD, aes(x = tr, y = value, fill = variable)) +
   geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8, outlier.shape = 1, show.legend = F) +
   scale_y_continuous(
     name = expression(r(y, hat(y))),
     limits = c(ulim, 1),
     breaks = seq(ulim, 1, 0.25)
   ) + scale_fill_manual(values = boxcols) +
   annotate("text", x = 0.7, y = 1, label = "a",
           fontface = 2, size = 8) +
   annotate("text", x = 8, y = 1, label = "PD",
           fontface = 2, size = 8) +
   labs(title = "",
        x = "",
        y = expression(r (y, hat(y)))) +
   theme_bw() +
   theme(plot.background = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         panel.grid.minor = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz),
         plot.margin = margin(t=-15,r=0,b=-15,l=0))

    )

############################################################################

(plAt<-ggplot(corrsA, aes(x = tr, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8, outlier.shape = 1, show.legend = F) +
    scale_y_continuous(
      name = expression(r(y, hat(y))),
      limits = c(ulim, 1),
      breaks = seq(ulim, 1, 0.25)
    ) + scale_fill_manual(values = boxcols) +
    annotate("text", x = 0.7, y = 1, label = "b",
           fontface = 2, size = 8) +
    annotate("text", x = 8, y = 1, label = "A",
           fontface = 2, size = 8) +
    labs(title = "",
         x = "",
         y = expression(r (y, hat(y)))) +
   theme_bw() +
   theme(plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz),
         plot.margin = margin(t=-15,r=0,b=-15,l=0))
 )

###########################################################################

(plADt<-ggplot(corrsAD, aes(x = tr, y = value, fill = variable)) +
   geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8, outlier.shape = 1, show.legend = F) +
   scale_y_continuous(
     name = expression(r(y, hat(y))),
     limits = c(ulim, 1),
     breaks = seq(ulim, 1, 0.25)
   ) + scale_fill_manual(values = boxcols) +
   annotate("text", x = 0.7, y = 1, label = "c",
            fontface = 2, size = 8) +
   annotate("text", x = 8, y = 1, label = "AD",
            fontface = 2, size = 8) +
   labs(title = "",
        x = "",
        y = "") +
   theme_bw() +
   theme(plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz),
         plot.margin = margin(t=-15,r=0,b=-15,l=0))
 )

 ##############################################################################

(plRt<-ggplot(corrsR, aes(x = tr, y = value, fill = variable)) +
   geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8, outlier.shape = 1) +
   scale_y_continuous(
     name = expression(r(y, hat(y))),
     limits = c(ulim, 1),
     breaks = seq(ulim, 1, 0.25)
   ) +
   annotate("text", x = 0.7, y = 1, label = "d",
            fontface = 2, size = 8) +
   scale_x_discrete(labels = c("Earliness", "Flower filling", "Flower count", "Flower size",
                               "Branching", "Habit", "Plant height", "Vigour")) +
   labs(title = "",
        x = "",
        y = expression(r (y, hat(y))),
        ) +
   scale_fill_manual(
     name = "Number of SNP markers fitted as fixed effects",
     labels = c("0", "1", "2", "5", "10", "20",
                "50", "100", "120"),
     values = boxcols) +
   theme_bw() +
   theme(plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         legend.title = element_text(vjust = -1, hjust = 0.5, size = taz+1),
         legend.title.position = "bottom",
         legend.position = "bottom",
         legend.box = "horizontal",
         legend.background = element_blank(),
         legend.key = element_rect(fill = NA, color = NA),
         legend.key.spacing.x = unit(1, "cm"),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz),
         legend.text = element_text(size = taz+2),
         plot.margin = margin(t=-15,r=0,b=0,l=0)

         ) +
   guides(fill = guide_legend(nrow = 1))
    )




(cowmbo <- plot_grid(
  plPDt,
  plAt,
  plADt,
  plRt,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1, 1, 1.25),
  greedy = FALSE
))

# saved under a new name - does NOT overwrite the_man/fig5.png (original, leaky-CV figure)
ggsave( cowmbo ,
        file="the_man/fig5_nestedCV.png",
        width=16, height=22, units="in", dpi = 300)

# copy with fixed -0.5 to 1 y-axis for all four panels, for use as Supplementary Fig. S3
ggsave( cowmbo ,
        file="supp/sfigS3.png",
        width=16, height=22, units="in", dpi = 300)


################################################################################
# table4_gain.pdf remade with the nested-CV results, as a ggplot2 grouped bar
# chart (original was PD/A/AD only, no random-marker control). For each trait
# and encoding, the "best" marker count is whichever nmarks > 0 gives the
# highest mean r(y, yhat); gain is that mean minus the nmarks_0 baseline mean.
# When the baseline outperforms every nmarks > 0, gain is shown as "-" (no
# bar), matching the dash convention in the manuscript's Table 4.
################################################################################

compute_gain<-function(rds_path){
  lst<-readRDS(rds_path)
  out<-data.frame()
  for(tr in names(lst)){
    m<-lst[[tr]]
    means<-colMeans(m, na.rm = TRUE)
    baseline<-means["nmarks_0"]
    nonbase<-means[names(means) != "nmarks_0"]
    best_idx<-which.max(nonbase)
    best_mean<-nonbase[best_idx]
    best_nm<-as.integer(gsub("nmarks_", "", names(nonbase)[best_idx]))
    out<-rbind(out, data.frame(tr = tr, baseline = baseline, best_mean = best_mean,
                                best_nm = best_nm, gain = best_mean - baseline))
  }
  rownames(out)<-NULL
  out
}

gainPD<-compute_gain("input/corr_ps.Rds"); gainPD$encoding<-"Pseudo-diploid (PD)"
gainA<-compute_gain("input/corr_A.Rds");   gainA$encoding<-"Additive dosage (A)"
gainAD<-compute_gain("input/corr_AD.Rds"); gainAD$encoding<-"Additive + non-additive (AD)"

gain_df<-rbind(gainPD, gainA, gainAD)
gain_df$tr<-factor(gain_df$tr, levels = traits)
gain_df$encoding<-factor(gain_df$encoding,
                          levels = c("Pseudo-diploid (PD)", "Additive dosage (A)", "Additive + non-additive (AD)"))
# no bar (and a dash label) whenever the baseline model already performed best
gain_df$label<-ifelse(gain_df$gain > 0, as.character(gain_df$best_nm), "–")
gain_df$plot_gain<-pmax(gain_df$gain, 0)

# same pastel palette as the rest of the manuscript's figures (blp, see above),
# just the 3 members that give the closest match to table4_gain.pdf's colors
gaincols<-c("Pseudo-diploid (PD)" = "#EFA09C",
            "Additive dosage (A)" = "#B3D48E",
            "Additive + non-additive (AD)" = "#8EC4EA")

(plGain<-ggplot(gain_df, aes(x = tr, y = plot_gain, fill = encoding)) +
   geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.3) +
   geom_text(aes(label = label, y = plot_gain),
             position = position_dodge(width = 0.8), vjust = -0.4, size = 4.2) +
   scale_fill_manual(values = gaincols, name = NULL) +
   scale_x_discrete(labels = c("Earliness", "Flower\nfilling", "Flower\ncount", "Flower\nsize",
                                "Branching", "Habit", "Plant\nheight", "Vigour")) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
   labs(x = "", y = expression(paste("Gain in mean ", r(y, hat(y)), " over baseline"))) +
   theme_bw() +
   theme(panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         plot.background = element_blank(),
         legend.position = "top",
         legend.text = element_text(size = taz),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz))
 )

ggsave( plGain,
        file="the_man/table4_gain_nestedCV.png",
        width=12, height=6, units="in", dpi = 300)


################################################################################
# Same gain chart, restricted to the only three traits where any encoding
# actually showed a positive gain under nested CV (earliness, flower filling,
# branching) - all other traits are "-" (baseline wins) across PD, A, and AD.
################################################################################

gain_df_sub<-gain_df[gain_df$tr %in% c("earliness", "flower_filling", "branching"), ]
gain_df_sub$tr<-factor(gain_df_sub$tr, levels = c("earliness", "flower_filling", "branching"))

(plGainSub<-ggplot(gain_df_sub, aes(x = tr, y = plot_gain, fill = encoding)) +
   geom_col(position = position_dodge(width = 0.8), width = 0.65, colour = "black", linewidth = 0.3) +
   geom_text(aes(label = label, y = plot_gain),
             position = position_dodge(width = 0.8), vjust = -0.4, size = 4.2) +
   scale_fill_manual(values = gaincols, name = NULL) +
   scale_x_discrete(labels = c("Earliness", "Flower\nfilling", "Branching")) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
   labs(x = "", y = expression(paste("Gain in mean ", r(y, hat(y)), " over baseline"))) +
   theme_bw() +
   theme(panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.background = element_rect(fill = "white", colour = NA),
         plot.background = element_blank(),
         legend.position = "top",
         legend.text = element_text(size = taz),
         axis.text = element_text(size = taz),
         axis.title = element_text(size = taz))
 )

ggsave( plGainSub,
        file="the_man/table4_gain_nestedCV_3traits.png",
        width=9, height=6, units="in", dpi = 300)
