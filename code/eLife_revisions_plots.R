# This script visualizes additional reviewer analyses for the paper entitled
# "Canonical neuro-developmental trajectories of structural and functional
# manifolds" (2025), published in eLife. Written by Alicja Monaghan in September
# -October 2025, at the MRC Cognition and Brain Sciences Unit (University of
# Cambridge).

rm(list = ls())
library(ggplot2)
library(dplyr)
library(scales)
library(ggpubr)
library(ggseg)
library(ggsegSchaefer)
library(ggcorrplot)
library(ggtext)
library(cowplot)
library(forcats)
library(ggsignif)
library(R.matlab)
library(tidyverse)
library(paletteer)

setwd('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients/')

# Load the relationship between head motion, age, and alignment of the principal
# communicability and functional connectivity gradients in CALM.
procrustes_alignment_sc = read.csv('data/sensitivity/calm_principal_sc_alignment_with_age.csv')
procrustes_alignment_fc = read.csv('data/sensitivity/calm_principal_fc_alignment_with_age.csv')
# Bind together, and add a factor level to indicate modality
procrustes_alignment = rbind(
  transform(procrustes_alignment_sc, modality=factor("sc")),
  transform(procrustes_alignment_fc, modality=factor("fc"))) %>%
  mutate(predictor = recode(X, scan_age = "Age", meanFWD = "Motion", scan_age = "Age"))
# Load the relationship between the variance explained of individual-level
# functional connectomes explained by the group-level gradient in CALM
calm_fc_var_explained = read.csv('data/calm/individual_variance_explained_by_group_CALM_fc.csv')
# Load the standard deviation of structural and functional gradients across
# development in CALM
calm_std_sc = read.csv('data/calm/sc_std_development.csv') 
calm_std_fc = read.csv('data/calm/fc_std_development.csv') 
# Subset by the columns we want to plot, and bind into a single data frame!
columns_to_retain = c('std_raw', 'scan_age')
calm_std_development = rbind(
  transform(calm_std_sc[columns_to_retain], modality=factor("sc")),
  transform(calm_std_fc[columns_to_retain], modality=factor("fc"))) 
# Load the GLM nodal results to assess the effect of referral type on nodal
# gradient values. We focus on communicability, as this had the more effects 
# that survived FDR correction compared to functional connectivity.
communicability_nodal_value_comparison = read.csv('data/nodal_value_comparison_structural.csv') %>%
  # Add a column to signify whether the region passed FDR correction
  mutate(significance = factor(ifelse(FDR_p < 0.05, "Significant", "NS"))) %>%
  # We're only going to fill in regions which survived FDR correction
  mutate(fill_value = ifelse(FDR_p < 0.05, t, NA))
# Load the nodal coupling data frame across the 4 measures
nodal_coupling = read.csv('data/structure.function/nodal.coupling.comparison.csv') %>%
  # Make column names aesthetically pleasing, and force the labels and dates to 
  # be on separate lines. 
  rename(c("eLife\n(2025)" = "eLife.measure", "PNAS\n(2020)" = "baum",
           "PNAS\n(2019)" = "vazquez", "PNAS\n(2022)" = "park"))
# Load the relationship between the sensorimotor-association axis and each of 
# the 4 coupling metrics
coupling_SA = read.csv('data/structure.function/coupling_SA_alignment.csv') %>%
  # Create a factor level to indicate which measure is novel
  mutate(novelty = factor(ifelse(measure == "eLife\n(2025)", "novel", "prior")))
# Load the child-adolescent gradient correlations!
child_adolescent_corr = read.csv('data/sensitivity/child_adolescent_gradient_comparison.csv') %>%
  mutate(dataset = factor(dataset, labels = c("calm", "nki"))) %>%
  mutate(modality = factor(modality, labels = c("sc", "fc"))) %>%
  mutate(gradient = factor(gradient, labels = c("Primary", "Secondary")))
# Load the prediction accuracies of structure-function coupling
prediction_accuracies = read.csv('data/coupling.phenotypic/empirical_prediction_accuracies_all.csv') %>%
  # Re-order the coupling divisions according to increasing mean accuracy
  reshape2::melt(value.name = "accuracy") %>% mutate(variable = fct_reorder(variable, accuracy, .fun = mean)) %>%
  # Capitalize 'G' in 'global' coupling 
  mutate(variable = fct_recode(variable, "Global" = "global"))
# Load the GAMM predicted fit data for working memory predicting coupling within
# the default-mode network, as well as the cognitive factor scores and raw
# coupling values!
gamm_fit_dmn = read.csv("data/structure.function/GAMM_fit_DMN.csv")
# Separate into two data frames based on high or low Factor 2 scores
low_gamm_fit_dmn = gamm_fit_dmn %>% subset(group_level == "low")
high_gamm_fit_dmn = gamm_fit_dmn %>% subset(group_level == "high")
cognition.structure.function.df = read.csv("data/structure.function/cognition.structure.function.df") %>%
  # Split participants into high or low scorers on Factor 2
  mutate(factor2_grouping = ifelse(Factor2 < median(Factor2), "low", "high"))
# Load the dominance analysis results for SC and FC in CALM
calm_sc_dominance = read.csv('data/calm/calm_sc_dominance.csv')
calm_fc_dominance = read.csv('data/calm/calm_fc_dominance.csv')
# Combine together...
calm_dominance = bind_rows(calm_sc_dominance, calm_fc_dominance, .id = c("modality")) %>%
  mutate(modality = recode_factor(modality, `1` = "sc", `2` = "fc")) %>%
  mutate(Variable_name = recode_factor(
    Variable_name, 'X1' = 'Weighted clustering', 'X2' = 'Within-module degree',
    'X3' = 'Participation coefficient', 'X4' = 'Degree centrality')) %>%
  rename(c("measure" = "Variable_name"))

# Plot the distribution of neuro-diversity effects (t-statistics from a GLM
# framework) for nodal communicability gradient values. 
nodal_communicability_surface_plot = ggplot() + 
  geom_brain(data=communicability_nodal_value_comparison, 
             atlas=schaefer7_200, mapping=aes(fill=fill_value), color="black") +
  scale_fill_gradientn(colors= c("#8a0f63", "#FFFFFF","#f4b100"), na.value = "white") +
  guides(color=FALSE, fill=guide_colorbar(ticks.colour=NA)) + theme_void() +
  theme(legend.position = "bottom", legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(2.3, 'cm'), legend.text = element_blank(),
        legend.title = element_blank())
ggsave('data/calm/neurodiversity_tstat.png', nodal_communicability_surface_plot, dpi=700,
       width=10, height=5)

# Plot the concordance between of primary and secondary gradients derived from 
# childhood or adolescence, across modalities and data sets. 
child_adolescent_corr_plot = 
  ggplot(data=child_adolescent_corr, mapping=aes(
    x=factor(gradient), y=corr, color=dataset, fill=dataset, alpha=modality)) +
  geom_point(position=position_dodge(width=0.5), size=5, shape=16) +
  labs(x = "Gradient", y = "Child-adolescence concordance (*r*)") +
  theme(axis.text = element_text(size=15), axis.title = element_markdown(size=15),
        axis.line = element_line(), panel.background = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0.96, 1), labels = c("0.96", "0.97", "0.98", "0.99", "1")) +
  scale_colour_manual(values = c("#440154FF", "#1F968BFF")) +
  scale_alpha_manual(values=c("sc"=1, "fc"=0.5))

# Plot the regression coefficients (age, sex, and mean FWD) for the Procrustes
# SC and FC analysis in CALM. Add the 95% confidence intervals!
procrustes_coefficients_CALM = 
  ggplot(data=procrustes_alignment, aes(x=estimates, y=reorder(predictor, estimates, decreasing=TRUE), color=modality)) +
  geom_pointrange(mapping=aes(xmin=CI_Lower, xmax=CI_Upper), position=position_dodge(0.5)) +
  labs(x=bquote(beta~"(95% CI) of pre-Procrustes alignment"), y="") +
  geom_vline(xintercept=0, linetype="dashed", color="grey") + 
  theme(panel.background=element_blank(), legend.position = "none", axis.text=element_text(size=15),
        axis.line = element_line(), axis.title=element_text(size=15)) +
  scale_x_continuous(labels=c("-0.50", "-0.25", "0", "0.25")) +
  scale_color_manual(values=c("#8a0f63", "#f4b100")) +
  coord_cartesian(clip = "off")

# Plot the developmental increase in variance explained of individual-level
# functional connectomes by the group-level gradient in CALM.
sequential_colour_palette = c("#323280FF", "#323280CC", "#32328099", "#32328066", "#32328033")
calm_fc_var_explained_plot = 
  ggplot(data=calm_fc_var_explained, mapping=aes(x=scan_age, y=G1, colour=scan_age)) + 
  geom_point(size=3) + geom_smooth(method="lm", colour="black") + 
  labs(x = "Age (years)", y = "% Variance explained in individual FC\n by group-level G1") +
  theme_minimal() + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.position="none") +
  scale_x_continuous(limits=c(6, 19), breaks = seq(6, 18, 3)) +
  scale_y_continuous(limits=c(2, 12), breaks = seq(2, 12, 2)) + 
  scale_fill_gradientn(colors=sequential_colour_palette) + 
  scale_color_gradientn(colors=sequential_colour_palette)

# Plot the relationship between age and the normalized standard deviation of 
# communicability and functional connectivity gradients in the baseline session
# of CALM.
calm_std_development_plot = 
  ggplot(data=calm_std_development, aes(x=scan_age, y=std_raw, colour=scan_age)) +
  labs(x = "Age (years)", y = "Gradient variation") +
  facet_wrap(~modality, scales="free_y", nrow = 2) + geom_point(size=3) + theme_minimal() +
  theme(strip.text = element_blank(), axis.text = element_text(size=15),
        axis.title = element_text(size=15)) +
  geom_smooth(method="lm", color="black") +
  scale_x_continuous(limits = c(6, 19), breaks = seq(6, 18, 3)) +
  scale_color_gradientn(colors = sequential_colour_palette) +
  guides(color=guide_colorbar(ticks.colour=NA))
 
# Assemble the sub-plots for Figure 1. Include a blank plot in which we will
# add the mediating effect of motion for the relationship between age and FC
# alignment in NKI.
figure1_compiled = plot_grid(
  child_adolescent_corr_plot, calm_fc_var_explained_plot,
  calm_std_development_plot,procrustes_coefficients_CALM, ncol=4, axis = "h")
ggsave('data/Figure1_additional_plots_compiled.png', dpi=700, height=5, width=20, plot=figure1_compiled)

# To compare the 4 structure-function coupling metrics among themselves, create
# a correlation matrix using the Spearman's rank correlation coefficient.
coupling_comparison_corr = cor(nodal_coupling[,1:4], method="spearman") %>% round(digits=2)
coupling_corrplot = 
  ggcorrplot(coupling_comparison_corr, lab=TRUE, lab_col = "white", type="lower",
           outline.color = "white", ggtheme = theme_void(), lab_size=10) +
  scale_fill_gradient2(low="white", mid ="#D282B7FF", high = "#791A60FF", midpoint=0.5, limit=c(0,1)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.height=unit(1.5, 'cm'), legend.key.width = unit(2, 'cm'),
        legend.text = element_blank(), axis.text.x = element_text(angle=0, vjust=0.5, hjust=0.5, size=25),
        axis.text.y = element_text(size=25)) + guides(fill=guide_colorbar(ticks.colour=NA)) +
  scale_y_discrete(position="right")
ggsave('data/structure.function/coupling_measure_comparison_corrplot.png', dpi=700, 
       plot=coupling_corrplot, height=5, width=10)

# We will correlate each of the 4 structure-function coupling measures with the 
# sensorimotor-association axis, and plot these coefficients as a lolly plot.
coupling_SA_plot = 
  ggplot(coupling_SA, aes(x = reorder(measure, corr), y=corr)) +
  geom_segment(aes(x=measure, xend=measure, yend=corr, y=0), linewidth=0.7) +
  geom_point(aes(fill=factor(novelty)), size=15, shape=22, color=alpha("white", 0)) +
  geom_hline(yintercept=0, linetype="dashed", size=2, color = "grey") +
  labs(x = "Structure-function coupling measure", y = "SA axis alignment (Spearman's *r*)") +
  theme(legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(), axis.text = element_text(size=25),
        axis.title = element_markdown(size=25)) +
  scale_fill_manual(values = c("#9898BF", "#672975")) +
  scale_y_continuous(
    breaks=seq(0, -0.7, -0.1), labels=c("0", "-0.1", "-0.2", "-0.3", "-0.4", "-0.5", "-0.6", "-0.7"))
ggsave('data/structure.function/coupling_SA_comparison.png', dpi=700,
       plot=coupling_SA_plot, height=8, width=8)

# Plot the distribution of prediction accuracies for working memory from 
# structure-function coupling! A plotting error appears about excluding certain
# values beyond the scale range, but this is due to not visualizing all of the 
# pairwise comparisons. 
prediction_accuracies_plot = 
  ggplot(data=prediction_accuracies, mapping=aes(x=variable, y=accuracy)) +
geom_boxplot(outlier.shape=3, color="#2D708EFF") +
  stat_compare_means(
    hide.ns = TRUE, method = "t.test", tip.length = 0, label = "p.signif",
    symnum.args = list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label.y = c(0.50, NA, 0.55, 0.60),
    comparisons = list(c("SomMot", "DorsAttn"), c("DorsAttn", "Vis"), c("Vis", "Global"), c("Global", "Default"))) +
labs(x = "Network", y="Prediction accuracy of working memory<br>for structure-function coupling (*r*)") +
theme(panel.background=element_blank(), axis.line=element_line(), axis.text.y = element_text(size=13),
      axis.title = element_markdown(size=13), axis.text.x = element_text(size=13, color=c("black", "black", "black", "black", "#2D708EFF"))) +
  scale_y_continuous(labels=c("-0.25", "0", "0.25", "0.50"), breaks=c(seq(-0.25, 0.50, 0.25)))

# Plot the developmental trajectories of structure-function coupling within the 
# default-mode network, as a function of high or low scores for the factor
# representing working memory. 
dmn_coupling_cognition = 
  ggplot(data=cognition.structure.function.df, aes(x=scan_age, y=Default)) +
  geom_point(alpha=0.85, size=2, aes(colour=Factor2)) + geom_line(aes(group=id, colour=Factor2)) +
  scale_colour_gradient2(low="#20A387FF", high="#33638DFF", mid = "white", midpoint = 0) +
  geom_smooth(data=low_gamm_fit_dmn, aes(x=scan_age, y=fit), color="#20A387FF") +
    geom_ribbon(data=low_gamm_fit_dmn, aes(ymin = lwr, ymax = upr), fill="#20A387FF", alpha=0.2) +
    geom_smooth(data=high_gamm_fit_dmn, aes(x=scan_age, y=fit), color="#33638DFF") +
    geom_ribbon(data=high_gamm_fit_dmn, aes(ymin = lwr, ymax = upr), fill="#33638DFF", alpha=0.2) +
  labs(x = "Age (years)", y = "<span style='color:#33638DFF;'>Default</span> network structure-function coupling") +
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text=element_text(size=13), 
        axis.title = element_markdown(size=13)) +
   scale_x_continuous(limits=c(6, 20), breaks=seq(6, 18, 3)) +
  guides(color=guide_colorbar(ticks.colour=NA))
  scale_y_continuous(labels=c("0", "0.10", "0.20", "0.30", "0.40"))

figure5_compiled = plot_grid(prediction_accuracies_plot, dmn_coupling_cognition, ncol=2, axis="h")
ggsave('data/structure.function/Figure5_compiled.png', height=5, width=10, dpi=700, plot=figure5_compiled)

# Plot the relative contribution of 4 graph theory metrics to structural and
# functional connectivity gradients, through a dominance analysis!
dominance_analysis_lolly = 
  ggplot(data=calm_dominance, mapping=aes(x=Dominance_percentage, y=measure)) +
  geom_errorbarh(aes(xmin=0, xmax=Dominance_percentage, group=modality),
                 color="black", height=0, position=position_dodge(width=0.7)) +
  geom_point(aes(fill=modality, color=modality), position=position_dodge(width=0.7), shape=22, size=10) +
  labs(x = "CALM dominance %", y = "Measure of integration/segregation") +
  theme(panel.background=element_blank(), axis.line=element_line(), legend.position = "none",
        axis.text = element_text(size=20), axis.ticks.y = element_blank(),
        axis.title = element_text(size=20), plot.margin = margin(.5,.5,.5,.5, "cm")) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 50)) +
  scale_fill_manual(values = c("#9898BF", "#672975")) +
  scale_color_manual(values = c("#9898BF", "#672975")) 
ggsave('data/calm/dominance_analysis_lolly.png', dpi=700, height=6, width=9, plot=dominance_analysis_lolly)



