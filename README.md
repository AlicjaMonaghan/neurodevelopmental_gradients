Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
# Neurodevelopmental Gradients
Data and code supporting 'Canonical Neurodevelopmental Trajectories of Structural and Functional Manifolds'. We focus on two data sets of children and adolescents. The first is the Centre for Attention, Learning, and Memory (CALM), aged 6 to 17 years old, referred to the service as having problems with one or more of these cognitive domains, by teachers, special educational needs coordinators, and clinical practioners. The second is the Nathan Kline Institute (NKI) Rockland Sample Longitudinal Discovery of Brain Development Trajectories, aged 6 to 19 years old, and a community-ascertained sample from the US. Qualified researchers can access these harmonized, longitudinal datasets spanning the neurotypical-neurodivergent continuum and with common dimensions of psychopathology and cognition upon reasonable request. Any questions should be forwarded to alicja.monaghan@outlook.com and Duncan.Astle@mrc-cbu.cam.ac.uk.

![alt text](https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/Figure1.png)

# `connectome construction`
Internal scripts used to clean rsfMRI data and construct functional connectomes, as well as calculating communicability for structural connectomes. Included here for reference. Also includes a MATLAB script assessing the effects of COMBAT harmonisation on CALM data. 

# `code` 
* [Summary statistics][code0] - Summary statistics for cross-sectional and longitudinal portions of CALM and NKI.
* [Cognitive and psychopathology data for CALM and NKI][code11] - Pulls age-standardised scores spanning cognitive and psychopathology in preparation for PCA.
* [Communicability function][code1a] - Weighted communicability function, used to derive structural gradients, based on a formulation by Estrada and Hatano (2008) and adapted by Crofts and Higham (2009).
* [Deriving group and individual level gradients][code1] - This describes conducting diffusion-map embedding (DME) to derive group-level and individual-level structural and functional gradients for CALM and NKI. Since both CALM and NKI are managed-access datasets, we cannot provide raw connectomes.
* [Collecting DME outputs][code2] - This pulls DME output from the above code, and formats it with co-variates required for statistical modelling. The formatted outputs are already provided for the user.
* [Figure 1][code3] - Code for plotting group-level DME eigenvectors on the cortical surface, exploring variability in the percentage variance accounted for by the first structural and functional components, and examining individual differences in variance explained.
* [Spatial permutation test of brain maps][code3a] - Function to generate a p-value for the spatial correlation between parcellated cortical surface maps, developed by Frantisek Vasa[FVref]. Used in Figure 1 to compare gradients across cohorts, and used in Figure 4 to relate structure-function coupling with the sensorimotor-association axis.
* [Gradient developmental stability][code3b] - Conducts all the sensitivity analyses in the final segment of Figure 1 (panels f-j), including effect of group vs individual-level analysis, mediating effect of motion on the relationship between age and alignment, alongside effects of Procrustes rotation.
* [Figure 2][code4] - Plots the relationship between mean manifold eccentricity and global graph theory metrics.
* [Dominance function][code4a] - Dominance analysis MATLAB function written by [Shiyong Yu][shiyongdominance], used for analysis in Figure 2. 
* [Dominance analysis of manifold eccentricity][code4b] - Relating manifold eccentricity with graph theory measures in a dominance analysis (Figure 2c)
* [Generalised additive mixed models (GAMM)][code5] - Contains all GAMM functions, including handling factor-smooth interactions, obtaining GAMM predictions for plotting, and plotting partial age effects on manifolds.
* [Structural and functional GAMMs][code8] - GAMMs to assess developmental effects on structural and functional manifold eccentricity, adjusting for covariates.
* [Computing structure-function coupling and deriving dimensions of psychopathology and cognition][code9] - Updates the CALM and NKI meta-data sheets with correct ages, then computes structure-function coupling and behavioural dimensions. Also assesses impact of missingness. 
* [Developmental effects of structure-function coupling][code8b] - Analyses and visualisation for panels a, b, c, and e of Figure 4.
* [Individual differences in and developmental trajectories of structure-function coupling][code10] - Conducts GAMMs to test developmental effects on structure-function coupling. Visualises individual differences in structure-function coupling as a function of data set and network.
* [Methodological structure-function coupling sensitivity][code8d] - This compares the newly-derived structure-function coupling gradient measure with established measures, alongside testing for neurotypicality effects in coupling magnitude, and calculating coupling for the neurotypical portion of CALM. The developmental effects on neurotypical coupling in CALM are tested in the code below.
* [Neurotypical structure-function coupling][code8c] - This tests for developmental effects in structure-function coupling within the neurotypical portion of CALM, acting as a sensitivity analysis for Figure 4d.
* [Linking psychopathology/cognition with structure-function coupling][code12] - Conducts GAMMs to test developmental relationship between 3 dimensions of psychopathology and 2 of cognition, with structure-function coupling.
* [GAMM cross-validation][code12a] - Cross-validation for the GAMMs predicting structure-function coupling using the second cognitive dimension (working memory). Train-test splits are shuffled by participant ID, with all observations belonging to each participant remaining in each fold. Visualised in Figure 5c.
* [eLife reviewer plots][code12b] - Plots for eLife reviewers, such as visualising the dominance analysis (Figure 2b). Also includes all visualisation in the final line of Figure 1, including the proportion of variance in individual-level FC explained by the group template (Figure 1g), gradient variation across modalities (Figure 1h), and generalised linear model coefficients relating motion and age with individual-level alignment (Figure 1i). Plots the distribution of prediction accuracies relating structure-function coupling and working memory (Figure 5c), as well as their interaction (Figure 5d).

[code0]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/nkir_calm_descriptives_v2.R
[code1a]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/communicability.m
[code1]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/deriving_group_and_individual_gradients_v2.py
[code2]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/dme_and_metadata.R
[code3]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/Figure_1_open.access.R
[code3a]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/perm.sphere.p.R
[FVref]: https://github.com/frantisekvasa/rotate_parcellation
[code3b]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/gradient_developmental_sensitivity.py
[code4]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/Figure_2_open.access.R
[code4a]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/dominance.m
[shiyongdominance]: https://uk.mathworks.com/matlabcentral/fileexchange/168806-dominance-analysis
[code4b]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/eccentricity_graph_theory.m
[code5]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/GAMM.functions.v3.R
[code7]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/structural_and_functional_gamms.R
[code8]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/structure_function_coupling_and_psychopathology.py
[code8b]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/developmental_effects_on_structure_function_relationships.R
[code8c]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/neurotypical.coupling.development.R
[code8d]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/structure_function_coupling_sensitivity.py
[code9]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/structure_function_coupling_individual_differences.R
[code11]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/NKI_and_CALM_cognitive_behavioural_V3.py
[code12]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/psychopathology_cognition_dimensions_structure_function.R
[code12a]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/gamm_coupling_crossvalidation.R
[code12b]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/eLife_revisions_plots.R



