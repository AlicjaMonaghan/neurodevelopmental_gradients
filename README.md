Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
# Neurodevelopmental Gradients
Data and code supporting 'Canonical Neurodevelopmental Trajectories of Structural and Functional Manifolds'. We focus on two data sets of children and adolescents. The first is the Centre for Attention, Learning, and Memory (CALM), aged 6 to 17 years old, referred to the service as having problems with one or more of these cognitive domains, by teachers, SEN-coordinators, and clinical practioners. The second is the Nathan Kline Institute (NKI) Rockland Sample Longitudinal Discovery of Brain Development Trajectories, aged 6 to 19 years old, and a community-ascertained sample from the US. 

## Code 
* [Deriving group and individual level gradients][code1] - This describes conducting diffusion-map embedding to derive group-level and individual-level structural and functional gradients for CALM and NKI. Since both CALM and NKI are managed-access datasets, we cannot provide raw connectomes.
* [Collecting DME outputs][code2] - This pulls DME output from the above code, and formats it with co-variates required for statistical modelling. The formatted outputs are already provided for the user.
* [Figure 1][code3] - Code for plotting group-level DME eigenvectors on the cortical surface, exploring variability in the percentage variance accounted for by the first structural and functional components, and examining individual differences in variance explained.
* 

[code1]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/deriving_group_and_individual_gradients_v2.py
[code2]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/dme_and_metadata.R
[code3]: https://github.com/AlicjaMonaghan/neurodevelopmental_gradients/blob/main/code/Figure_1_open.access.R

