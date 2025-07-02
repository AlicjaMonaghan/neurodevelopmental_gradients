# This script visualizes additional reviewer analyses for the paper entitled
# "Canonical neurodevelopmental trajectories of structural and functional
# manifolds" (2025), published in eLife. Written by Alicja Monaghan in June-July
# 2025, at the MRC Cognition and Brain Sciences Unit (University of Cambridge).

library(ggplot2)
library(stringr)

setwd('Desktop/neurodevelopmental_gradients/')

# Load the relationship between age and alignment with the principal structural
# and functional gradients in CALM, without Procrustes rotation.
calm_fc_alignment = read.csv('data/sensitivity/calm_principal_fc_alignment_with_age.csv')
calm_sc_alignment = read.csv('data/sensitivity/calm_principal_sc_alignment_with_age.csv')







# Load the structure-function coupling estimates from our manifold eccentricity 
# measure and the field-standards.
coupling_df = read.csv('data/structure.function/nodal.coupling.comparison.csv')
# Load the group-level (child vs adolescent) structural and functional gradients
# from baseline NKI
nki_child_fc = read.csv('data/sensitivity/nki_child_fc.csv')
nki_adolescent_fc = read.csv('data/sensitivity/nki_adolescent_fc.csv')
developmental_stages = c("child", "adolescent")
modalities = c("sc", "fc")
for (i in 1:length(modalities)){
  for (j in 1:length(developmental_stages)){
    nki_developmental_df = get(paste0("nki_", developmental_stages[j], "_", modalities[i]))
    ggplot(data = nki_developmental_fc, mapping = aes(x=G1, y=G0, colour=label)) +
      geom_point(size = 2) +
      labs(x = "Gradient 2", y = "Gradient 1", title = str_to_title(developmental_stages[j]))
  }
}



ggplot(data=coupling_df, mapping=aes(x=eLife.measure)) +
  geom_smooth(aes(y=baum.2020), method="lm", se=F, color="#F66D7A") +
  geom_smooth(aes(y=vazquez.2019), method="lm", se=F, color="#F3B584") +
  labs(x = "Manifold structure-function coupling", y = "Existing measures") +
  theme(panel.background = element_blank(), axis.line = element_line()) 


  
