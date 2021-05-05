# PSID_education_eq_k6
#This repository hosts resources for an ongoing project asking the question: "To what extent does employment quality (EQ) mediate the effect of educational attainment on mental health?".

#In order to answer this, the parametric mediational g-formula approach outlined by Lin et al. in their 2017 paper "Parametric Mediational g-Formula Approach to Mediation Analysis with #Time-varying Exposures, Mediators, and Confounders" is used. This approach can decompose the total effect of an exposure on on outcome into randomised interventional analogues of #natural direct and indirect effects. 

#For our project, this approach seems appropriate given 1) EQ is a dynamic factor over the working life (i.e. it is time-varying), 2) education likely affects mental health through #pathways separate from EQ that themselves confound the EQ-mental health relationship (i.e. there is exposure-induced mediator-outcome confounding), and 3) EQ cannot easily be #categorised in a meaningful way without over-simplification. This third factor makes the parametric approach by Lin et al. favourable over the IPW approach by VanderWeele and Tchetgen #Tchetgen (2017).

#For our project, however, existing software operationalising the mediational g-formula cannot support situations with a time-fixed exposure, as educational attainment is for most #individuals during their working life.

#This repo therefore contains:

#1) An R script assigning each individual with sufficient information a PCA-based EQ score.   
#2) An R script creating a final analysis dataset in wide- and long-format for the maximum study population.
#3) An R script containing the mediational g-formula for marginal comparisons.
#4) An R script containing a simplified version of the mediational g-formula for general use.
#5) A .jpeg file for the assumed causal DAG used in this analysis.  