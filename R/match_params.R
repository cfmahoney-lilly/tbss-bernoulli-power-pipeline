# Propensity score formula and matching arguments
# Additional arguments in the match_args list will be passed to the glm().
#  This script contains only matching parameters and is not intended as a target
#  on its own
#  Replace the propensity score formula with covariates relevant to your study
#  design; the sample covariates here correspond to categories of the Elixhauser
#  comorbidity score and are for representative purposes only

# propensity score formula
ps_form <- treatment_group ~ gender + age + chf + carit + valv + pcd + pvd + hypunc +
  hypc + para + ond + cpd + diabunc + diabc + hypothy + rf + ld +
  pud + aids + lymph + metacanc + solidtum + rheumd + coag + obes +
  wloss + fed + blane + dane + alcohol + drug + psycho + depre 

# specify matching parameters, suggest these:
match_method <- "full"
match_distance <- "glm"
match_link <- "logit"
match_caliper <- 0.1

match_args <- list(
  formula = ps_form,
  method = match_method,
  distance = match_distance,
  caliper = match_caliper
)
