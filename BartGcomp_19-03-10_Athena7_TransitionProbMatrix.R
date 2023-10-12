# transition probabilities

de = newdat[!is.na(outcome_2_lag1) & outcome_2_lag1==1,] #previously engaged
dde = newdat[!is.na(outcome_2_lag1) & outcome_2_lag1==0,] #previously disengaged

det = de[, .N ,by = outcome]
det$N/sum(det$N)

ddet = dde[, .N ,by = outcome]
ddet$N/sum(ddet$N)
