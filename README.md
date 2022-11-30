## ATM (document in progress)
Age-dependent topic modelling (ATM) is inferring comorbidity profiles for individuals at Biobank Scale. Details of the Method is available at https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2 
ATM assigns to each individual topic weights for several disease topics; each disease topic reflects a set of diseases that tend to co-occur as a function of age, quantified by age-dependent *topic loadings* for each disease. The model assumes that for each disease diagnosis, a topic is sampled based on the individual’s *topic weights* (which sum to 1 across topics, for a given individual), and a disease is sampled based on the individual’s age and the age-dependent *topic loadings* (which sum to 1 across diseases, for a given topic at a given age). The model generalises the latent dirichlet allocation (LDA) model by allowing topic loadings for each topic to vary with age. 
![My Image](ATM_schematic.png)

# Generative process of ATM
![My Image](ATM_generative_process.png)

## Estimate individual comorbidity weight using the comorbidity profiles inferred from UK Biobank.(XXXinsert code example here; use simulated data to show)

## Inferring comorbidity profiles from individual diagnostic data set. 
(XXX this is computational expensive, therefore not only provide R package, also complie an executable)
