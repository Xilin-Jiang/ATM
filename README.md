## ATM

Age-dependent topic modelling (ATM) is a method for inferring comorbidity profiles for individuals at Biobank Scale. Details of the Method is available in the paper [Age-dependent topic modelling of comorbidities in UK Biobank identifies disease subtypes with differential genetic risk](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2).

ATM assigns to each individual topic weights for several disease topics; each disease topic reflects a set of diseases that tend to co-occur as a function of age, quantified by age-dependent *topic loadings* for each disease. The model assumes that for each disease diagnosis, a topic is sampled based on the individual’s *topic weights* (which sum to 1 across topics, for a given individual), and a disease is sampled based on the individual’s age and the age-dependent *topic loadings* (which sum to 1 across diseases, for a given topic at a given age). The model generalises the latent dirichlet allocation (LDA) model by allowing topic loadings for each topic to vary with age. 
![My Image](ATM_schematic.png)

For bug reports, please email: <xilinjiang@hsph.harvard.edu>. 

## Installation

You can install the development version of ATM from [GitHub](https://github.com/ATM) with:

```r
# install.packages("devtools")
devtools::install_github("Xilin-Jiang/ATM")
```

## Quick start

Run ATM on diagnosis data to infer topic loadings and topic weights from diagnosis data. Note that runing ATM on 100K individuals would take ~30min (default number of inference is 5 runs; the function will pick the best fit). If the data set is too small for inferring its disease topics and the goal is to infer patient-level *topic weights* (i.e. assign comorbidity profiles to individuals based on the set of diseases they have), please use `loading2weights`. The input data should be format data as `HES_age_example`; first column is individual ID, second column is the disease code; third column is the age-at-diagnosis. 

Note for each individual, we only keep the first onset of each diseases. Therefore, if there are recurrent incidences of the same disease code for the same individual, the rest will be ignored.

```r
library(ATM)
# head(HES_age_example)
ATM <- wrapper_ATM(HES_age_example, 10, CVB_num = 1)
```

If the goal is obtaining the  *topic weights* for a group of individuals to learn about their comorbidity profile, there is no need to infer the comorbidity *topic loadings*. Following code below to map the example diagnosis history (example data `HES_age_example`) to the optimal disease topics inferred from UK Biobank HES data. Details are in [Inferring comorbidity profiles for individuals](#inferring-comorbidity-profiles-for-individuals) section.

```r
new_weights <- loading2weights(HES_age_example, ds_list = UKB_349_disease, topics = UKB_HES_10topics)
```

`UKB_HES_10topics` is an internal dataset containing topic loadings inferred from 349 diseases in the UK Biobank HES data. You could substitute it to disease topics inferred from other populations, with the same data format (a tensor of shape $age \times disease number \times topic number$). The output will be the topic weights of each individual in two formats: (1) `new_weights$topic_weights` representing the how much weight each individual have profile (sum to one across topics for each individual); (2) `new_weights$incidence_weight_sum` representing the cumulative weights across diseases (sum across topics equals the number of diagnosis for each individual). 

To visualise the topic loadings, use `plot_age_topics` function. Details are provided in [Visualise the comorbidity topic loadings](#visualise-the-comorbidity-topic-loadings) section. 

```r
disease_list <- UKB_349_disease %>% 
  left_join(disease_info_phecode_icd10, by = c("diag_icd10"="phecode" )) %>% 
  pull(phenotype)
topic_id <- 1 # topic id
plot_age_topics(disease_names = disease_list,
        trajs = UKB_HES_10topics[30:80,,topic_id])
```

## Internal data example

We provide example simulated data in the package. `UKB_349_disease` is the list of 349 diseases (Phecode) that have more than 1000 incidences in the UK Biobank HES data. `HES_age_example` is an example data simulated using the comorbidity distribution in UK Biobank; for inferring disease topics using ATM, you should format the data as `HES_age_example`, which requires individual id, disease diagnosis, and age-at-diagnosis.  `UKB_HES_10topics` is the optimal disease topic from UK Biobank HES data set, using the 349 diseases.

Though ATM could be run on any valid coding system, we recommend using Phecode for ATM to reduce coding redundancy in systems such as ICD-10. To map from ICD-10 code to Phecode, use function `icd2phecode`. `icd2phecode` make use of ICD-10 to phecode mapping which are saved as internal data in ATM package: `phecode_icd10cm` maps between ICD-10-CM to Phecode; `phecode_icd10` maps ICD-10 to Phecode; `disease_info_phecode_icd10` saves the disease names of 1755 Phecodes, use `UKB_349_disease %>% left_join(disease_info_phecode_icd10, by = c("diag_icd10"="phecode" ))`.

## Data preparation

ATM inference is based on age-at-diagnosis information of many diseases. We use the long format to encode the patient id, disease code, and age-at-diagnosis information are provided, instead of using a data matrix where each row is an individual and each column is a disease. This data format save memory as only a small proportion of diseases are diagnosed for each individual. `HES_age_example` is a data example, where each entry contains one diagnosis entry, with individual, disease, and age information. 

The default disease encoding of many biobanks are ICD-10; ATM support any coding system but we recommend using Phecode system which groups ICD-10 codes that represent the same disease. To map data from ICD-10 codes to Phecode, use `icd2phecode` function:  

```r
phecode_data <- icd2phecode(HES_icd10_example)
```

`icd2phecode` maps ICD-10 or ICD-10-CM codes to the Phecodes; when one ICD-10/ICD-10-CM is mapped to multiple Phecodes, it will choose the Phecode that collects the largest number of ICD-10 codes (to reduce the number of Phecodes in the data, which is always good for comorbidity analysis). Before use the function, you should remove all marks such as period in the ICD-10/ICD-10CM coes and only keeps number and capital letters. For example, "I25.1" should be changed to "I251".  

#### All Of Us users (SNOMED)
We note that [All Of Us](https://databrowser.researchallofus.org) and several primary care EHRs release their disease conditions in [SNOMED](https://www.researchallofus.org/faq/what-is-snomed/). For analysis using *All Of Us* alone, you could directly run ATM using the SNOMED coding. However, for those aiming to perform cross-population comparison of comorbidity topics, we provide an internal function  `ATM:::snomed2phecode`, which is used in the same way as `icd2phecode`. This function could map SNOMED code to Phecode, making it comparable to analysis on data sets using ICD-10 coding (mapping ICD-10 codes to Phecodes using `icd2phecode`; e.g. UK Biobank).   

## Inferring disease topics using diagnosis data

If you have an EHR data set with age-at-diagnosis information across many diseases, you could use ATM to infer topic loadings and topic weights. Inferring ATM topic loadings is computational expensive, and the inferred topic loadings usually represents the pattern for the specific data set and should not be extended to other populations, unless they were inferred from large comprehensive biobanks and mapped to populations with similar healthcare system. If the data set is small and the goal is to infer patient-level topic weights (i.e. assign comorbidity profiles to individuals based on their diseases), please use `loading2weights` in the next section. The input data should be formatted as `HES_age_example`; first column is individual ids, second column is the disease code; third column is the age-at-diagnosis. 

One reason that ATM inference is computational expensive is that you need to run multiple models to choose the best number of disease topics in the dataset. ATM does not automatically choose the best number of topics as each model (of different topic numbers) should be run in parallel and you should compare the [ELBO](https://en.wikipedia.org/wiki/Evidence_lower_bound) or *prediction odds ratio* (discussed below) to choose the best fit. In the following example, `topic_num` is the number of topics ( $K$ in [math details](#generative-process-of-atm) section), which for a common EHR data you should choose between 5 to 15; `CVB_num` is the number of random model initialisations, where multiple ATM inferences will be performed and the best model fit will be returned; you are recommended to choose larger number for this parameter if computational power permitting (default is 10);  `ATM_results$multiple_run_ELBO_compare` in the following section provides the ELBOs of all the runs (i.e. for `CVB_num=10` you will get 10 ELBOs), the run with highest ELBO is kept. Use `?wrapper_ATM` to get the details of the function. 

```r
# head(HES_age_example)
ATM_results <- wrapper_ATM(rec_data=HES_age_example, topic_num = 10, CVB_num = 1)
print(ATM_results$multiple_run_ELBO_compare)
```

To choose the optimal model structure that fits the data, running `wrapper_ATM` for each model structure (number of topics and parametric form of curves) and comparing the  `multiple_run_ELBO_compare` for each model structure. The optimal model should has the highest average ELBO across runs. 

An more rigorous way to choose the best model structure is using *prediction odds ratio* defined in the [ATM paper](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2). To perform this analysis, first split the data into training data and testing data, based on individual ids (a common mistake is splitting the diagnosis, where diagnoses of the same individual are presented in both the training and testing data; this would cause using testing data at training stage). Code below provides an example where 20% of the individuals are used sampled as testing data and rest as training. The *topic loadings* inferred from training data is used to compute the *prediction odds ratio* on the testing data.  

```r
testing_ids <- HES_age_example %>% group_by(eid) %>% slice(1) %>% ungroup() %>% sample_frac(size = 0.2) 
training_data <- HES_age_example %>% anti_join(testing_ids, by = "eid")
testing_data <- HES_age_example %>% semi_join(testing_ids, by = "eid")
ATM_results <- wrapper_ATM(rec_data=training_data, topic_num = 10, CVB_num = 1)
testing_prediction_OR <- prediction_OR(testing_data = testing_data, ds_list = ATM_results$ds_list, topic_loadings = ATM_results$topic_loadings)
# print(testing_prediction_OR$OR_top1, testing_prediction_OR$OR_top2, testing_prediction_OR$OR_top5) 
```

Note `ds_list` is a require input of `prediction_OR` as it specify the disease order of the `topic_loadings` input as well as their prevalence in the training data for computing the baseline prediction odds. All diseases in the `testing_data` that are not included in `ds_list` will be discarded. The prediction odds ratio is the odds predicted by ATM versus the odds from a naive prediction using disease prevalence in the training data. Since ATM predict the probability of all diseases simultaneously (multinomial probability), we compute the odds ratio using the probability that the target disease is within the top 1%, 2%, or 5% (i.e. `testing_prediction_OR$OR_top1`,`testing_prediction_OR$OR_top2`, `testing_prediction_OR$OR_top5` ) of the disease codes predicted by ATM respectively. 

## Inferring comorbidity profiles for individuals

In many scenarios, we are not interested in inferring a new set of topics. Instead, for a new set of individuals with medical history, we want to obtain their individual comorbidity profiles (i.e. having CVD related comorbidities). ATM provides *topic weights* which encode comorbidity profile at patient level. To be more specific, using disease topics from [ATM paper](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2), if one individual has elevated CVD topic weight, this individual has excess comorbidities related to cardiovascular diseases. `loading2weights` function provides an easy handle for inferring topic weights, where the input `rec_data` has the same format as in [previous section](#inferring-disease-topics-using-diagnosis-data) and the default comorbidity topics are 10 topics inferred from UK Biobank common diseases `UKB_HES_10topics`.  Code below maps diagnosis history (contained in the example data `HES_age_example`) to the default disease topics inferred from UK Biobank HES data. 

```r
new_weights <- loading2weights(rec_data=HES_age_example, ds_list = UKB_349_disease, topics = UKB_HES_10topics)
patient_topic_weights <- new_weights$topic_weights
cummulative_disease_weights <- new_weights$incidence_weight_sum
```

The outputs of `loading2weights` has two elements: `topic_weights` are patient-level topic weights, referred to as *topic weights* in the  [ATM paper](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2); this should be used when the aim is to understand individual comorbidity profiles. `incidence_weight_sum` is the sum of *diagnosis-specific topic weights* $z_{sn}$ in the ATM paper; this is a more useful metric for prediction, as it represent the cumulative comorbidity burden of each individual. 

## Visualise the comorbidity topic loadings
To visualize the disease topics inferred using ATM, use the `plot_age_topics` functions. To use this function, you need th specify `disease_names`, disease topics `trajs`, title of the plot `plot_title`, and an optional age starts point `start_age` to adjust the x-axis labels in case the supplied topic loadings does not start from age 0. You could also specify how many disease you want to show using `top_ds`.

`disease_names` and `trajs` could be extracted from the output variable from `ATM_wrapper`. Assuming `ATM_results` is output of `ATM_wrapper`, you could  use`ATM_results$ds_list$diag_icd10` as the input for `disease_names`, which provides the names of the diseases in the order of the *topic loadings*. While you could directly use the disease code, we recommend you use meaningful disease description as it will directly shows on the output figure. Similarly, `trajs` could be obtained by slice one matrix from `ATM_results$topic_loading`. For example `ATM_results$topic_loading[,,3]` is the topic loading from the third topic.    

If not clear about the inputs above, please test and check code below, which provides an example of disease topic visualization. 
```r
disease_list <- UKB_349_disease %>%
  dplyr::left_join(disease_info_phecode_icd10, by = c("diag_icd10"="phecode" )) %>%
  dplyr::pull(phenotype)
topic_id <- 1 # plot the first topic
plot_age_topics(disease_names = disease_list,
        trajs = UKB_HES_10topics[30:80,,topic_id],
        plot_title = paste0("topic ", topic_id),
        start_age = 30,
        top_ds = 7)
```

## Analysing comorbidity subtype of diseases
We reported in the [ATM paper](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2) that patients of the same diseases that have different comorbidities has distinct genetic profiles, even after accounting for the genetic backgrounds of these comorbidities (i.e. after correcting for collider effects). We define these as comorbidity subtypes based on *topic weight*.  *topic weight* could be extracted from either from the output of `wrapper_ATM` function (see [Inferring disease topics using diagnosis data](#inferring-disease-topics-using-diagnosis-data) section) or the output of `loading2weights` function (see [Inferring comorbidity profiles for individuals](#inferring-comorbidity-profiles-for-individuals) section). For example, for analysing comorbidity subtypes from the `wrapper_ATM` outputs, users could extract the *topic weights* which is a $\mathcal{R}^{M \times K}$ matrix. The columns are ordered as the *topic loadings*, while the rows are ordered as the `ATM_results$patient_list`.  
```r
ATM_results <- wrapper_ATM(rec_data=HES_age_example, topic_num = 10, CVB_num = 1)
print(ATM_results$multiple_run_ELBO_compare)
subtypes_atm <- data.frame(individual_id = ATM_results$patient_list, topic_weights = ATM_results$topic_weights)
```
Users could then select the set of patients (by subsetting based on `individual_id` of `subtypes_atm`) or which comorbidity subtypes they want to analysis (by selecting columns of the `subtypes_atm`). As an example, in Figure 5 of the [ATM paper](https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2), we analysed the correlation of PRS with topic weights; PRS enriched in the type 2 diabete patients with elevated CVD topic suggested the CVD subtype patients are more driven by genetic risk factors. 

The `topic_weights` of the output of `loading2weights` function has the same format as `subtypes_atm` from example above. 

## Generative process of ATM

![My Image](ATM_generative_process.png)

We constructed a Bayesian hierarchical model to infer latent risk profiles for common diseases.  In summary, the model assumes there exist a few disease topics that underlie many common diseases.  Each topic is age-evolving and contain risk trajectories for all diseases considered. An individual's risk for each diseases is determined by the weights of all topics. The indices in this note are as follows:
$$s= 1,...,M;$$
$$n= 1,...,N_s;$$
$$i= 1,...,K;$$
$$j= 1,...,D;$$
where $M$ is the number of subjects, $N_s$ is the number of records within $s^{th}$ subject, $K$ is number of topics, and $D$ is the total number of diseases we are interested in.
The generative process (Supplementary Figure 1) is as follows: 

- $\theta \in \mathcal{R}^{M \times K}$ is the topic weight for all individuals, each row of which ( $\theta_s \in \mathcal{R}^{K}$ ) is assumed to be sampled from a Dirichlet distribution with parameter $\alpha$. $\alpha$ is set as a hyper parameter.
  $$\theta_s \sim Dir(\alpha).$$

- $\mathbf{z} \in \{1,2,...,K\}^{\sum_s N_s}$ is the topic assignment for each diagnosis $\mathbf{w}  \in \{1,2,...,D\}^{\sum_s N_s}$. Note the total number of diagnoses across all patients are $\sum_s N_s$. The topic assignment for each diagnosis is generated from a multinoulli distribution with parameter equal to $s^{th}$ individual topic weight. 
  $$z_{sn} \sim Multi(\theta_s).$$

- $\beta(t) \in \mathcal{F}(t)^{K \times D}$ is the topic which is $K \times D$ functions of age $t$. $\mathcal{F}(t)$ is the class of functions of $t$. At each plausible $t$, the following is satisfied:
  $$\sum_j \beta_{ij}(t) = 1.$$
  In practice we use softmax function to ensure above is true and add smoothness by constrain $\mathcal{F}(t)$ to be spline or polynomial functions:
  $$\beta_{ij}(t) = \frac{\exp(p_{ij}^T \phi (t))}{\sum_{j} \exp(p_{ij}^T \phi (t))},$$
  where $p_{ij} = ( p_{ij,d} ), \; d = 1,2,...,P$ is the vector of parameter for the topic loading functions; $P$ is the degree of freedom than controls the smoothness; $\phi (t)$ is polynomial and spline basis for age $t$.

- $w \in \{1,2,...,D\}^{\sum_s N_s}$ are observed diagnoses. The $n^{th}$ diagnosis of $s^{th}$ individual $w_{sn}$ is sampled from the topic $\beta_{z_{sn}}(t)$ chosen by $z_{sn}$:
  $$w_{sn} \sim Multi(\beta_{z_{sn}}(t_{sn})),$$
  here $t_{sn}$ is the age of the observed age-at-onset of the observed diagnosis $w_{sn}$.
