## ATM (document in progress)
Age-dependent topic modelling (ATM) is inferring comorbidity profiles for individuals at Biobank Scale. Details of the Method is available at https://www.medrxiv.org/content/10.1101/2022.10.23.22281420v2 
ATM assigns to each individual topic weights for several disease topics; each disease topic reflects a set of diseases that tend to co-occur as a function of age, quantified by age-dependent *topic loadings* for each disease. The model assumes that for each disease diagnosis, a topic is sampled based on the individual’s *topic weights* (which sum to 1 across topics, for a given individual), and a disease is sampled based on the individual’s age and the age-dependent *topic loadings* (which sum to 1 across diseases, for a given topic at a given age). The model generalises the latent dirichlet allocation (LDA) model by allowing topic loadings for each topic to vary with age. 
![My Image](ATM_schematic.png)

# Generative process of ATM
![My Image](ATM_generative_process.png)

We constructed a Bayesian hierarchical model to infer latent risk profiles for common diseases.  In summary, the model assumes there exist a few disease topics that underlie many common diseases.  Each topic is age-evolving and contain risk trajectories for all diseases considered. An individual's risk for each diseases is determined by the weights of all topics. The indices in this note are as follows:
$$s= 1,...,M;$$
$$n= 1,...,N_s;$$
$$i= 1,...,K;$$
$$j= 1,...,D;$$
where $M$ is the number of subjects, $N_s$ is the number of records within $s^{th}$ subject, $K$ is number of topics, and $D$ is the total number of diseases we are interested in.
The generative process (Supplementary Figure \ref{fig:generative}) is as follows: 

- $\theta \in \mathcal{R}^{M \times K}$ is the topic weight for all individuals, each row of which ($\in \mathcal{R}^{K}$) is assumed to be sampled from a Dirichlet distribution with parameter $\alpha$. $\alpha$ is set as a hyper parameter.
$$\theta_s \sim Dir(\alpha).$$

- $\mathbf{z} \in \{1,2,...,K\}^{\sum_s N_s}$ is the topic assignment for each diagnosis $\mathbf{w}  \in \{1,2,...,D\}^{\sum_s N_s}$. Note the total number of diagnoses across all patients are $\sum_s N_s$. The topic assignment for each diagnosis is generated from a multinoulli distribution with parameter equal to $s^{th}$ individual topic weight. 
$$z_{sn} \sim Multi(\theta_s).$$

- $\beta(t) \in \mathcal{F}(t)^{K \times D}$ is the topic which is $K \times D$ functions of age $t$. $\mathcal{F}(t)$ is the class of functions of $t$. At each plausible $t$, the following is satisfied:
$$\sum_j \beta_{ij}(t) = 1.$$
In practice we use softmax function to ensure above is true and add smoothness by constrain $\mathcal{F}(t)$ to be spline or polynomial functions:
$$\beta_{ij}(t) = \frac{\exp(\boldsymbol{p}_{ij}^T \phi (t))}{\sum_{j = 1}^D \exp(\boldsymbol{p}_{ij}^T \phi (t))},$$
where $\boldsymbol{p}_{ij} = \{ p_{ijd} \}, \; d = 1,2,...,P$; $P$ is the degree of freedom than controls the smoothness; $\phi (t)$ is polynomial and spline basis for age $t$.

- $w \in \{1,2,...,D\}^{\sum_s N_s}$ are observed diagnoses. The $n^{th}$ diagnosis of $s^{th}$ individual $w_{sn}$ is sampled from the topic $\beta_{z_{sn}}(t)$ chosen by $z_{sn}$:
$$w_{sn} \sim Multi(\boldsymbol{\beta}_{z_{sn}}(t_{sn})),$$
here $t_{sn}$ is the age of the observed age-at-onset of the observed diagnosis $ w_{sn}$.

## Estimate individual comorbidity weight using the comorbidity profiles inferred from UK Biobank.(XXXinsert code example here; use simulated data to show)

## Inferring comorbidity profiles from individual diagnostic data set. 
(XXX this is computational expensive, therefore not only provide R package, also complie an executable)
