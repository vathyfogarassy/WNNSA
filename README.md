# WNNSA

# Structure of the folder
- **Datasets**: This folder contains the datasets for three different scenarios for control group selection. In each scenario 100 datasets are available.

  - *Scenario I*: ÁTRNI -> The second scenario models such studies in which fewer descriptive variables are available. In this scenario, each individual is characterized by 1 ordinal and 5 binary variables (). The ordinal variable () represents, for example, 5 age groups, while the binary variables () may represent, for example, the gender of the subject or various diagnoses. In this scenario, 700 individuals are simulated in each dataset and the ratio of the candidate subjects to the treated individuals (case group) in the 100 datasets is between 2.0 and 3.1.
  - *Scenario II*: ÁTRNI -> The third scenario is similar to the second one regarding the attributes of individuals and the total number of subjects in each dataset. However, it simulates a more difficult control group selection problem. Although each dataset still containes 700 individuals, the number of treated individuals (case group) in case of the third scenario is higher than in the second one. While in scenario II, the size of the treated group varies between 24.5 and 33.0 percent of the dataset, in the case of the third scenario, it is between 31.7 and 40.7 percent (the ratio of the candidate subjects to the treated individuals is between 1.5 and 2.2).

  *Attributes*:

  - _id: index of the individual
  - x1, ..., xn: attributes (covariates) of the individuals
  - treated: the value of this attribute indicates if an individual is the member of the case group (treated=1) or if it is an individual which can be selected into the control group (treated=0)
  - ps: propensity score value of the individual
  More details of the scenarios are given in the abovementioned article.

- **Results**: This folder contains the result files of different matching algorithms for Scenario I, II and III. Each output contains the following attributes:

  - _id: index of the individual
  - x1, ..., xn: attributes (covariates) of the individuals
  - treated: the value of this attribute indicates if an individual is the member of the case group (treated=1) or if it is an individual which can be selected into the control group (treated=0)
  - ps: propensity score value of the individual
  - ss_pair: index (_id) of the matched pair resulted by the stratified matching method
  - wnnem_pair: index (_id) of the matched pair resulted by the WNNEM method
  - nn_euk_pair_numberoftrial: index (_id) of the matched pair resulted by the nearest neighbour matching based on the Euclidean distances of the individuals
  - nn_mah_pair_numberoftrial: index (_id) of the matched pair resulted by the Mahalanobis matching (nearest neighbour matching based on Mahalanobis distances)
  - psm_02_pair_numberoftrial: index (_id) of the matched pair resulted by the greedy 1:1 PSM method applying caliper size set as 0.2 of the standard deviation of the logit of the propensity scores.
  - psm_dyn_numberoftrial: index (_id) of the matched pair resulted by the greedy 1:1 PSM method applying caliper size determined dynamically and set at the minimal value for which 1:1 matching can be performed.
