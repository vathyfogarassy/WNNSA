# WNNSA

# Structure of the folder
- **Datasets**: This folder contains the datasets for two different scenarios for control group selection. In each scenario 100 datasets are available.

  - *Scenario I*: ÁTRNI -> In the first scenario, each individual is characterized by 10 binary variables. The binary variables may represent, for example, the gender of the subject or various diagnoses. In this scenario, 1000 individuals are simulated in each dataset and the ratio of the candidate subjects to the treated individuals (case group) in the 100 datasets is around 25%.
  - *Scenario II*: ÁTRNI -> The second scenario simulates a more difficult control group selection problem where each individual is characterzied by 8 variables: 2 ordinal, 5 binary and 2 continuous. The weights are a mix of neagtive and positive effects. Each dataset containes 600 individuals. The number of treated individuals (case group) is around 19%.

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
