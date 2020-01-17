# A penalized inference approach to stochastic block modelling of community-structure in the Italian Parliament

This folder contains R scripts related to the publication of Signorelli, M., Wit, E. C. (2018). A penalized inference approach to stochastic block modelling of community-structure in the Italian Parliament. *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 67 (2), 355-369. 
You can read the paper (with open access) here: https://journals.sagepub.com/doi/full/10.1177/1471082X19871128

# About this repository
This repository contains the data and code to reproduce the simulations and data analyses presented in Signorelli and Wit (2018).

The material is divided into two folders: 'data_analysis' and 'simulations'.

The folder 'data_analysis' allows to reproduce the results presented in Section 5 of the paper. To reproduce the results of the analysis, run script '1-model-estimation.R'. The reduced graphs can be reproduced by running '2-create-reduced-graphs.R'.

The folder 'simulations' contains the scripts to reproduce the results of the simulations in the Supplementary Material. The networks can be generated running script '1_network_generation.R', whereas model estimation, model selection and the evaluation of accuracies can be done with script '2-model_selection.R'
