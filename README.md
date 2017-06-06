# penalized-blockmodels

This folder contains datasets and scripts associated to the paper

Signorelli, M., Wit, E. C. (2017). A penalized inference approach to stochastic blockmodelling of community-structure in the Italian Parliament. Journal of the Royal Statistical Society: Series C.

The material is divided into two folders: 'data_analysis' and 'simulations'.

The folder 'data_analysis' allows to reproduce the results presented in Section 5 of the paper. To reproduce the results of the analysis, run script '1-model-estimation.R'. The reduced graphs can be reproduced by running '2-create-reduced-graphs.R'.

The folder 'simulations' contains the scripts to reproduce the results of the simulations in the Supplementary Material. The networks can be generated running script '1_network_generation.R', whereas model estimation, model selection and the evaluation of accuracies can be done with script '2-model_selection.R'