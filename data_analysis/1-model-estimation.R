#############################################
# This script allows to reproduce the data
# analysis in Section 5 of the paper.
# Use the variable case (= 1, 2, 3, 4) to
# select the legislature of interest, as 
# detailed below
#############################################

library(Matrix); library(glmnet)

# REMARK:
# SET case = 1 for the analysis of the 14th legislature;
# SET case = 2 for the analysis of the 15th legislature;
# SET case = 3 for the analysis of the 16th legislature;
# SET case = 4 for the analysis of the 17th legislature.
case = 4

if (case == 1) load('data\\dataset-XIV-legislature.RData')
if (case == 2) load('data\\dataset-XV-legislature.RData')
if (case == 3) load('data\\dataset-XVI-legislature.RData')
if (case == 4) load('data\\dataset-XVII-legislature.RData')

# Xcov contains the covariates in Table 1
# Xregion contains the dummies for the regional effects
# in Table 2 of the Supplementary material
# Xparty contains the transformed dummies T_r and T_rs
# described in Equation (6), Section 3.1

# standardization of the quantitative covariates (age difference and triads)
Xcov[,4] = Xcov[,4]/sqrt(var(Xcov[,4]))
Xcov[,9] = Xcov[,9]/sqrt(var(Xcov[,9]))

# create the design matrix
X = cbind(Xcov, Xregion, Xparty)

# compute weights for the adaptive Lasso
mle = glm(y ~ as.matrix(X), family = poisson)
mlecoef = mle$coefficients
adlassowei = 1/(mlecoef^2)

# define the penalty factors for glmnet
penalize = adlassowei[-1] # no intercept
penalize[(dim(Xcov)[2]+dim(Xregion)[2]+1):(dim(Xcov)[2]+dim(Xregion)[2]+nparties-1)] = 0 # no penalty on main effects (alpha)

# model fit
nlambda = 100
fit <- glmnet(x = X, y = y, intercept=TRUE, family = "poisson", standardize = FALSE,
               penalty.factor=penalize, nlambda = nlambda)

# BIC:
N = v*(v-1)/2
BIC = (fit$df+1)*log(N) + deviance(fit)

estimates = coef(fit, s = fit$lambda[BIC==min(BIC)])

# parameter estimates for the variables in Table 1:
covariates = estimates[1:(dim(Xcov)[2]+1)]
cbind(c('int', colnames(Xcov)), round(covariates,3))

# parameter estimates for the regional effects in Table 2
# of the Supplementary Material
regional_eff = estimates[(dim(Xcov)[2]+2):(dim(Xcov)[2]+dim(Xregion)[2]+1)]
round(regional_eff, 3)

if (case == 1) save.image('data\\model-XIV-legislature.RData')
if (case == 2) save.image('data\\model-XV-legislature.RData')
if (case == 3) save.image('data\\model-XVI-legislature.RData')
if (case == 4) save.image('data\\model-XVII-legislature.RData')
