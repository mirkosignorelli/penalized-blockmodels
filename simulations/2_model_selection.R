###############################################
# Choose the simulation (A, B, C or D) and the 
# number of vertices (50, 100, 150, ..., 500)
# in order to reproduce the results in Figure
# 1 of the Supplementary Material
###############################################

source('0_functions.R')
library(Matrix); library(glmnet)

simulation = 'D' # choose the simulation
v = 500 # choose among v = 50, 100, 150, 200, ..., 500

if (simulation == 'A') load('simA.RData')
if (simulation == 'B') load('simB.RData')
if (simulation == 'C') load('simC.RData')
if (simulation == 'D') load('simD.RData')

adj = as.matrix(netobj$adj[1:v,1:v])
adj = adj+t(adj) # fill in the lower triangle (saved with 0 to save memory)
group = netobj$group[1:v]

# check that there are 10 groups and each of them is made by 2 units
if (length(table(group) == 1) < p | sum(table(group) == 1) > 0) {
  adj = as.matrix(netobj$adj[(v+1):(2*v),(v+1):(2*v)])
  adj = adj+t(adj)
  group = netobj$group[(v+1):(2*v)]
}

###########################################
###### CREATION OF THE DESIGN MATRIX ######
###########################################

nodes_ord = order(group)
adj_ord = adj[nodes_ord, nodes_ord]
rm(adj, netobj)
group_ord = group[order(group)]
adj_ord[lower.tri(adj_ord, diag = T)] = 0
adj_ord = as(adj_ord, 'sparseMatrix')

data = matrix(nrow = 0, ncol = 2 + (p-1) + p*(p-1)/2)
# first and second column: y and weights
# remaining columns: X

for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    inodes = (1:v)[group_ord == i]
    jnodes = (1:v)[group_ord == j]
    subadj = as.matrix(adj_ord[inodes,jnodes])
    freq = as.data.frame(table(subadj))
    freq = cbind(as.numeric(levels(freq[,1])),as.numeric(freq[,2]))
    l = dim(freq)[1]
    if (l == 0) warning(paste('something went wrong with groups', i, j))
    if (l >0) {
      add = cbind(freq, matrix(rep(dtrasfvec(i,j,p),l), nrow = l, byrow = T))
      data = rbind(data, add)
    }
  }
}
rm(i,j)

for (i in 1:(p)) {
  inodes = (1:v)[group_ord == i]
  subadj = as.matrix(adj_ord[inodes,inodes])
  subadj = subadj[upper.tri(subadj, diag = F)]
  freq = as.data.frame(table(subadj))
  freq = cbind(as.numeric(levels(freq[,1])),as.numeric(freq[,2]))
  l = dim(freq)[1]
  if (l == 0) warning(paste('not enough observations for group', i))
  if (l >0) {
    add = cbind(freq, matrix(rep(dtrasfvec(i,i,p),l), nrow = l, byrow = T))
    data = rbind(data, add)
  }
}
rm(i)

y = data[,1]
weights = data[,2]
X = data[,-c(1,2)]


######################################
########## MODEL ESTIMATION ##########
########### ADAPTIVE LASSO ###########
######################################
mle = glm(y ~ X, weights = weights, family = poisson)
relev = abs(mle$coefficients[(p+1):(p+p*(p-1)/2)])
penalty_weights = 1/(relev^2)

penalty_factor = numeric()
penalty_factor[1:(p-1)] = 0
penalty_factor[(p):(p-1+p*(p-1)/2)] = penalty_weights

nlambda = 100
lmr = 0.00001

sol_path <- glmnet(x = X, y = y, weights = weights, intercept=TRUE, family ="poisson", 
               penalty.factor=penalty_factor, nlambda = nlambda, standardize = FALSE,
               lambda.min.ratio = lmr)

############################################
############ MODEL SELECTION ###############
############################################
N = v*(v-1)/2
pseq = sol_path$df

# cross-validation:
nfolds = 10
cv = cv.glmnet(x = X, y = y, weights = weights, intercept=TRUE, family ="poisson", 
               type.measure='deviance', penalty.factor=penalty_factor, nlambda = nlambda,
               standardize = FALSE, lambda.min.ratio = lmr, nfolds = nfolds)
cvsol = coef(cv, s = "lambda.min")
cvopt = sum(cvsol!=0)

# AIC:
AIC = 2*pseq + deviance(sol_path)
laic = sol_path$lambda[(1:length(AIC))[AIC==min(AIC)]]
aicsol = coef(cv, s = laic)
aicopt = sum(aicsol!=0)

# BIC:
BIC = pseq*log(N) + deviance(sol_path)
lbic = sol_path$lambda[(1:length(BIC))[BIC==min(BIC)]]
bicsol = coef(cv, s = lbic)
bicopt = sum(bicsol!=0)

# GIC in Fan and Tang (2013))
GIC = pseq*log(log(N))*log(pseq) + deviance(sol_path)
lgic = sol_path$lambda[(1:length(GIC))[GIC==min(GIC)]]
gicsol = coef(cv, s = lgic)
gicopt = sum(gicsol!=0)

# MBIC in Chand (2012)
MBIC = pseq*log(N)*sqrt(N)/pseq + deviance(sol_path)
lmbic = sol_path$lambda[(1:length(MBIC))[MBIC==min(MBIC)]]
mbicsol = coef(cv, s = lmbic)
mbicopt = sum(mbicsol!=0)

phi_vec = as.numeric(phi[lower.tri(phi)])
true_phi = (phi_vec!=0)
indexes = ((p+1):(p*(p+1)/2))
cvsubs = (cvsol[indexes] !=0)
aicsubs = (aicsol[indexes] !=0)
bicsubs = (bicsol[indexes] !=0)
gicsubs = (gicsol[indexes] !=0)
mbicsubs = (mbicsol[indexes] !=0)

accu = function(truth, est) {
  tp = sum(truth == T & est == T)
  fp = sum(truth == F & est == T)
  fn = sum(truth == T & est == F)
  tn = sum(truth == F & est == F)
  acc = (tp+tn)/(tp+fp+fn+tn)
  return(round(acc,3))
}

# compute accuracy of CV, AIC, BIC, GIC and MBIC
acc_cv = accu(true_phi, cvsubs)
acc_aic = accu(true_phi, aicsubs)
acc_bic = accu(true_phi, bicsubs)
acc_gic = accu(true_phi, gicsubs)
acc_mbic = accu(true_phi, mbicsubs)

accseq = numeric()
for (i in 1:length(sol_path$df)) {
  lam = sol_path$lambda[i]
  sol = coef(cv, s = lam)
  solsubs = (sol[indexes] !=0)
  accseq = c(accseq, accu(true_phi, solsubs))
}

# maximum achievable accuracy over the grid
acc_max = max(accseq)

#################################
########### RESULTS: ############
#################################

cat('Simulation ', simulation, ', number of vertices = ', v, '\n', 'Accuracy:', '\n',
    'Cross-validation: ', acc_cv, '\n', 'AIC: ', acc_aic, '\n', 'BIC: ', acc_bic, '\n',
    'GIC: ', acc_gic, '\n', 'MBIC: ', acc_mbic, '\n',
    'Max achievable accuracy: ', acc_max, sep = '')
