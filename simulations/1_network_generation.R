##############################################
# This script allows to generate the networks
# considered in simulations A, B, C and D.
# In order to reproduce the results reported
# in the Supplementary Material, do not alter 
# the seed (123) and the order in which
# the code appears. 
##############################################

source('0_functions.R')
library(Matrix)

p = 10
n = 5000
set.seed(123)

##### SIMULATION B
nzeros = 20
phi = phi_matrix(p, nzeros, 0.2, 0.5)
netobj = create_network(n, p, phi, theta0 = 0.7, alphar = 0.3) 
save.image('simB.RData')

##### SIMULATION A
nzeros = 10
phi = phi_matrix(p, nzeros, 0.2, 0.5)
netobj = create_network(n, p, phi, theta0 = 0.7, alphar = 0.3) 
save.image('simA.RData')

##### SIMULATION C
nzeros = 30
phi = phi_matrix(p, nzeros, 0.2, 0.5)
netobj = create_network(n, p, phi, theta0 = 0.7, alphar = 0.3) 
save.image('simC.RData')

##### SIMULATION D
nzeros = 20
phi = phi_matrix(p, nzeros, 0.1, 0.5)
netobj = create_network(n, p, phi, theta0 = 0.7, alphar = 0.3) 
save.image('simD.RData')
