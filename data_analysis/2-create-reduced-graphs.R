##################################################
# This script allows to reproduce the reduced 
# graphs of collaborations and repulsions 
# reproduced in Figure 2 of the paper and
# in Figure 2 of the Supplementary Material.
# Use the variable case (= 1, 2, 3, 4) to select 
# the legislature of interest, as detailed below
##################################################

library(igraph)

# REMARK:
# SET case = 1 for the analysis of the 14th legislature;
# SET case = 2 for the analysis of the 15th legislature;
# SET case = 3 for the analysis of the 16th legislature;
# SET case = 4 for the analysis of the 17th legislature.
case = 4

if (case == 1) load('data\\model-XIV-legislature.RData')
if (case == 2) load('data\\model-XV-legislature.RData')
if (case == 3) load('data\\model-XVI-legislature.RData')
if (case == 4) load('data\\model-XVII-legislature.RData')

normalize = function(x) (x-min(x))/(max(x)-min(x))

alphatemp = estimates[(dim(Xcov)[2]+dim(Xregion)[2]+2):(dim(Xcov)[2]+dim(Xregion)[2]+nparties)]
phitemp = estimates[(dim(Xcov)[2]+dim(Xregion)[2]+nparties+1):(dim(X)[2]+1)]

alpha = rep(0, nparties)
alpha[1] = -sum(alphatemp)
alpha[2:nparties] = alphatemp

phi = matrix(0, nparties, nparties)
k=1
for (i in 1:(nparties-1)) {
  for (j in (i+1):nparties) {
    phi[i,j] = phitemp[k]; phi[j,i] = phi[i,j]
    k=k+1
  }
}
diag(phi) = -apply(phi,1,sum)

# display the estimates of phi_rs coefficients:
colnames(phi) = partyshort; rownames(phi) = partyshort
round(phi,2)

collaborations = phi>0
repulsions = phi<0

posgraph = graph.adjacency(collaborations, mode='undirected', weighted = TRUE)
V(posgraph)$label = partyshort
V(posgraph)$size = 15 + 15*normalize(alpha)
neggraph = graph.adjacency(repulsions, mode='undirected', weighted = TRUE)
V(neggraph)$label = partyshort
V(neggraph)$size = 15 + 15*normalize(alpha)

if (case == 1) {
  V(posgraph)$color = 'white'; V(neggraph)$color = 'white';
  V(posgraph)$color[7] = 'gray50'; V(neggraph)$color[7] = 'gray50'
  shape = c(rep('square',4),rep('circle',2),c('circle','circle'))
}

if (case == 2) {
  V(posgraph)$color = 'white'; V(neggraph)$color = 'white';
  V(posgraph)$color[6] = 'gray50'; V(neggraph)$color[6] = 'gray50'
  shape = c(rep('square',3),'circle','square',rep('circle',6),c('square','circle'))
}

if (case == 3) {
  V(posgraph)$color = 'white'; V(neggraph)$color = 'white';
  V(posgraph)$color[c(1,5,8)] = 'gray50'; V(neggraph)$color[c(1,5,8)] = 'gray50'
  shape = c('square','square','circle','square','circle','circle','square','square')
}

if (case == 4) {
  V(posgraph)$color = 'white'; V(neggraph)$color = 'white';
  V(posgraph)$color[c(5,9)] = 'gray50'; V(neggraph)$color[c(5,9)] = 'gray50';
  V(posgraph)$color[6] = 'gray75'; V(neggraph)$color[6] = 'gray75'
  shape = c(rep('square',4),rep('circle',4),'square','circle')
}

V(posgraph)$shape = shape; V(neggraph)$shape = shape

main = c('XIV legislature (2001-2006)', 'XV legislature (2006-2008)', 'XVI legislature (2008-2013)',
         'XVII legislature (2013-2015)')

V(posgraph)$label.color = 'black'
V(neggraph)$label.color = 'black'

E(posgraph)$color='black'
E(neggraph)$color='black'

if (case == 3) V(posgraph)$label[c(2,7)] = c('PDL','P&T')
if (case == 3) V(neggraph)$label[c(2,7)] = c('PDL','P&T')
if (case == 2) set.seed(3)

# REDUCED GRAPH DISPLAYING COLLABORATIONS (Figure 2):
plot(posgraph, layout=layout.fruchterman.reingold, main = main[case])

# REDUCED GRAPH DISPLAYING REPULSIONS (Figure 2 in Supplementary Material):
plot(neggraph, layout=layout.fruchterman.reingold, main = main[case])
