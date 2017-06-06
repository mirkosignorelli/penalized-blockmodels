################################################
# This script contains functions that are used
# in scripts '1_network_generation.R' and
# '2_model_selection.R', and are collected
# here for convenience. Do not run this script,
# source it at the beginning of scripts 1 and 2
# with the command: source('0_functions.R')
################################################


phi_matrix = function(p, nzeros, cmin, cmax) {
  phi_matrix = matrix(0, nrow = p, ncol = p)
  zeropos = sample(1:(p*(p-1)/2), size = nzeros, replace = F)
  k = 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      phi_matrix[i,j] = sample(x = c(-1,1), size = 1) * runif(1, min = cmin, max = cmax)
      if (k %in% zeropos) phi_matrix[i,j] = 0
      phi_matrix[j,i] = phi_matrix[i,j]
      k = k + 1
    }
  }
  diag(phi_matrix) = -apply(phi_matrix,1,sum)
  return(phi_matrix)
}

create_network = function(n, p, phi_matrix, theta0, alphar) {
  group = sample(1:p, size = n, replace = T)
  p_r = runif(p, min = -alphar, max = alphar)
  mu = matrix(NA, nrow = p, ncol = p)
  k=1
  for (i in 1:p) {
    for (j in i:p) {
      w1 = runif(1)
      linear = theta0 + p_r[i] + p_r[j] + phi_matrix[i,j]
      mu[i,j] = exp(linear)
      k=k+1
    }
  }
  minmax = function(a,b) {
    if (a<=b) out = c(a,b)
    else out = c(b,a)
    return(out)
  }
  adj = matrix(0, nrow=n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      ordgr = minmax(group[i],group[j])
      adj[i,j] = rpois(1, lambda = mu[ordgr[1],ordgr[2]])
    } 
    if (i%%1000==0) print(i)
  }
  adjM = as(adj, 'sparseMatrix')
  phi_vec = phi_matrix[lower.tri(phi_matrix, diag = T)]
  return( list('adj' = adjM, 'group' = group, 'mu' = mu, 'phi_vec' = phi_vec) )
}

colmatr = function(p) {
  out = matrix(NA, nrow=p, ncol = p)
  k=1
  for (i in 1:p) {
    for (j in i:p) {
      out[i,j] = k
      k = k+1
    }
  }
  return(out)
}

dvec = function(gri, grj, p) {
  dcolumn = colmatr(p)
  dvec = rep(0, p + p*(p+1)/2)
  if (gri != grj) {
    dvec[gri] = 1
    dvec[grj] = 1
  }
  else if (gri == grj) dvec[gri] = 2
  if (gri <= grj) dvec[p + dcolumn[gri,grj]] = 1
  else dvec[p + dcolumn[grj,gri]] = 1
  return(dvec)
}

dtrasfvec = function(gri, grj, p) {
  dcolumn = colmatr(p)
  dvec = dvec(gri, grj, p)
  dtrasfvec = rep(0, length = p-1 + p*(p-1)/2)
  dtrasfvec[1:(p-1)] = dvec[2:p] - dvec[1]
  k = p
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      dtrasfvec[k] = dvec[p+dcolumn[i,j]] - dvec[p+dcolumn[i,i]] - dvec[p+dcolumn[j,j]]
      k = k+1
    }
  }
  return(dtrasfvec)
}
