HPCA_bootstrap <- function(B,       # number of bootstrap samples (scalar)
                           func,    # functional domain grid (vector)
                           long,    # longitudinal domain grid (vector)
                           reg,     # regional domain grid (vector)
                           HPCA,    # HPCA output from HPCA_decomp.R
                           region,  # region to consider for hypothesis test (scalar)
                           D,       # number of groups (scalar)
                           N,       # group sample sizes (vector)
                           K,       # number of leading regional marginal eigenvectors used to generate outcomes  (scalar)
                           L,       # number of leading functional marginal eigenfunctions used to generate outcomes  (scalar)
                           M        # number of leading longitudinal marginal eigenfunctions used to generate outcomes  (scalar)
                          ){
#############################################################################
## Description: Function for performing group-level inference for a given
##              region as described in "Hybrid Principal Components Analysis For  
##              Region-Referenced Longitudinal Functional EEG Data" by Scheffler et al. (2017).
##              The null hypothesis assumes that groups share a common region shift for
##              a given region.
## Args:        see above
## Returns:     p-value: A region-specific p-value for the observed test 
##              statistic under the null hypothesis.
##               
#############################################################################    
  
# Install missing packages
list.of.packages = c("refund", "caTools")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load libraries
library(refund)
library(caTools)
  
# Numerical integration function
num_int <- function(x){
  int = trapz(func, x)
  return(int)
}

# Indexing function
eig_index <- function(x){
  ind = (x[1] - 1) + K * (x[2] - 1) + K * L * (x[3] - 1) + 1;
  return(ind)
}
#############################################################################
# 0. Setup Parametric Bootstrap Components
#############################################################################

# Define global variables
f_tot = length(func)  # Total grid points in functional domain
l_tot = length(long)  # Total grid points in longitudinal domain
r_tot = length(reg)  # Total grid points in regional domain
 
# Load fixed effects from HPCA decomposition results
mu = HPCA[[1]]
eta = HPCA[[2]]
  
# Obtain region shift under null hypothesis
eta_null = matrix(0, nrow = f_tot, ncol = l_tot)
for(d in 1:D){
  eta_null = eta_null + eta[[d]][[region]] / D
}
  
for(d in 1:D){
  for(r in region:region){
      eta[[d]][[r]] = eta_null
  }
}
  
# Load eigencomponents from HPCA decomposition results
phi = matrix(list(), nrow = D)
psi = matrix(list(), nrow = D)
nu = matrix(list(), nrow = D)
  
for(d in 1:D){
  nu[[d]]<-HPCA[[5]][[d]][[2]]$vectors
  phi[[d]]<-HPCA[[3]][[d]][[2]]$vectors
  psi[[d]]<-HPCA[[4]][[d]][[2]]$vectors
}
  
# Random effects error variance
lambda = matrix(list(), nrow = D)
for(d in 1:D){
  lambda[[d]]  = diag(HPCA[[6]][[d]][[2]])
}  

# Measurement error variance
sigma = rep(NA, D)
for(d in 1:D){
  sigma[d]  =  HPCA[[6]][[d]][[3]]
}  

print("0. Setup Parametric Bootstrap Components (completed)")
  
#############################################################################
# 1. Bootstrap Procedure 
#############################################################################

# Create list to store group-region shifts for each bootstrap iteration
eta_b = matrix(list(), nrow = D)
for(d in 1:D){
  eta_b[[d]] = matrix(list(), nrow = B)
}

# Create list to store region shifts under the null for each bootstrap iteration
eta_b_null = matrix(list(), nrow = B)

# Parametric bootstrap procedure  
for(b in 1:B){
  # Generate data
  y = matrix(list(), nrow = D)
  for(d in 1:D){
    y[[d]] = matrix(list(), nrow = N[d])
    for(i in 1:N[d]){
      y[[d]][[i]] = matrix(list(), nrow = r_tot)
      for(r in 1:r_tot){
        noise = matrix(rnorm(f_tot * l_tot, 0, sqrt(sigma[[d]])), nrow = f_tot, ncol = l_tot)
        y[[d]][[i]][[r]] = mu + eta[[d]][[r]] + noise
        for(k in 1:K){
          for(l in 1:L){
            for(m in 1:M){
              RE_sim = rnorm(1, 0, sqrt(lambda[[d]][eig_index(c(k, l, m))]))
              y[[d]][[i]][[r]] = y[[d]][[i]][[r]] + RE_sim * nu[[d]][r, k] * matrix(kronecker(phi[[d]][, l],
                                 psi[[d]][, m]), nrow = f_tot, ncol = l_tot)
              }
            }
          }
        }
      }
    }
  
    # Estimate mean surface
    obs_overall = matrix(NA, nrow = f_tot, ncol = l_tot * r_tot * sum(N))
    j=1
    for(d in 1:D){
      for(i in 1:N[d]){
        for(r in 1:r_tot){
          ll = 1 + (j - 1) * f_tot
          ul = j * f_tot
          obs_overall[, c(ll:ul)] = y[[d]][[i]][[r]]
          j = j + 1
        }
      }
    }
    
    long_obs = rep(long, r_tot * sum(N))
    mean_surf = fbps(obs_overall, covariates = list(x = t(func), z = long_obs)) 
    newdata = list(x = rep(func, l_tot), z = rep(long, each = f_tot))
    prediction = predict(mean_surf, newdata = newdata)
    mu_bs = matrix(prediction$fitted.values, nrow = f_tot)
    
    # Estimate region shifts 
    for(d in 1 :D){
      obs_overall<-matrix(NA, nrow=50, ncol = l_tot * N[d])
      j = 1
      for(i in 1:N[d]){
        for(r in region:region){
          ll = 1 + (j - 1) * f_tot
          ul = j * f_tot
          obs_overall[, c(ll:ul)] = y[[d]][[i]][[r]] - mu_bs
          j = j + 1
        }
      }
      long_obs = rep(long, N[d])
      eta_surf = fbps(obs_overall, covariates = list(x = func, z = long_obs)) #3.151 seconds
      newdata = list(x = rep(func, l_tot), z = rep(long, each = f_tot))
      prediction = predict(eta_surf, newdata = newdata)
      eta_bs = matrix(prediction$fitted.values, nrow = f_tot)
      eta_b[[d]][[b]] = eta_bs
      rm(eta_bs)
      rm(obs_overall)
    }
    
    # Calculate region shift under the null hypothesis
    eta_b_null[[b]] = matrix(0, nrow = f_tot, ncol= l_tot)
    for(d in 1:D){
      eta_b_null[[b]] = eta_b[[d]][[b]] / D + eta_b_null[[b]]
    }
    
    #end bootstrap
    print(paste("Bootstrap...", as.character(b)))
    
}


#############################################################################
# 3. Calculate p-values   
#############################################################################  
# Find distribution of test statistic under the null hypothesis
dist <- matrix(NA, nrow = B)
for(b in 1:B){
  dist[b] = 0
  for(d in 1:D){
  dist[b] = trapz(long, apply((eta_b[[d]][[b]] - eta_b_null[[b]])^2, 1, num_int)) + dist[b]
  }
  dist[b] = sqrt(dist[b])
}

# Calculate test statistic for observed data
ts = 0
for(d in 1:D){
  ts = trapz(long,apply((HPCA[[2]][[d]][[region]] - eta[[d]][[region]])^2, 1, num_int)) + ts
}
ts = sqrt(ts)

# Obtain p-value
p_val = sum(dist>ts) / B

return(p_val)
}