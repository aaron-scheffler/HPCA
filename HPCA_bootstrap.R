HPCA_bootstrap <- function(B,       # number of bootstrap samples (scalar)
                           func,    # functional domain grid (vector)
                           long,    # longitudinal domain grid (vector)
                           reg,     # regional domain grid (vector)
                           HPCA,    # HPCA output from HPCA_decomp.R
                           region,  # region to consider for hypothesis test (scalar)
                           D,       # list of groups to include in bootstrap (vector)
                           N,       # group sample sizes (vector),
                           kFE      # knots for smoothing of fixed effects (vector)
                          ){
#############################################################################
## Description: Function for performing group-level inference for a given
##              region as described in "Hybrid Principal Components Analysis For  
##              Region-Referenced Longitudinal Functional EEG Data" by Scheffler et al. (2017).
##              The null hypothesis assumes that groups share a common region shift for
##              a given region.
## Args:        (see above)
## Returns:     p_boot: A region-specific p-value for the observed test 
##              statistic under the null hypothesis.
##              eta_b: The common region shift under the null hypothesis from
##              each run of the bootstrap procedure.
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
  int = trapz(long, x)
  return(int)
}

# Indexing function
eig_index <- function(x){
  ind = (x[1] - 1) + K[d] * (x[2] - 1) + K[d] * L[d] * (x[3] - 1) + 1;
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
eta_null <- matrix(list(), nrow = length(D))

for(d in D){ # storage
  eta_null[[d]] <- matrix(list(), nrow = r_tot)
  for(r in 1:r_tot){
    eta_null[[d]][[r]] <- matrix(0, nrow = f_tot, ncol = l_tot)
  }
}

for(d in D){
  for(r in 1:r_tot){
    if(r %in% region){
      for(d1 in D){
        eta_null[[d]][[r]] = eta[[d1]][[r]] / length(D) + eta_null[[d]][[r]] 
      }
    } else{
      eta_null[[d]][[r]] <- eta[[d]][[r]]
      }
  }
}

# Load eigencomponents from HPCA decomposition results
phi = matrix(list(), nrow = length(D))
psi = matrix(list(), nrow = length(D))
nu = matrix(list(), nrow = length(D))
  
for(d in D){
  nu[[d]] <- HPCA[[5]][[d]][[2]]$vectors
  phi[[d]] <- HPCA[[3]][[d]][[2]]$vectors
  psi[[d]] <- HPCA[[4]][[d]][[2]]$vectors
}
 
# Exctract number of eigencomponents to include
K <- rep(NA, length(D))    
L <- rep(NA, length(D)) 
M <- rep(NA, length(D)) 

for(d in D){
  K[d] <- HPCA[[5]][[d]][[2]]$FVE
  L[d] <- HPCA[[3]][[d]][[2]]$FVE
  M[d] <- HPCA[[4]][[d]][[2]]$FVE
}  

# Random effects error variance
lambda = matrix(list(), nrow = length(D))
for(d in D){
  lambda[[d]]  = diag(HPCA[[6]][[d]][[2]])
}  

# Measurement error variance
sigma = rep(NA, length(D))
for(d in D){
  sigma[d] = HPCA[[6]][[d]][[3]]
}  

print("0. Setup Parametric Bootstrap Components (completed)")
  
#############################################################################
# 1. Bootstrap Procedure 
#############################################################################

# Create list to store group-region shifts for each bootstrap iteration
eta_b = matrix(list(), nrow = length(D))
for(d in D){
  eta_b[[d]] = matrix(list(), nrow = B)
}

# Create list to store region shifts under the null for each bootstrap iteration
eta_b_null = matrix(list(), nrow = B)

# Parametric bootstrap procedure  
for(b in 1:B){
 
   # Create storage for group-region specific shifts and null
   for(d in D){
    eta_b[[d]][[b]] = matrix(list(), nrow = r_tot)
   }
   
   eta_b_null[[b]] = matrix(list(), nrow = r_tot)
    
  # Generate data
  y = matrix(list(), nrow = length(D))
  for(d in D){
    y[[d]] = matrix(list(), nrow = N[d])
    for(i in 1:N[d]){
      y[[d]][[i]] = matrix(list(), nrow = r_tot)
      for(r in 1:r_tot){
        noise = matrix(rnorm(f_tot * l_tot, 0, sqrt(sigma[[d]])), nrow = f_tot, ncol = l_tot)
        y[[d]][[i]][[r]] = mu + eta_null[[d]][[r]] + noise
        for(k in 1:K[d]){
          for(l in 1:L[d]){
            for(m in 1:M[d]){
              RE_sim = rnorm(1, 0, sqrt(lambda[[d]][eig_index(c(k, l, m))]))
              y[[d]][[i]][[r]] = y[[d]][[i]][[r]] + RE_sim * nu[[d]][r, k] * matrix(kronecker(psi[[d]][, m],
                                 phi[[d]][, l]), nrow = f_tot, ncol = l_tot)
              }
            }
          }
        }
      }
    }
  
    # Estimate mean surface
    obs_overall = matrix(list(), nrow = length(D))
    mu_bs <- matrix(0, nrow = f_tot, ncol = l_tot)
    
    for(d in D){
      obs_overall[[d]] <- matrix(NA, nrow = f_tot, ncol = l_tot * r_tot * sum(N[d]))
      j=1 # initialize indexing
      for(i in 1:N[d]){
        for(r in 1:r_tot){
          ll = 1 + (j - 1) * l_tot
          ul = j * l_tot
          obs_overall[[d]][, c(ll:ul)] = y[[d]][[i]][[r]]
          j = j + 1
        }
      }
    long_obs = rep(long, r_tot * sum(N[d]))
    mean_surf = fbps(obs_overall[[d]], covariates = list(x = t(func), z = long_obs), knots = kFE, m = 2, p = 3) 
    newdata = list(x = rep(func, l_tot), z = rep(long, each = f_tot))
    prediction = predict(mean_surf, newdata = newdata)
    mu_bs = mu_bs +  matrix(prediction$fitted.values, nrow = f_tot)/length(D)
    }
    
    
    # Estimate region shifts 
    for(d in D){
      for(r in region){
        obs_overall<-matrix(NA, nrow=50, ncol = l_tot * N[d])
        j = 1
        for(i in 1:N[d]){
            ll = 1 + (j - 1) * l_tot
            ul = j * l_tot
            obs_overall[, c(ll:ul)] = y[[d]][[i]][[r]] - mu_bs
            j = j + 1
        }
        long_obs = rep(long, N[d])
        eta_surf = fbps(obs_overall, covariates = list(x = func, z = long_obs), knots = kFE, m = 2, p = 3)  
        newdata = list(x = rep(func, l_tot), z = rep(long, each = f_tot))
        prediction = predict(eta_surf, newdata = newdata)
        eta_b[[d]][[b]][[r]] = matrix(prediction$fitted.values, nrow = f_tot)
        rm(obs_overall)
      }
    }
    
    # Calculate region shift under the null hypothesis
    for(r in region){
      eta_b_null[[b]][[r]] = matrix(0, nrow = f_tot, ncol= l_tot)
      for(d in D){
        eta_b_null[[b]][[r]] = eta_b[[d]][[b]][[r]] / length(D) + eta_b_null[[b]][[r]]
      }
    }
    
    #end bootstrap
    print(paste("Bootstrap...", as.character(b)))
    
}


#############################################################################
# 3. Calculate p-values   
#############################################################################  


# Find distribution of test statistic under the null hypothesis
dist <- matrix(list(), nrow = B)
for(b in 1:B){
  dist[[b]] <- matrix(0, nrow = r_tot)
  for(r in region){
    for(d in D){
     dist[[b]][[r]] = trapz(func, apply((eta_b[[d]][[b]][[r]] - eta_b_null[[b]][[r]])^2, 1, num_int)) + dist[[b]][r]
    }
  }
  dist[[b]]  = sqrt(dist[[b]])
}
    

# Calculate test statistic for observed data
ts <- matrix(0, nrow = r_tot)
for(r in region){
   for(d in D){
    ts[r] <- trapz(func, apply((HPCA[[2]][[d]][[r]] - eta_null[[d]][[r]])^2, 1, num_int)) + ts[r]
  }
}
ts <- sqrt(ts)


# Obtain p-value
dist_mat <- matrix(unlist(dist), nrow = length(reg))
p_val <- matrix(NA, nrow = r_tot)
for(r in region){
    p_val[r] <- sum(dist_mat[r, ] > ts[r]) / B
}
 

#############################################################################
# 4. Output results  
############################################################################# 
p_boot <- matrix(unlist(p_val), nrow = length(reg))[region]
boot <- list(p_boot, eta_b) # return a list with the p-value and common region shift under
                            # the null from each run of the bootstrap procedure
return(boot)
}