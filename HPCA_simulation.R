HPCA_simulation <- function(c,       # tuning parameter for signal-to-noise ratio (scalar)  
                            N,       # group sample size (scalar)
                            sparse,  # 1 produces longitudinal sparsity, 0 otherwise (scalar)
                            ll,      # minimum number of longitudinal observations (scalar)
                            ul       # maximum number of longitudinal observations (scalar)
                           ){
  
#############################################################################
## Description: Function for simulating data for the HPCA_decomp function as described
##              in the supplementary materials.
## Args: see above
## Returns: data, data.frame with columns c("ID", "group" , "y", "reg", "long", "func")
#############################################################################  

#############################
# 1. Define model components
#############################

# Domain
reg = c(1:9) # Regional domain
reg_vec = seq(from=0, to = 1, length.out=9)  # Regional domain for construction of eigenvectors
func = seq(from = 0, to = 1, length.out = 50)  # Functional domain
long = seq(from = 0, to = 1, length.out = 50)  # Longitudinal domain
func_long = expand.grid(func, long)  # Longitudinal-funcitonal domain

# Global variables
K = 3; L = 3; M = 3  # Number of eigencomponents
D = 2  # Number of groups

# Fixed Effects
mu_vec = 5 * sqrt(1 - (func_long[, 1] - .5)^2 - (func_long[, 2] - .5)^2)  # Overall mean function 

eta_vec = matrix(list(), nrow = D)  # Group-region-specific shifts
for(d in 1:D){
  eta_vec[[d]] = matrix(list(), nrow = length(reg))
  for(r in reg){
    eta_vec[[d]][[r]] = (- 1)^d * (sqrt(3) * func_long[, 1])^2
  }
}

# Eigenfunctions
# Reginonal marginal eigenvectors
v_k = matrix(list(), nrow = 3)
v_k[[1]] = sin(1 * reg_vec * pi) / 2
v_k[[2]] = sin(2 * reg_vec * pi) / 2
v_k[[3]] = sin(3 * reg_vec * pi) / 2

# Functional marginal eigenfunctions
phi_l = matrix(list(), nrow = 3)
phi_l[[1]] = sqrt(2) * sin(1 * func * pi)
phi_l[[2]] = sqrt(2) * sin(2 * func * pi)
phi_l[[3]] = sqrt(2) * sin(3 * func * pi)

# Longitudinal marginal eigenfunctions
psi_m = matrix(list(), nrow = 3)
psi_m[[1]] = sqrt(2) * cos(1 * long * pi)
psi_m[[2]] = sqrt(2) * cos(2 * long * pi)
psi_m[[3]] = sqrt(2) * cos(3 * long * pi)


# Multidimensional orthonormal basis
eigenfunc_vec = matrix(list(), nrow = L)
for(l in 1:L){
  eigenfunc_vec[[l]] = matrix(list(), nrow = M)
  for(m in 1:M){
    eigenfunc_vec[[l]][[m]] = kronecker(psi_m[[m]], phi_l[[l]])  # Form basis functions for longitudinal/functional domain
  }
}

# Variance components
tau_klm = matrix(list(), nrow = 3)
for(k in 1:K){
  tau_klm[[k]] = matrix(list(), nrow = 3)
  for(l in 1:L){
    tau_klm[[k]][[l]] = matrix(list(), nrow = 3)   
    for(m in 1:M){
      tau_klm[[k]][[l]][[m]] = sqrt(1 / (k * l * m))  # Define variance components for subject-specific scores     
    }
  }
}

# Measurement error variance 
sig_eps = 5/c

#############################
# 2. Simulate data
#############################

subj_length = length(func) * length(long) * length(reg)

# Create data.frame to store observations
data = data.frame(matrix(NA, nrow=subj_length * D * N, ncol=6))
colnames(data) = c("ID", "group" , "y", "reg", "long", "func")
data$group = rep(c(1:D), each = subj_length * N)
data$ID = rep(c(1:(D * N)), each = subj_length)
data$reg = rep(rep(reg,each=length(func) * length(long)), D * N)
data$func = rep(rep(rep(func,length(long)), length(reg)), D * N)
data$long = rep(rep(rep(long,each=length(func)), length(reg)), D * N)

# Simulate subject-specific scores
scores<-matrix(list(), nrow = D)
for(d in 1:D){
  scores[[d]] = matrix(list(), nrow = N)
  for(i in 1:N){
    scores[[d]][[i]] = matrix(list(), nrow = K)
    for(k in 1:K){
      scores[[d]][[i]][[k]] = matrix(list(), nrow = L)
      for(l in 1:L){
        scores[[d]][[i]][[k]][[l]] = matrix(list(), nrow = M)
        for(m in 1:M){
          scores[[d]][[i]][[k]][[l]][[m]] = rnorm(1, 0, sqrt(tau_klm[[k]][[l]][[m]]))
        }
      }
    }
  }
}

# Simulate subject-specific trajectories
for(d in 1:D){
  for(i in 1:N){
    for(r in 1:length(reg)){
      subj_traj = mu_vec + eta_vec[[d]][[r]]
      for(k in 1:K){
        for(l in 1:L){
          for(m in 1:M){
            subj_traj = subj_traj + scores[[d]][[i]][[k]][[l]][[m]] * eigenfunc_vec[[l]][[m]] * v_k[[k]][r]
          }
        }
      }
      ind = 1 + (d - 1) * N * subj_length + (i - 1) * subj_length + (r -1) * length(func) * length(func)
      data[ind:(ind + length(func) * length(long) - 1), 3] = subj_traj + rnorm(length(func) * length(long), 0, sqrt(sig_eps))
    }
  }
} 

# Introduce longitudinal sparsity

if(sparse == 1){
  # Introduce sparsity in the longitudinal dimension 
  long_logical = matrix(1, nrow = 2 * N, ncol = 50)
  long_num = sample(ll:ul, 2 * N , replace = T)  # Simulate number of longitudinal functional observations
  
  for(i in 1:length(long_num)){
    if(long_num[i] < 50){
      long_logical[i, ] = sample(cbind(t(rep(0, long_num[i])), t(rep(1, 50 - long_num[i]))))
    }else{long_logical[i, ] = rep(0, 50)}
  }
  
  for(i in unique(data$ID)){  # Remove unobserved longitudinal functional observations
    if(sum(long_logical[i, ]) > 0){
      data = data[-which(data$ID == i & data$long %in% long[as.logical(long_logical[i, ])]), ]
    }
  }
}

return(data)

}