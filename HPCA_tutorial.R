## HPCA_tutorial.R
#############################################################################
## Description: A step-by-step implementation of HPCA_decomposition.R and the associated  
## procedures described in "Hybrid Principal Components Analysis For Region-Referenced Longitudinal 
## Functional EEG Data" by Scheffler et al. (2017). The procedure assumes that data is densely
## observed in the regional and functional dimensions and observed either with or without
## sparsity in the longitudinal dimension.
#############################################################################
## Functions implemented: 
## HPCA_decomposition.R, HPCA_simulation.R, HPCA_bootstrap.R
#############################################################################
## Tutorial Outline:
## 1. Simulate region-referenced longitudinal functional data (HPCA_simulation.R)
## 2. Perform HPCA decomposition (HPCA_decomposition.R)
## 3. Likelihood-ratio test on the correlation structure of subject-specific scores
## 4. Group-level inference via bootstrap (HPCA_bootstrap.R)
## 5. Visualization of HPCA decomposition results
#############################################################################

# Install missing packages
list.of.packages <- c("refund", "data.table", "fANCOVA", "mgcv", "Matrix", "caTools", "nlme", "dplyr", "glmmTMB")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
    
# Load packages  
library(refund)
library(data.table)
library(fANCOVA)
library(mgcv)
library(Matrix)
library(caTools)  
library(nlme)
library(dplyr) 


#############################################################################
# 1. Simulate region-referenced longitudinal functional data
#############################################################################

# Simulate data from the dense simulation design at nd = 15 and high SNR 
data = HPCA_simulation(c = 10, N = 15, sparse = 0 , ll = 20, ul = 50)  # HPCA_simulation.R

#############################################################################
# 2. Perform HPCA decomposition
#############################################################################

# NOTE: Performing the HPCA decomposition will take approximately seven minutes.

# Domain
reg = c(1:9) # Regional domain
func = seq(from = 0, to = 1, length.out = 50)  # Functional domain
long = seq(from = 0, to = 1, length.out = 50)  # Longitudinal domain

# Perform HPCA decomposition 
HPCA = HPCA_decomp(data = data, # HPCA_decomposition.R
                          out1 = func, 
                          out2 = long, 
                          out3 = reg,
                          kFE = c(35, 35),
                          kCOV = 10,
                          reControl = lmeControl(returnObject = TRUE, opt="nlminb", 
                                      niterEM = 25, minAbsParApVar = .0001)
                          ) 

#############################################################################
# 3. Likelihood-ratio test on the correlation structure of subject-specific scores
#############################################################################

# NOTE: Performing the likelihood-ratio test may take several hours
#       depending on the size of the data set and processor speed. 

# Load glmmTMB for fitting linear mixed effects models 
library(glmmTMB)

# Specify linear mixed effects models
re_mat = paste("x", 1:27, sep="")
cs_form = as.formula(paste('y_c2 ~', '-1+ cs(0 +', paste(re_mat, collapse= "+"), '|ID)'))
diag_form = as.formula(paste('y_c2 ~', '-1+ diag(0 +', paste(re_mat, collapse= "+"), '|ID)'))

# Fit linear mixed effects models for group 1
fit.cs = glmmTMB(cs_form, # Heterogeneous compound symmetry structure
                 data = HPCA[[7]][which(HPCA[[7]]$group == 1), ])  
fit.diag = glmmTMB(diag_form, # Diagonal covariance strucure
                   data = HPCA[[7]][which(HPCA[[7]]$group == 1), ])  

# Perform likelihood-ratio test, 
csq_group_1 = 2 * (logLik(fit.cs) - logLik(fit.diag))
covar_p = (1 - pchisq(csq_group_1[1], df = 1))[1] # p-value for likelihood-ratio test

#############################################################################
# 4. Group-level inference via bootstrap 
#############################################################################

# NOTE: Fitting the bootstrap procedure will take approximately one hour. 

# Produce p-value from bootstrap procedure
boot_p = HPCA_bootstrap(B = 200, # HPCA_bootstrap.R
                        func = func, 
                        long = long, 
                        reg = reg, 
                        HPCA, 
                        region = 5, 
                        D = c(1, 2), 
                        N = c(15, 15), 
                        kFE = c(35, 35))  

#############################################################################
# 5. Visualization of HPCA decomposition results
#############################################################################  

# Install missing packages
list.of.packages <- c("refund", "caTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries 
library(refund) 
library(caTools)

eig_index<-function(x){ # Indexing function
  ind<-(x[1]-1)+3*(x[2]-1)+3*3*(x[3]-1)+1;
  return(ind)
}

##############################
# True model components from simulated data

reg_vec = seq(from=0, to = 1, length.out=9)  # Regional domain for construction of eigenvectors
func_long = expand.grid(func, long)  # Longitudinal-functional domain

# Overall mean function 
mu = matrix(5 * sqrt(1 - (func_long[, 1] - .5)^2 - (func_long[, 2] - .5)^2), nrow = 50)  

# Group-region shifts
eta = matrix(list(), nrow = 2)  
for(d in 1:2){
  eta[[d]] = matrix(list(), nrow = length(reg))
  for(r in reg){
    eta[[d]][[r]] = matrix((- 1)^d * (sqrt(3) * func_long[, 1])^2, nrow = 50)
  }
}

# Regional marginal eigenvectors
v_k = matrix(list(), nrow = 3)
for(k in 1:3){
  v_k[[k]] = sin(k * reg_vec * pi) / 2
}

# Functional marginal eigenfunctions
phi_l = matrix(list(), nrow = 3)
for(l in 1:3){
  phi_l[[l]] = sqrt(2) * sin(l * func * pi)
}

# Longitudinal marginal eigenfunctions
psi_m = matrix(list(), nrow = 3)
for(m in 1:3){
  psi_m[[m]] = sqrt(2) * cos(m * long * pi)
}

# Store estimated model components (group = 1) 
mu_hat = HPCA[[1]]  # Overall mean function
eta_hat = HPCA[[2]]  # Group-region-specific shifts
v_k_hat = HPCA[[5]][[1]][[2]]$vectors  # Regional marginal eigenvectors
phi_l_hat = HPCA[[3]][[1]][[2]]$vectors  # Functional marginal eigenfunctions
psi_m_hat = HPCA[[4]][[1]][[2]]$vectors  # Longitudinal marginal eigenfunctions

# Align eigenfunctions
for(i in 1:3){
  if(trapz(func, (phi_l_hat[, i] - phi_l[[i]])^2) > 
     trapz(func, ( -phi_l_hat[, i] - phi_l[[i]])^2)){
    phi_l_hat[, i] = -phi_l_hat[, i]
    }
  if(trapz(func, (psi_m_hat[, i] - psi_m[[i]])^2) > 
     trapz(func, ( -psi_m_hat[, i] - psi_m[[i]])^2)){
    psi_m_hat[, i] = -psi_m_hat[, i]
    }  
}

##############################
# Visualize model components

# Plot overall mean function
par(mfrow = c(1, 2))
persp(mu, col = "blue", theta = -45, expand = 0.8, 
      r = 4 , main = "True mean function",
      xlab = "functional", ylab = "longitudinal", zlab = "")
persp(mu_hat, col = "red", theta = -45,  
      expand = 0.8, r = 4 , main = "Estimated mean function",
      xlab = "functional", ylab = "longitudinal", zlab = "")

# Plot group-region shift (group = 1, region = 1)
par(mfrow = c(1, 2))
persp(eta[[1]][[5]], col = "blue", theta = 45,  expand = 0.8, 
      r = 2 , main = "True group-region shift",
      xlab = "functional", ylab = "longitudinal", zlab = "")
persp(eta_hat[[1]][[5]], col = "red", theta = 45, 
      expand = 0.8, r = 2 , main = "Estimated group-region shift",
       xlab = "functional", ylab = "longitudinal", zlab = "")

##############################
# Plot functional marginal eigenfunctions (group = 1)
par(mfrow = c(1, 3))
# phi_1
plot(func, phi_l_hat[, 1], type = "l", ylab = expression(phi[d1](omega)), 
     xlab =  expression(paste("functional ", (omega))), 
     lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, phi_l[[1]], type = "l")
legend(0, -1.5, c("True", "Estimated"), lty=c(1,2))
# phi_2
plot(func, phi_l_hat[, 2], type = "l", ylab = expression(phi[d2](omega)), 
     xlab = expression(paste("functional ", (omega))), 
     lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, phi_l[[2]], type = "l")
# phi_3
plot(func, phi_l_hat[, 3], type = "l", ylab = expression(phi[d3](omega)), 
     xlab =  expression(paste("functional ", (omega))),
     lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, phi_l[[3]], type = "l")
title("Functional Marginal Eigenfunctions", outer=TRUE, line = - 2)


##############################
# Plot longitudinal marginal eigenfunctions (group = 1)
par(mfrow = c(1, 3))
# psi_1
plot(func, psi_m_hat[, 1], type = "l", ylab = expression(psi[d1](s)), 
     xlab = expression(paste("longitudinal ", (s))), lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, psi_m[[1]], type = "l")
legend(0, -1.5, c("True", "Estimated"), lty=c(1,2))
# psi_2
plot(func, psi_m_hat[, 2], type = "l", ylab = expression(psi[d2](s)), 
     xlab = expression(paste("longitudinal ", (s))), lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, psi_m[[2]], type = "l")
# psi_3
plot(func, psi_m_hat[, 3], type = "l", ylab = expression(psi[d3](s)),
     xlab = expression(paste("longitudinal ", (s))), lty = 2, ylim = c(-2, 2), xaxs="i")
lines(func, psi_m[[3]], type = "l")
title("Longitudinal Marginal Eigenfunctions", outer=TRUE, line = - 2)

##############################
# Subject-level trajectories  (group = 1, ID = 1, region = 1)

# Observed subject-level trajectory
y_obs = matrix(data[which((data$group == 1) & (data$ID == 1) & (data$reg == 1)), 3], 
               nrow = 50)

# Predicted subject-level trajectory
y_pred = HPCA[[1]] + HPCA[[2]][[1]][[1]]                        
for(k in 1:3){
  for(l in 1:3){
    for(m in 1:3){
      y_pred = y_pred + HPCA[[6]][[1]][[1]][5, eig_index(c(k,l,m))] * 
        v_k[[k]][1] * matrix(kronecker(psi_m[[m]], phi_l[[l]]), nrow = 50)
    }
  }
}
 
# Plot a subject-specific trajectory
par(mfrow = c(1, 2))
persp(y_obs, col = "blue", theta = 30, zlim = c(-1, 7),
      expand = 0.8, r = 4 , main = "Observed Subject-specific Trajectory",
      xlab = "functional", ylab = "longitudinal", zlab = "")
persp(y_pred, col = "red", theta = 30, zlim = c(-1, 7), 
      expand = 0.8, r = 4 , main = "Predicted Subject-specific Trajectory",
      xlab = "functional", ylab = "longitudinal", zlab = "")
