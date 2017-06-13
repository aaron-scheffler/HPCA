HPCA_decomposition <- function(data,  # data.frame in long format with six labeled columns (described below)
                                      # and row length equal to the length of the vectorized region-referenced
                                      # observations across all subjects and groups
                                          # DATA.FRAME COLUMNS:
                                          # ID: subject IDs (vector)
                                          # group: subject group (vector)
                                          # y: region-referenced longitudinal functional data (vector)
                                          # reg: regional argument (vector)
                                          # long: longitudinal argument (vector)
                                          # func: functional argument (vector)
                              out1,  # functional domain grid (vector)
                              out2,  # longitudinal domain grid (vector)
                              out3   # regional domain grid (vector)
                        ){
#############################################################################
## Description: Function for performing HPCA decomposition described in "Hybrid Principal Components 
##              Analysis For Region-Referenced Longitudinal Functional EEG Data" by Scheffler
##              et al. (2017), including estimation of fixed effects, marginal covariance functions,
##              marginal eigencomponents, subject-specific scores, variance components, and measurement 
##              error variance. The HPCA estimation procedure mirrors the algorithm described Appendix B
##              of the aforementioned paper. 
## Args:        see above
## Returns:     list()
##              mu: overall mean function (matrix)
##              eta: group-region-specific shifts (list) 
##              func marg cov: functional marginal covariance function and eigencomponents (list)
##              long marg cov: longitudinal marginal covariance function and eigencomponents (list)
##              reg marg cov:  regional marginal covariance function and eigencomponents (list)
##              RE model: subject-specific scores, variance components, 
##              and measurement error estimated from linear mixed effects model (list)
##              data: data.table in HPCA_decomp.R (data.frame)
## HPCA Decomposition Outline:
##              0. Format data and create return list  
##              1. Estimation of Fixed Effects  
##              2. Estimation of Marginal Covariances and Measurement Error Variance 
##              3. Estimation of Marginal Eigencomponents 
##              4. Estimation of Variance Components and Subject-specific Scores  
#############################################################################

# Install missing packages
list.of.packages <- c("refund", "data.table", "fANCOVA", "mgcv", "Matrix", "caTools", "nlme", "dplyr")
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
# 0. Format data and create return list
#############################################################################

# Format data and sort  
data = data.table(data)  # Convert data matrix to data.table format 
data = data[order(group, ID, reg, long, func)]  # Sort data.table 

# Define global variables
id = unique(data$ID)  # subject ID's
groups = unique(data$group)  # Groups
g_tot = length(groups)  # Total groups
f_tot = length(out1)  # Total grid points in functional domain
l_tot = length(out2)  # Total grid points in longitudinal domain
r_tot = length(out3)  # Total grid points in regional domain

# Create dimension indices for merging 
data[, func_ind := (match(func, out1))]  # Functional dimension index
data[, long_ind := (match(long, out2))]  # Longitudinal dimension index
data[, reg_ind := (match(reg, out3))]  # Regional dimension index 
data[, index := ((reg_ind - 1) + r_tot * (func_ind - 1) + r_tot * f_tot * (long_ind - 1) + 1)]  # Regional, functional, and longitudinal index
data[, rindex := (((group - 1) * (r_tot * f_tot * l_tot)) + ((reg_ind - 1)  # Group, regional, functional, and longitudinal index 
                 + r_tot * (func_ind - 1) + r_tot * f_tot * (long_ind - 1) + 1))]
  
# Create list for output
HPCA = matrix(list(), 6, 1)
for(h in 2:8){
  HPCA[[h]] = matrix(list(), 5, 1)
}
for(g in 1:g_tot){
  HPCA[[2]][[g]] = matrix(list(), r_tot , 1)
}
for(e in 3:6){
  HPCA[[e]] = matrix(list(), g_tot, 1)
}

print("O. Data formatting (completed)")

#############################################################################
# 1. Estimation of Fixed Effects
#############################################################################

# a. Calculate overall mean function

# Create subject-specific longitudinal index 
long_obs<-NULL
for(i in 1:length(id)){
  temp = data[which(data$ID == id[i]), ]
  long = as.matrix(t(unique(temp$long)))
  long = t(as.matrix(rep(long, r_tot))) 
  long_obs = cbind(long_obs, long)
}

# Obtain overall mean function 
mean_surf = fbps(matrix(data$y, nrow = f_tot), covariates = list(x = t(out1), z = long_obs))  # Apply sandwich smoother to all data
newdata = list(x = rep(out1, l_tot), z = rep(out2, each = f_tot))  # Form 2-D grid from longitudinal and functional domains
HPCA[[1]] = matrix(predict(mean_surf, newdata = newdata)$fitted.values, nrow = f_tot)  # Predict overall mean function for 2-D grid

# Center data by overall mean function
mean_vec = data.table(cbind(rep(as.vector(HPCA[[1]]), r_tot),  # Vectorize overall mean function  
           rep(c(1:f_tot), l_tot * r_tot),  # Functional index
           rep(rep(c(1:l_tot), each = f_tot), r_tot),  # Longitudinal index
           rep(c(1:r_tot), each = f_tot*l_tot)))  # Regional index
colnames(mean_vec) = c("mean_fit", "func_ind", "long_ind", "reg_ind")  # Name columns
mean_vec[, index := ((reg_ind - 1) +r_tot * (func_ind - 1) + r_tot * f_tot * (long_ind - 1) + 1)]  # Regional, functional, and longitudinal index
mean_vec = mean_vec[, c(1,5)]  
data = inner_join(data, mean_vec, by="index")  # Merge vectorized overall mean function with data   
data = data.table(data)  # Reformat as data.table following inner_join 
data[, y_c1 := (y - mean_fit)]  # Center data by overall mean function

# b. Calculate group-region-specific shifts

# Obtain group-region-specific shifts
for(g in groups){ 
  for(r in out3){
    temp = data[which((data$reg == r) & (data$group == g)), ]  # Subset data by group-region
    reg_long=NULL
    for(i in unique(temp$ID)){  # Create subject-specific region-long index
      temp_id = temp[which(temp$ID == i), ]
      reg_long = cbind(reg_long, as.matrix(t(unique(temp_id$long))))
    }
    reg_model = fbps(matrix(temp$y_c1, nrow = f_tot), covariates = list(x = t(out1), z = reg_long))  # Apply sandwich smoother to group-region data
    HPCA[[2]][[g]][[r]] = matrix(predict(reg_model, # Predict group-region-specific shifts for 2-D grid
                                 newdata = newdata)$fitted.values, nrow = f_tot)  
  }
}

# Center data by group-region-specific shifts
reg_vec = data.table(cbind(unlist(HPCA[[2]]),  # Vectorize group-region-specific shifts
           rep(rep(c(1:f_tot), l_tot * r_tot),g_tot),  # Functional index
           rep(rep(rep(c(1:l_tot), each = f_tot), r_tot), g_tot),  # Longitudinal index
           rep(rep(c(1:r_tot), each = f_tot * l_tot), g_tot),  # Regional index
           rep(groups, each = r_tot * f_tot * l_tot)))  # Regional, functional, and longitudinal index
colnames(reg_vec) = c("reg_fit", "func_ind", "long_ind", "reg_ind", "g_ind")  # Name columns
reg_vec[, rindex := (((g_ind - 1 ) * (r_tot * f_tot * l_tot))  # Group, regional, functional, and longitudinal index 
                     + ((reg_ind - 1 ) + r_tot * (func_ind - 1) + r_tot * f_tot * (long_ind-1) + 1))] 
reg_vec= reg_vec[, c(1, 6)] 
data = inner_join(data, reg_vec, by="rindex")  # Merge vectorized group-region-specific shifts with data 
data = data.table(data)  # Reformat as data.table following inner_join 
data[, y_c2 := (y_c1 - reg_fit)]  # Center data by group-region-specific shifts

print("1. Estimation of Fixed Effects (completed)")

#############################################################################
# 2. Estimation of Marginal Covariances and Measurement Error Variance
#############################################################################

# Function for obtaining marginal covariances across regional, functional, or longitudinal dimensions
marginal_covar <- function(data,      # data.table from preceding estimation procedures
                           g,         # a group indicator 
                           dom_marg,  # marginal domain 
                           ind_marg,  # marginal domain index 
                           ind_dim1,  # fixed domain index 1 
                           ind_dim2,  # fixed domain index 2 
                           lm,        # number of grid points in marginal domain 
                           ld1,       # number of grid points in fixed domain 1
                           ld2,       # number of grid points in fixed domain 2 
                           smooth     # 1 to turn on covariance smoothing, 0 otherwise
                           ){
  data = data[which(data$group == g), ]    # Subset data by group
  group_subj = as.matrix(unique(data$ID))  # Obtain unique group ID's
  lid = length(group_subj)  # Number of subjects in group
  data[, row_ind := ((match(ID, group_subj) - 1) # Assign an index to the columns of the marginal design matrix
                    + lid * (eval(ind_dim1) - 1) + lid * ld1 * (eval(ind_dim2) - 1) + 1)]
  sparse_mat = sparseMatrix(i = data[, row_ind],  # Form sparse marginal design matrix
                            j = data[, eval(ind_marg)], x = data[, y_c2], dims = c(lid * ld1 * ld2, lm)) 
  sparse_mat_ind = sparseMatrix(i = data[,row_ind], # Form sparse marginal indicator matrix
                                j = data[, eval(ind_marg)], x = 1, dims = c(lid * ld1 * ld2, lm)) 
  cov_mat=as.matrix(crossprod(sparse_mat, sparse_mat) # Calculate covariance matrix 
                    /crossprod(sparse_mat_ind, sparse_mat_ind))
  if(smooth == 1){  # Marginal functional or longitudinal covariance
    cov_vec = data.frame(as.vector(cov_mat),  # Vectorize covariance for smoothing
                         w1 = (rep(dom_marg, lm)), w2 = (rep(dom_marg, each = lm))) 
    colnames(cov_vec) = c("covar", "w1", "w2")
    cov_vec = cov_vec[complete.cases(cov_vec), ]  # Remove entries without any observations
    cov_vec_d = cov_vec[which(cov_vec$w1 != cov_vec$w2), ] # Remove diagonal entries
    cov_vec_s = gam( covar ~ te(w1, w2, k = 10), data = cov_vec_d)  # Smooth the pooled sample covariance
    newdata = data.frame(w1 = (rep(dom_marg, lm)), w2 = (rep(dom_marg, each = lm)))  # Form 2-D grid to predict covariance function
    cov_mat_s = matrix(predict(cov_vec_s, newdata = newdata), nrow = lm)  # Estimated marginal covariance function
    cov_mat_s = (cov_mat_s + t(cov_mat_s)) / 2  #  Symmetrize covariance function
    
    # Calculate measurement error variance
    marg_diag = cov_vec[which( cov_vec$w1 == cov_vec$w2), ]  # Extract diagonal entries from pooled sample covariance
    loess_diag = suppressWarnings(loess.as(marg_diag[, 2], marg_diag[, 1], degree = 1, 
                                           criterion = "gcv", user.span = NULL, plot = F))  # Smooth diagonal entries
    sigma = mean(predict(loess_diag, dom_marg) - diag(cov_mat_s))  # Calculate measurement error variance
    return(list(cov_mat, cov_mat_s, sigma))
  } else {  # Marginal regional covariance
    cov_mat_s = cov_mat-diag(((HPCA[[3]][[g]][[1]][[3]] + HPCA[[4]][[g]][[1]][[3]]) / 2), r_tot)  # Subtract measurement error variance from diagonals
    cov_mat_s = (cov_mat_s + t(cov_mat_s)) / 2  #  Symmetrize covariance function
    return(list(cov_mat, cov_mat_s))
  }
}

# Estimate marginal covariances for each group
for(g in groups){ 
  HPCA[[3]][[g]] = matrix(list(),2,1)
  HPCA[[4]][[g]] = matrix(list(),2,1)
  HPCA[[5]][[g]] = matrix(list(),2,1)
  
  # a. Calculate functional marginal covariance
  HPCA[[3]][[g]][[1]] = marginal_covar(data, g, out1, quote(func_ind), quote(long_ind), quote(reg_ind), f_tot, l_tot, r_tot, 1)
  # b. Calculate longitudinal marginal covariance 
  HPCA[[4]][[g]][[1]] = marginal_covar(data, g, out2, quote(long_ind), quote(func_ind), quote(reg_ind), l_tot, f_tot, r_tot, 1)
  # c. Calculate regional marginal covariance 
  HPCA[[5]][[g]][[1]] = marginal_covar(data, g, out3, quote(reg_ind), quote(func_ind), quote(long_ind), r_tot, f_tot, l_tot, 0)
}

print("2. Estimation of Marginal Covariances and Measurement Error Variance (completed)")

#############################################################################
# 3. Estimation of Marginal Eigencomponents
#############################################################################

# Function for obtaining marginal eigencomponents that explain 85% FVE
eigenfunctions <- function(covariance,  # marginal covariance function
                           dom_marg,    # marginal domain
                           type         # 1 indicates smooth covariance function, 0 indicates covariance matrix
                           ){ 
  eigen_temp = eigen(covariance, symmetric = TRUE)
  eigen_temp$values = eigen_temp$values[which(eigen_temp$values > 0)]  # Obtain positive eigenvalues
  eigen_temp$vectors = eigen_temp$vectors[, 1:length(eigen_temp$values)]  # Obtain eigenvectors associated with positive eigenvalues
  if(type == 1){
    for(e in 1:length(eigen_temp$values)){  # Normalize the eigenvalues over the domain
      eigen_temp$vectors[, e] = eigen_temp$vectors[, e] / sqrt(trapz(dom_marg, eigen_temp$vectors[, e]^2))
    }
  }
  K = length(which(cumsum(eigen_temp$values) / sum(eigen_temp$values) < .85)) + 1  # Select components that explain at least 85% FVE
  return(list("values" = eigen_temp$values, "vectors" = eigen_temp$vectors,"FVE" = K))  
}

# Estimate eigencomponents for each group
for(g in groups){ 
  # a. Employ PCA on regional marginal covariance
  HPCA[[5]][[g]][[2]] = eigenfunctions(HPCA[[5]][[g]][[1]][[2]], out3, 0)
  # b. Employ FPCA on functional marginal covariance
  HPCA[[3]][[g]][[2]] = eigenfunctions(HPCA[[3]][[g]][[1]][[2]], out1, 1) 
  # c. Employ FPCA on longitudinal marginal covariance
  HPCA[[4]][[g]][[2]] = eigenfunctions(HPCA[[4]][[g]][[1]][[2]], out2, 1)
}

print("3. Estimation of Marginal Eigencomponents (completed)")

#############################################################################
# 4. Estimation of Variance Components and Subject-specific Scores
#############################################################################

# Function for indexing the multidimensional orthonormal basis
eig_index <- function(x){
  ind = (x[1] - 1)+ HPCA[[5]][[g]][[2]]$FVE * (x[2] - 1) + HPCA[[5]][[g]][[2]]$FVE * HPCA[[3]][[g]][[2]]$FVE * (x[3] - 1) + 1;
  return(ind)
}

# Form vectorized versions of the multidimensional orthonormal basis 
prod_surf = matrix(list(), nrow=g_tot, 1)
for(g in groups){
  prod_surf[[g]] = data.table(matrix(NA,  # Define dimensions of group=specific matrix
                                     nrow = r_tot * f_tot * l_tot,
                                     ncol = HPCA[[3]][[g]][[2]]$FVE * HPCA[[4]][[g]][[2]]$FVE * HPCA[[5]][[g]][[2]]$FVE + 4))
  prod_surf[[g]][, 1] = as.vector(rep(kronecker(c(1:r_tot), matrix(1, nrow = l_tot, ncol = 1)), f_tot))  # Regional index
  prod_surf[[g]][, 2] = kronecker(c(1:f_tot), matrix(1, nrow = l_tot * r_tot))  # Functional index
  prod_surf[[g]][, 3] = rep(c(1:l_tot), r_tot * f_tot)  # Longitudinal index
  prod_surf[[g]][, 4 := ((V1 - 1) + r_tot * (V2 - 1) + r_tot * f_tot * (V3 - 1) + 1)]  # Regional, functional, and longitudinal index 
  for(k in 1:HPCA[[5]][[g]][[2]]$FVE){
    for(l in 1:HPCA[[3]][[g]][[2]]$FVE){
      for(m in 1:HPCA[[4]][[g]][[2]]$FVE){
        v_k = HPCA[[5]][[g]][[2]]$vectors[, k]  # kth region marginal eigenvector
        phi_l = HPCA[[3]][[g]][[2]]$vectors[, l]  # lth functional marginal eigenfunction
        psi_m = HPCA[[4]][[g]][[2]]$vectors[, m]  # mth longitudinal marginal eigenfunction
        prod_surf[[g]][, (4 + eig_index(c(k, l, m, g))) := (v_k[V1] * phi_l[V2] * psi_m[V3])]  # Multdimensional orthonormal basis
      }
    }
  }
  prod_surf[[g]] = prod_surf[[g]][, -c(1:3)]
  colnames(prod_surf[[g]]) = c("index", paste("x", 1:(dim(prod_surf[[g]])[2] - 1), sep = ""))
}

# Create group-specific design matrix 
data_lmm = matrix(list(), nrow=g_tot, 1)
for(g in groups){
  data_lmm[[g]] = data[which(data$group == g), ]
  data_lmm[[g]] = inner_join(data_lmm[[g]], prod_surf[[g]], by = "index")
}

# a-b. Calculate variance components, subject-specific scores, and measurement 
#      error variance by linear mixed effects model for each group
for(g in groups){
  HPCA[[6]][[g]] = matrix(list(), 4)
  xnam = paste("x", 1:(dim(prod_surf[[g]])[2] - 1), sep = "")
  formula = as.formula(paste('~0+', paste(xnam, collapse = "+"), '|ID'))
  RE_struct = reStruct(formula, pdClass = 'pdDiag', data = data_lmm[[g]])  # Specify linear mixed effects model
  RE_model = lme(y_c2 ~ -1, random = RE_struct, data = data_lmm[[g]], control = lmeControl(returnObject = TRUE, opt="optim"))  
  HPCA[[6]][[g]][[1]] = random.effects(RE_model)  # Extract subject-specific scores
  HPCA[[6]][[g]][[2]] = getVarCov(RE_model)  # Extract variance components
  HPCA[[6]][[g]][[3]] = RE_model$sigma^2  # Extract measurement error variance
  HPCA[[6]][[g]][[4]] = c("subj_score", "tau_klm", "sigma")
}

print("4. Estimation of Variance Components and Subject-specific Scores (completed)")

#############################################################################

HPCA[[7]] = do.call("rbind", data_lmm)

HPCA[[8]] = c("mu","eta","func marg covar","long marg covar","reg marg covar","RE model","data")

return(HPCA)

}
