CONTENTS OF THIS FOLDER
——————————————	

* HPCA_tutorial.R : 			 A step-by-step implementation of HPCA_decomposition.R and the associated  
   					procedures described in "Hybrid Principal Components Analysis For Region-Referenced Longitudinal 
 					Functional EEG Data" by Scheffler et al. (2018).

* HPCA_decomposition.R : 	Function for performing HPCA decomposition including estimation of fixed effects, marginal covariance functions,
            			           marginal eigencomponents, subject-specific scores, variance components, and measurement 
			                      error variance.

* HPCA_simulation.R :		Function for simulating data for the HPCA_decomp function and the associated procedures referred to above.

* HPCA_bootstrap.R : 		Function for performing group-level inference for a given scalp region. The null hypothesis assumes that groups share 
					a common region shift for a given region.
           				
		 
INTRODUCTION
——————————————	

The contents of this folder allow for implementation of the HPCA decomposition described in  "Hybrid Principal Components Analysis
For Region-Referenced Longitudinal Functional EEG Data" by Scheffler et al. (2018). Users can simulate a sample data frame (HPCA_simulation.R) and 
apply the proposed HPCA decomposition (HPCA_decomposition.R). Further, we include tools to perform group-level inference via a bootstrap procedure
(HCPA_bootstrap.R), allowing users to test whether the longitudinal functional stochastic process in a fixed region varies among groups. Detailed instructions
on how to perform the aforementioned procedures, visualize results, and check the assumption of weak separability via a likelihood-ratio test on the 
correlation structure of the random effects are included in HPCA_tutorial.R. 

REQUIREMENTS
——————————————	

The included R programs require R 3.3.2 (R Core Team, 2016) and the packages listed in HPCA_tutorial.R. 


INSTALLATION
——————————————	

Load the R program files into the global environment and install required packages using commands in HPCA_tutorial.R.
