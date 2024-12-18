## Description

R code to conduct simulation studies corresponding to small (n,d,p) from Jeong and Bandyopadhyay (2024+) are made available. One can easily modify the code for larger (n,d,p). The fortran .dll files were created from some fortran codes that implement the SBF methods. The .dll files will later integrated within a forthcomong R package. 


Simulation 1:
--

To get the prediction errors of proposed method, please follow the steps:
--

Step 1: Download 

- SBF_vector_valued_L2_new_grid.dll
- CBS_vector_valued_L2_new_grid.dll
- Code_for_simulation1.R
- Code_for_simulation1_proposed.R

Step 2: Open Code_for_simulation1.R and change the addresses of SBF_vector_valued_L2_new_grid.dll and CBS_vector_valued_L2_new_grid.dll there

Step 3: Open Code_for_simulation1_proposed.R and change the address of Code_for_simulation1.R there

Step 4: Run Code_for_simulation1_proposed.R to get the result for (n,d)=(50,2)

Step 5: Change n in Code_for_simulation1_proposed.R to get the results for different n

To get the prediction errors of global frechet regression, please follow the steps:
--

Step 1: Download Code_for_simulation1_global_frechet.R

Step 2: Run Code_for_simulation1_global_frechet.R to get the result for (n,d)=(50,2)

Step 3: Change n in Code_for_simulation1_global_frechet.R to get the results for different n

-------------
To get the prediction errors of local constant frechet regression, please follow the steps:
--

Step 1: Download Code_for_simulation1_local_frechet.R

Step 2: Run Code_for_simulation1_local_frechet.R to get the result for (n,d)=(50,2)

Step 3: Change n in Code_for_simulation1_local_frechet.R to get the results for different n

Simulation 2:
--

To perform the proposed method, please follow the steps:
--

Step 1: Download 

- SBF_vector_valued_Euclidean_new_grid_PL.dll
- SBF_vector_valued_Euclidean_new_grid_PL_tilde_out.dll
- CBS_vector_valued_Euclidean_new_grid_PL.dll
- Code_for_simulation2.R
- Code_for_simulation2_proposed.R
- Code_for_simulation2_proposed_testing.R

Step 2: Open Code_for_simulation2.R and change the addresses of SBF_vector_valued_Euclidean_new_grid_PL.dll, SBF_vector_valued_Euclidean_new_grid_PL_tilde_out.dll and CBS_vector_valued_Euclidean_new_grid_PL.dll there

Step 3: Open Code_for_simulation2_proposed.R and Code_for_simulation2_proposed_testing.R, and change the addresses of Code_for_simulation2.R there

Step 4: Run Code_for_simulation2_proposed.R to get the prediction result for (n,d,p)=(50,2,1)

Step 5: Run Code_for_simulation2_proposed_testing.R to get the testing result for (n,d,p)=(50,2,1)

Step 6: Change n in Code_for_simulation2_proposed.R and Code_for_simulation2_proposed_testing.R to get the results for different n

To get the prediction errors of global frechet regression, please follow the steps:
--

Step 1: Download Code_for_simulation2_global_frechet.R

Step 2: Run Code_for_simulation2_global_frechet.R to get the result for (n,d,p)=(50,2,1)

Step 3: Change n in Code_for_simulation2_global_frechet.R to get the results for different n

## Reference: 

Jeon JM and Bandyopadhyay D. (2024+). Structured nonparametric regression for functional DTI data, (Under Review). 
