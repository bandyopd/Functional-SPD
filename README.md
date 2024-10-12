R codes for simulation studies with small n, d and p can be found here. Detailed descriptions of codes can be found in the codes. One can easily modify the codes for larger n, d and p. The fortran dll files were created from some fortran codes that implement the SBF methods. Unfortunately, the fortran codes will not be distributed until an R package is made.


Simulation 1:
--

To get the prediction errors of proposed method, please follow the steps:
--

Step 1: Download SBF_vector_valued_L2_new_grid.dll, CBS_vector_valued_L2_new_grid.dll, Code_for_simulation1.R and Code_for_simulation1_proposed.R

Step 2: Set the addresses of SBF_vector_valued_L2_new_grid.dll and CBS_vector_valued_L2_new_grid.dll in Code_for_simulation1.R 

Step 3: Set the address of Code_for_simulation1.R in Code_for_simulation1_proposed.R

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

Simulation 2
--

To perform the proposed method, please follow the steps:

Step 1: Download SBF_vector_valued_Euclidean_new_grid_PL.dll, SBF_vector_valued_Euclidean_new_grid_PL_tilde_out.dll, CBS_vector_valued_Euclidean_new_grid_PL.dll, Code_for_simulation2.R, Code_for_simulation2_proposed.R and Code_for_simulation2_proposed_testing.R

Step 2: Set the addresses of SBF_vector_valued_Euclidean_new_grid_PL.dll, SBF_vector_valued_Euclidean_new_grid_PL_tilde_out.dll and CBS_vector_valued_Euclidean_new_grid_PL.dll in Code_for_simulation2.R

Step 3: Set the address of Code_for_simulation2.R in Code_for_simulation2_proposed.R and Code_for_simulation2_proposed_testing.R

Step 4: Run Code_for_simulation2_proposed.R to get the prediction result for (n,d,p)=(50,2,1)

Step 5: Run Code_for_simulation2_proposed_testing.R to get the testing result for (n,d,p)=(50,2,1)

Step 6: Change n in Code_for_simulation2_proposed.R and Code_for_simulation2_proposed_testing.R to get the results for different n

To get the prediction errors of global frechet regression, please follow the steps:
--

Step 1: Download Code_for_simulation2_global_frechet.R

Step 2: Run Code_for_simulation2_global_frechet.R to get the result for (n,d,p)=(50,2,1)

Step 3: Change n in Code_for_simulation2_global_frechet.R to get the results for different n
