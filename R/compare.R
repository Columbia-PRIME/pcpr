##################################################
# Compare R output with original MATLAB functions#
# (Sanity Check) #################################
##################################################

# library(tidyverse)
# library(R.matlab)
# options(scipen = 999)

# Read in Boston air pollution data
# mixture <- readMat("./Data/mixtures_data.mat")
#
# mix_data <- as.data.frame(mixture) %>% as_tibble() %>%
#   select(Al, As, Ba, bc, Br, Ca, Cl,
#          Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
#          Ti,  V, Zn) %>%
#   drop_na()
#
# summary(mix_data)
#
# # Create version with 50% lowest values for each variable as below the LOD
# mix_data_50 <- mix_data %>%
#   mutate_all(~ifelse(. <= quantile(., probs = .50), -1, .))
#
# # For MATLAB
# # write_csv(mix_data_50, "./Data/mix_data_50.csv")
#
# # LOD quantiles at scalar, vector, and matrix
#
# # vector LOD
# vector50 <- mix_data %>% summarise_all(quantile, probs = .50) %>% as_vector()
#
# # matrix LOD
# empty = matrix(1, nrow = nrow(mix_data), ncol = ncol(mix_data))
# matrix50 = t(t(empty) * vector50)
#
# # For MATLAB
# # write_csv(as_tibble(matrix50), "./Data/matrix50.csv")
#
# # Run PCP-LOD on 3 LODs
# lambda <- 1/sqrt(nrow(mix_data_50)) # default lambda
# mix_data_50 <- as.matrix(mix_data_50)
#
# # scalar LOD
# out_s <- pcp_lod(mix_data_50, lambda, 10, 0.01)
# low_s <- out_s[[1]]
# sparse_s <- out_s[[2]]
#
# # vector LOD
# out_v <- pcp_lod(mix_data_50, lambda, 10, vector50)
# low_v <- out_v[[1]]
# sparse_v <- out_v[[2]]
#
# # matrix LOD
# out_m <- pcp_lod(mix_data_50, lambda, 10, matrix50)
# low_m <- out_m[[1]]
# sparse_m <- out_m[[2]]
#
# # Load matlab results
#
# mat_low_s <- readMat("./Data/lowrank_50s.mat")[[1]]
# mat_low_v <- readMat("./Data/lowrank_50v.mat")[[1]]
# mat_low_m <- readMat("./Data/lowrank_50m.mat")[[1]]
#
# mat_sparse_s <- readMat("./Data/sparse_50s.mat")[[1]]
# mat_sparse_v <- readMat("./Data/sparse_50v.mat")[[1]]
# mat_sparse_m <- readMat("./Data/sparse_50m.mat")[[1]]
#
# # Compare
#
# # scalar LOD
# norm(low_s    - mat_low_s,    type = "F")
# norm(sparse_s - mat_sparse_s, type = "F")
# # Difference is basically zero, good!
#
# # vector LOD
# norm(low_v    - mat_low_v,    type = "F")
# norm(sparse_v - mat_sparse_v, type = "F")
# # Difference is basically zero, good!
#
# # matrix LOD
# norm(low_m    - mat_low_m,    type = "F")
# norm(sparse_m - mat_sparse_m, type = "F")
# # Difference is basically zero, good!
#
