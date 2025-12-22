library(RJSMC)

# Load your real data (adjust paths as needed)
abb_data <- readRDS("/Users/lele_phd/Desktop/di21_final.RDS")
Tvec <- abb_data$time
Yvec <- as.numeric(abb_data$label)

# Load all parameters
lambdamat <- readRDS("/Users/lele_phd/Desktop/parameter_values/lambdamat_par.RDS")
keyvec <- readRDS("/Users/lele_phd/Desktop/parameter_values/greenvec_par.RDS")
etavec <- readRDS("/Users/lele_phd/Desktop/parameter_values/etavec_par.RDS")
key0vec <- readRDS("/Users/lele_phd/Desktop/parameter_values/greenvec_empty_par.RDS")
eta0vec <- readRDS("/Users/lele_phd/Desktop/parameter_values/etavec_empty_par.RDS")
alphavec <- readRDS("/Users/lele_phd/Desktop/parameter_values/alphavec_par.RDS")
muvec <- readRDS("/Users/lele_phd/Desktop/parameter_values/muvec_par.RDS")
probvec_Z <- readRDS("/Users/lele_phd/Desktop/parameter_values/qvec_par.RDS")
probvec_V <- readRDS("/Users/lele_phd/Desktop/parameter_values/deltavec_par.RDS")
probvec_F <- readRDS("/Users/lele_phd/Desktop/parameter_values/epsivec_par.RDS")
probvec_Q <- c(0.5, 0.5)
P0 <- readRDS("/Users/lele_phd/Desktop/parameter_values/P0_par.RDS")
minimum_n <- 0

# Calculate derived values
K <- length(alphavec)
U <- length(probvec_V)
W <- 2
num_logs <- ncol(lambdamat)

# Create the dataset list (similar to english_words structure)
# Note: For real data, we don't have true values (Bvec, Vvec, etc.)
# so we set them to NULL or omit them
real_data <- list(
  # Time series data
  Tvec = Tvec,
  Yvec = Yvec,
  
  # Parameters
  U = U,
  W = W,
  K = K,
  num_logs = num_logs,
  minimum_n = minimum_n,
  
  # Probability vectors
  probvec_V = probvec_V,
  probvec_Z = probvec_Z,
  probvec_Q = probvec_Q,
  probvec_F = probvec_F,
  P0 = P0,
  
  # Parameter matrices and vectors
  lambdamat = lambdamat,
  alphavec = alphavec,
  muvec = muvec,
  keyvec = keyvec,
  etavec = etavec,
  key0vec = key0vec,
  eta0vec = eta0vec,
  
  # True values (set to NULL for real data since we don't know them)
  Bvec = NULL,
  Vvec = NULL,
  Zvec = NULL,
  Qvec = NULL,
  Fvec = NULL,
  Mvec = NULL,
  Nvec = NULL
)

# Save the dataset

usethis::use_data(real_data, overwrite = TRUE)

# Option 2: Using save() directly (alternative)
#save(real_data, file = "data/real_data.rda", compress = "xz")

# # Verify the structure
# cat("Dataset created with the following elements:\n")
# cat(paste(names(real_data), collapse = ", "), "\n")
# cat("\nDataset saved to: data/real_data.rda\n")
