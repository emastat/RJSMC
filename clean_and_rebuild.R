#!/usr/bin/env Rscript
# Script to clean all compiled cache and rebuild the package from scratch

cat("Cleaning compiled files...\n")

# Remove object files
system("find src -name '*.o' -delete", ignore.stdout = TRUE, ignore.stderr = TRUE)
system("find src -name '*.so' -delete", ignore.stdout = TRUE, ignore.stderr = TRUE)
system("find src -name '*.dylib' -delete", ignore.stdout = TRUE, ignore.stderr = TRUE)

# Remove RcppExports
if (file.exists("src/RcppExports.cpp")) {
  file.remove("src/RcppExports.cpp")
  cat("Removed src/RcppExports.cpp\n")
}
if (file.exists("src/RcppExports.R")) {
  file.remove("src/RcppExports.R")
  cat("Removed src/RcppExports.R\n")
}

# Clean DLLs using devtools
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::clean_dll()
  cat("Cleaned DLLs\n")
}

# Clean R cache
if (dir.exists(".Rproj.user")) {
  unlink(".Rproj.user", recursive = TRUE)
  cat("Removed .Rproj.user cache\n")
}

cat("\nNow rebuilding package...\n")
cat("Run: devtools::load_all() or devtools::install()\n")

