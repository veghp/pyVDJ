# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("powerTCR")
library(powerTCR)
# Load data (where we have at least 100 cells):
path_counts <- "count_vectors"
vdj_counts <- list()
for(f in dir(path_counts)) {
  f_values <- scan(file.path(path_counts, f))
  if (length(f_values) > 100) {
    sample <- gsub("count_vector_", "", f)
    sample <- gsub(".csv", "", sample)
    vdj_counts[[sample]] <- f_values
  }
}

fits <- list()
for(i in seq_along(vdj_counts)) {
  print(i)
  # Choose a sequence of possible u for the model fit
  thresholds <- unique(round(quantile(vdj_counts[[i]], c(0.7, 0.73, .75, 0.77, .8, 0.83, .85, 0.88, .9))))

  fits[[i]] <- fdiscgammagpd(
    vdj_counts[[i]],
    useq = thresholds,
    shift = min(vdj_counts[[i]]))
}
names(fits) <- names(vdj_counts)

thresholds <- c()
for(i in seq_along(fits)) {
  thresholds[i] <- get_mle(fits)[[i]]["thresh"]
}

thresholds_df <- data.frame(name = names(fits))
thresholds_df$threshold <- thresholds
thresholds_df
write.csv(thresholds_df, file = "threshold_df.csv", quote = F, row.names = F)
