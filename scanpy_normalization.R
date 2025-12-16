# R function to replicate scanpy's sc.pp.normalize_total() and sc.pp.log1p() 
# transformation pipeline commonly used in single-cell analysis

#' Scanpy-style normalization and log transformation
#' 
#' This function replicates the scanpy preprocessing pipeline:
#' 1. sc.pp.normalize_total(x) - normalizes counts to a target sum per cell
#' 2. sc.pp.log1p(x) - applies log(x + 1) transformation
#' 
#' @param count_matrix A numeric matrix or data.frame with features/genes as rows and cells as columns (scanpy convention)
#' @param target_sum Target sum for normalization. If NULL (default), uses median of total counts per cell (scanpy default behavior). If numeric, uses that value (e.g., 1e4)
#' @param copy Whether to return a copy or modify in place (default: TRUE)
#' @return Normalized and log-transformed matrix
#' 
#' @examples
#' # Example usage:
#' # raw_counts <- matrix(rpois(1000, 5), nrow=10, ncol=100)  # 10 genes, 100 cells
#' # normalized <- scanpy_normalize_log1p(raw_counts)
scanpy_normalize_log1p <- function(count_matrix, target_sum = NULL, copy = TRUE) {
  
  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }
  
  if (!is.null(target_sum) && (!is.numeric(target_sum) || target_sum <= 0)) {
    stop("target_sum must be NULL or a positive number")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Make copy if requested
  if (copy) {
    result <- count_matrix
  } else {
    result <- count_matrix
  }
  
  # Step 1: Normalize total counts per cell (sc.pp.normalize_total)
  # Calculate sum per column (cell) - features are on rows, cells on columns
  col_sums <- colSums(result)
  
  # Handle zero-sum columns to avoid division by zero
  zero_cols <- col_sums == 0
  if (any(zero_cols)) {
    warning(sprintf("Found %d cells with zero total counts. These will remain unchanged.", 
                   sum(zero_cols)))
  }
  
  # Determine target sum: if NULL, use median of column sums (scanpy default behavior)
  if (is.null(target_sum)) {
    # Only consider non-zero cells for median calculation
    non_zero_sums <- col_sums[col_sums > 0]
    if (length(non_zero_sums) == 0) {
      stop("All cells have zero counts. Cannot determine target sum.")
    }
    target_sum <- median(non_zero_sums)
    message(sprintf("Using median total counts as target_sum: %.2f", target_sum))
  }
  
  # Normalize each cell to target_sum
  # This replicates: each cell's counts / cell_total * target_sum
  for (j in 1:ncol(result)) {
    if (col_sums[j] > 0) {
      result[, j] <- result[, j] * (target_sum / col_sums[j])
    }
  }
  
  # Step 2: Apply log1p transformation (sc.pp.log1p)
  # This is log(x + 1) transformation
  result <- log1p(result)
  
  return(result)
}

#' Alternative vectorized implementation for better performance on large matrices
#' 
#' @param count_matrix A numeric matrix with features/genes as rows and cells as columns (scanpy convention)
#' @param target_sum Target sum for normalization. If NULL (default), uses median of total counts per cell (scanpy default). If numeric, uses that value
#' @return Normalized and log-transformed matrix
scanpy_normalize_log1p_fast <- function(count_matrix, target_sum = NULL) {
  
  # Convert to matrix if needed
  if (!is.matrix(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Calculate column sums (cells are columns)
  col_sums <- colSums(count_matrix)
  
  # Determine target sum: if NULL, use median of column sums (scanpy default behavior)
  if (is.null(target_sum)) {
    # Only consider non-zero cells for median calculation
    non_zero_sums <- col_sums[n_zero_sums) == 0) {
      stop("All cells have zero counts. Cannot determine target sum.")
    }
    target_sum <- median(non_zero_sums)
    message(sprintf("Using median total counts as target_sum: %.2f", target_sum))
  }
  
  # Create scaling factors
  scaling_factors <- target_sum / col_sums
  scaling_factors[col_sums == 0] <- 0  # Handle zero-sum columns
  
  # Vectorized normalization using sweep (sweep over columns)
  normalized <- sweep(count_matrix, 2, scaling_factors, "*")
  
  # Apply log1p transformation
  result <- log1p(normalized)
  
  return(result)
}

#' Check if the transformation matches expected scanpy behavior
#' 
#' @param count_matrix Input count matrix
#' @param target_sum Target sum for normalization. If NULL, uses median of total counts per cell
#' @return List with diagnostic information
check_scanpy_transformation <- function(count_matrix, target_sum = NULL) {
  
  result <- scanpy_normalize_log1p_fast(count_matrix, target_sum)
  
  # Back-transform to check normalization worked correctly
  # exp(result) - 1 should have column sums approximately equal to target_sum
  back_transformed <- expm1(result)  # expm1(x) = exp(x) - 1
  actual_sums <- colSums(back_transformed)
  
  # Check for cells that should sum to target_sum (non-zero original cells)
  original_sums <- colSums(count_matrix)
  non_zero_cells <- original_sums > 0
  
  # Determine what target_sum was used
  if (is.null(target_sum)) {
    target_sum <- median(original_sums[non_zero_cells])
  }
  
  list(
    original_sums = original_sums,
    target_sum = target_sum,
    actual_sums_after_backtransform = actual_sums,
    mean_deviation_from_target = mean(abs(actual_sums[non_zero_cells] - target_sum)),
    max_deviation_from_target = max(abs(actual_sums[non_zero_cells] - target_sum)),
    all_within_tolerance = all(abs(actual_sums[non_zero_cells] - target_sum) < 1e-10)
  )
}

#' cTPnet-style log normalization using geometric mean
#' 
#' This function replicates the cTPnet normalization:
#' np.log((x+1.0)/stats.gmean(x+1.0), axis=0)
#' 
#' @param count_matrix A numeric matrix or data.frame with cells as rows and genes/features as columns
#' @param axis Which axis to compute geometric mean along. 0 = columns (genes), 1 = rows (cells)
#' @param copy Whether to return a copy or modify in place (default: TRUE)
#' @return Log-normalized matrix using geometric mean scaling
#' 
#' @examples
#' # Example usage:
#' # raw_counts <- matrix(rpois(1000, 5), nrow=100, ncol=10)
#' # normalized <- ctpnet_log_normalize(raw_counts)
ctpnet_log_normalize <- function(count_matrix, axis = 0, copy = TRUE) {
  
  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }
  
  if (!axis %in% c(0, 1)) {
    stop("axis must be 0 (columns/genes) or 1 (rows/cells)")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Make copy if requested
  if (copy) {
    result <- count_matrix
  } else {
    result <- count_matrix
  }
  
  # Add 1 to all values (x + 1.0)
  x_plus_1 <- result + 1.0
  
  if (axis == 0) {
    # Compute geometric mean along columns (for each gene across all cells)
    # This matches numpy's axis=0 behavior
    
    # Calculate geometric mean for each column
    geom_means <- apply(x_plus_1, 2, function(col) {
      # Handle zeros by using only positive values for geometric mean
      positive_vals <- col[col > 0]
      if (length(positive_vals) == 0) {
        return(1)  # Avoid division by zero
      }
      exp(mean(log(positive_vals)))
    })
    
    # Divide each column by its geometric mean
    normalized <- sweep(x_plus_1, 2, geom_means, "/")
    
  } else if (axis == 1) {
    # Compute geometric mean along rows (for each cell across all genes)
    # This matches numpy's axis=1 behavior
    
    # Calculate geometric mean for each row
    geom_means <- apply(x_plus_1, 1, function(row) {
      # Handle zeros by using only positive values for geometric mean
      positive_vals <- row[row > 0]
      if (length(positive_vals) == 0) {
        return(1)  # Avoid division by zero
      }
      exp(mean(log(positive_vals)))
    })
    
    # Divide each row by its geometric mean
    normalized <- sweep(x_plus_1, 1, geom_means, "/")
  }
  
  # Apply log transformation
  result <- log(normalized)
  
  return(result)
}

#' Alternative implementation using built-in geometric mean function
#' 
#' @param count_matrix A numeric matrix with cells as rows and genes/features as columns
#' @param axis Which axis to compute geometric mean along. 0 = columns (genes), 1 = rows (cells)
#' @return Log-normalized matrix using geometric mean scaling
ctpnet_log_normalize_fast <- function(count_matrix, axis = 0) {
  
  # Convert to matrix if needed
  if (!is.matrix(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Add 1 to all values
  x_plus_1 <- count_matrix + 1.0
  
  if (axis == 0) {
    # Geometric mean for each column (gene)
    geom_means <- apply(x_plus_1, 2, function(x) exp(mean(log(x[x > 0]))))
    geom_means[is.na(geom_means) | geom_means == 0] <- 1  # Handle edge cases
    
    # Normalize and log
    result <- log(sweep(x_plus_1, 2, geom_means, "/"))
    
  } else {
    # Geometric mean for each row (cell)
    geom_means <- apply(x_plus_1, 1, function(x) exp(mean(log(x[x > 0]))))
    geom_means[is.na(geom_means) | geom_means == 0] <- 1  # Handle edge cases
    
    # Normalize and log
    result <- log(sweep(x_plus_1, 1, geom_means, "/"))
  }
  
  return(result)
}

#' Check cTPnet transformation behavior
#' 
#' @param count_matrix Input count matrix
#' @param axis Which axis to compute geometric mean along
#' @return List with diagnostic information
check_ctpnet_transformation <- function(count_matrix, axis = 0) {
  
  result <- ctpnet_log_normalize_fast(count_matrix, axis)
  
  # Back-transform to check normalization
  back_transformed <- exp(result) - 1
  
  # Calculate geometric means of original data + 1
  x_plus_1 <- count_matrix + 1.0
  
  if (axis == 0) {
    original_geom_means <- apply(x_plus_1, 2, function(x) exp(mean(log(x[x > 0]))))
    # After normalization, geometric mean should be 1 for each column
    normalized_geom_means <- apply(exp(result), 2, function(x) exp(mean(log(x[x > 0]))))
  } else {
    original_geom_means <- apply(x_plus_1, 1, function(x) exp(mean(log(x[x > 0]))))
    # After normalization, geometric mean should be 1 for each row
    normalized_geom_means <- apply(exp(result), 1, function(x) exp(mean(log(x[x > 0]))))
  }
  
  list(
    original_geom_means = original_geom_means,
    normalized_geom_means = normalized_geom_means,
    geom_means_close_to_1 = all(abs(normalized_geom_means - 1) < 1e-10, na.rm = TRUE),
    data_range_after_transform = range(result, na.rm = TRUE),
    any_infinite_values = any(is.infinite(result)),
    any_na_values = any(is.na(result))
  )
}

# Example usage and testing
if (FALSE) {  # Set to TRUE to run examples
  
  # Create example count matrix (genes as rows, cells as columns)
  set.seed(42)
  example_counts <- matrix(rpois(1000, lambda = 5), nrow = 10, ncol = 100)
  rownames(example_counts) <- paste0("Gene_", 1:10)
  colnames(example_counts) <- paste0("Cell_", 1:100)
  
  # Apply scanpy transformation
  normalized_scanpy <- scanpy_normalize_log1p(example_counts)
  
  # Apply cTPnet transformation (axis=0 for genes, axis=1 for cells)
  normalized_ctpnet_genes <- ctpnet_log_normalize(example_counts, axis = 0)
  normalized_ctpnet_cells <- ctpnet_log_normalize(example_counts, axis = 1)
  
  # Check the transformations
  diagnostics_scanpy <- check_scanpy_transformation(example_counts)
  diagnostics_ctpnet_genes <- check_ctpnet_transformation(example_counts, axis = 0)
  diagnostics_ctpnet_cells <- check_ctpnet_transformation(example_counts, axis = 1)
  
  print("Scanpy transformation diagnostics:")
  print(diagnostics_scanpy)
  
  print("cTPnet transformation (axis=0, genes) diagnostics:")
  print(diagnostics_ctpnet_genes)
  
  print("cTPnet transformation (axis=1, cells) diagnostics:")
  print(diagnostics_ctpnet_cells)
  
  # Compare with fast version
  normalized_fast <- scanpy_normalize_log1p_fast(example_counts)
  
  # Should be identical (within machine precision)
  print(all.equal(normalized_scanpy, normalized_fast))
  
  # Show before/after statistics
  cat("Original data:\n")
  cat("  Column sums range:", range(colSums(example_counts)), "\n")
  cat("  Mean counts per cell:", mean(colSums(example_counts)), "\n")
  
  cat("\nAfter scanpy normalization + log1p:\n")
  cat("  Data range:", range(normalized_scanpy), "\n")
  cat("  Mean value:", mean(normalized_scanpy), "\n")
  
  cat("\nAfter cTPnet normalization (genes):\n")
  cat("  Data range:", range(normalized_ctpnet_genes), "\n")
  cat("  Mean value:", mean(normalized_ctpnet_genes), "\n")
  
  cat("\nAfter cTPnet normalization (cells):\n")
  cat("  Data range:", range(normalized_ctpnet_cells), "\n")
  cat("  Mean value:", mean(normalized_ctpnet_cells), "\n")
}

#' Z-score normalization by batch (replicates scanpy's sc.pp.scale by batch)
#' 
#' This function performs z-score normalization (standardization) of features by batch,
#' similar to sciPENN's preprocessing step. For each batch, features are standardized
#' to have mean=0 and standard deviation=1.
#' 
#' @param count_matrix A numeric matrix with features/genes as rows and cells as columns
#' @param batch_vector A vector indicating batch membership for each cell (column)
#' @param copy Whether to return a copy or modify in place (default: TRUE)
#' @return Z-score normalized matrix
#' 
#' @examples
#' # Example usage:
#' # matrix with 100 genes, 50 cells
#' # raw_counts <- matrix(rnorm(5000, mean=5, sd=2), nrow=100, ncol=50)
#' # batch_info <- rep(c("batch1", "batch2"), each=25)
#' # normalized <- zscore_normalize_by_batch(raw_counts, batch_info)
zscore_normalize_by_batch <- function(count_matrix, batch_vector, copy = TRUE) {
  
  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }
  
  if (length(batch_vector) != ncol(count_matrix)) {
    stop("batch_vector length must match number of columns in count_matrix")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Make copy if requested
  if (copy) {
    result <- count_matrix
  } else {
    result <- count_matrix
  }
  
  # Get unique batches
  unique_batches <- unique(batch_vector)
  
  # Z-score normalize each batch separately
  for (batch in unique_batches) {
    # Get indices for this batch
    batch_indices <- which(batch_vector == batch)
    
    if (length(batch_indices) > 1) {
      # Extract data for this batch
      batch_data <- result[, batch_indices, drop = FALSE]
      
      # Calculate mean and standard deviation for each feature (row) in this batch
      row_means <- rowMeans(batch_data)
      row_sds <- apply(batch_data, 1, sd)
      
      # Handle features with zero variance (set to 0 after standardization)
      zero_var_features <- row_sds == 0
      if (any(zero_var_features)) {
        warning(sprintf("Found %d features with zero variance in batch '%s'. These will be set to 0.", 
                       sum(zero_var_features), batch))
        row_sds[zero_var_features] <- 1  # Avoid division by zero
      }
      
      # Z-score normalize: (x - mean) / sd
      for (i in 1:nrow(batch_data)) {
        result[i, batch_indices] <- (batch_data[i, ] - row_means[i]) / row_sds[i]
      }
      
      # Set zero-variance features to 0
      if (any(zero_var_features)) {
        result[zero_var_features, batch_indices] <- 0
      }
      
    } else {
      # Single cell in batch - set to 0 (can't calculate meaningful z-score)
      warning(sprintf("Batch '%s' has only one cell. Setting values to 0.", batch))
      result[, batch_indices] <- 0
    }
  }
  
  return(result)
}

#' Vectorized z-score normalization by batch for better performance
#' 
#' @param count_matrix A numeric matrix with features/genes as rows and cells as columns
#' @param batch_vector A vector indicating batch membership for each cell (column)
#' @return Z-score normalized matrix
zscore_normalize_by_batch_fast <- function(count_matrix, batch_vector) {
  
  # Convert to matrix if needed
  if (!is.matrix(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Initialize result matrix
  result <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
  rownames(result) <- rownames(count_matrix)
  colnames(result) <- colnames(count_matrix)
  
  # Get unique batches
  unique_batches <- unique(batch_vector)
  
  # Process each batch
  for (batch in unique_batches) {
    batch_indices <- batch_vector == batch
    
    if (sum(batch_indices) > 1) {
      # Extract batch data
      batch_data <- count_matrix[, batch_indices, drop = FALSE]
      
      # Calculate means and standard deviations for each feature
      row_means <- rowMeans(batch_data)
      row_sds <- apply(batch_data, 1, sd)
      
      # Handle zero variance features
      zero_var <- row_sds == 0
      row_sds[zero_var] <- 1  # Avoid division by zero
      
      # Vectorized z-score calculation
      batch_normalized <- sweep(batch_data, 1, row_means, "-")  # Subtract means
      batch_normalized <- sweep(batch_normalized, 1, row_sds, "/")  # Divide by sd
      
      # Set zero-variance features to 0
      batch_normalized[zero_var, ] <- 0
      
      # Store results
      result[, batch_indices] <- batch_normalized
      
    } else {
      # Single cell batch - set to 0
      result[, batch_indices] <- 0
    }
  }
  
  return(result)
}

#' Check z-score normalization results
#' 
#' @param normalized_matrix Z-score normalized matrix
#' @param batch_vector Batch vector used for normalization
#' @return List with diagnostic information per batch
check_zscore_normalization <- function(normalized_matrix, batch_vector) {
  
  unique_batches <- unique(batch_vector)
  results <- list()
  
  for (batch in unique_batches) {
    batch_indices <- batch_vector == batch
    batch_data <- normalized_matrix[, batch_indices, drop = FALSE]
    
    if (ncol(batch_data) > 1) {
      # Calculate statistics for this batch
      row_means <- rowMeans(batch_data)
      row_sds <- apply(batch_data, 1, sd)
      
      results[[batch]] <- list(
        n_cells = ncol(batch_data),
        mean_of_means = mean(row_means),
        mean_of_sds = mean(row_sds),
        max_abs_mean = max(abs(row_means)),
        mean_deviation_from_unit_variance = mean(abs(row_sds - 1)),
        properly_normalized = all(abs(row_means) < 1e-10) && all(abs(row_sds - 1) < 1e-10)
      )
    } else {
      results[[batch]] <- list(
        n_cells = 1,
        note = "Single cell batch - all values set to 0"
      )
    }
  }
  
  return(results)
}