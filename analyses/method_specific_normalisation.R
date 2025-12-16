

# many of the methods simply do scanpy style cell normalise and log normalise. 
normalise_scanpy <- function(count_matrix, cellnorm = T, lognorm = T, zscore = F, batch_vector = NULL, target_sum = NULL) {
  
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
  

  if(cellnorm){
    
    # Step 1: Normalize total counts per cell (sc.pp.normalize_total)
    # Calculate sum per column (cell) - features are on rows, cells on columns
    col_sums <- colSums(count_matrix)
    
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
    result <- count_matrix
    # Normalize each cell to target_sum
    # This replicates: each cell's counts / cell_total * target_sum
    for (j in 1:ncol(result)) {
      if (col_sums[j] > 0) {
        result[, j] <- result[, j] * (target_sum / col_sums[j])
      }
    }
  }else{
    result <- count_matrix
  }
  
  # Step 2: Apply log1p transformation (sc.pp.log1p)
  # This is log(x + 1) transformation
  result <- log1p(result)
  
  
  if(zscore){
    
    if (length(batch_vector) != ncol(count_matrix)) {
      stop("batch_vector length must match number of columns in count_matrix")
    }
    
    # Convert to matrix if data.frame
    if (is.data.frame(result)) {
      result <- as.matrix(result)
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
    
    
  }
  
  return(result)
}



#scTranslator does min/max adjustment

normalise_scTranslator <- function(x){
  MIN = min(x)
  MAX = max(x)
  newx = 1e-8 +  (x-MIN) / (MAX - MIN )  * (1 - 1e-8)
}


# cTPnet does a slightly different log transform
normalise_cTPnet <- function(count_matrix, axis = 1) {
  
  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }
  
  if (!axis %in% c(0, 1)) {
    stop("axis must be 0 (genes) or 1 (cells)")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  result <- count_matrix

  # Add 1 to all values (x + 1.0)
  x_plus_1 <- result + 1.0
  
  if (axis == 1) {
    # Compute geometric mean along columns 

    
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
    
  } else if (axis == 0) {
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

