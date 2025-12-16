# date created: 6/4/2022
# date modified: 14/7/2025
# author: Maren Büttner
# project: Convex Gating
# description: This R script runs the Hypergate algorithm to compute gates on 
# flow data, which come as adata objects from python

library(hypergate)
library(rhdf5)
library(jsonlite)

#' Load and preprocess HDF5 data for hypergate analysis
#' 
#' @param file_path Path to the HDF5 file
#' @param sample_key Key for sample information in the HDF5 file
#' @param cell_type_key Key for cell type information in the HDF5 file
#' @return List containing data_matrix, cellData, sample_levels, cell_types_levels
load_and_preprocess_data <- function(file_path, sample_key, cell_type_key) {
    # Read in adata object with the corresponding keys
    dset <- h5read(file_path, '/', compoundAsDataFrame = FALSE)
    
    # Process data to the correct format for hypergate
    # Get genes/cell ID
    cellID <- unlist(dset$obs$`_index`)
    featureID <- unlist(dset$var$`_index`)
    
    # Read data matrix
    data_matrix <- t(unlist(dset$X))
    colnames(data_matrix) <- featureID
    rownames(data_matrix) <- cellID
    
    # Get observations
    if (length(dset$obs[[sample_key]]) == 2) {
        sample <- unlist(dset$obs[[sample_key]]$codes)
        cell_types <- unlist(dset$obs[[cell_type_key]]$codes)
        sample_levels <- dset$obs[[sample_key]]$categories
        cell_types_levels <- dset$obs[[cell_type_key]]$categories
    } else {
        sample <- unlist(dset$obs[[sample_key]])
        cell_types <- unlist(dset$obs[[cell_type_key]])
        sample_levels <- dset$obs$`__categories`[[sample_key]]
        cell_types_levels <- dset$obs$`__categories`[[cell_type_key]]
    }
    
    # Convert observations into factors in a data frame
    cellData <- data.frame(
        sample = factor(sample_levels[sample + 1], levels = sample_levels),
        sample_code = sample + 1,
        cell_types = factor(cell_types_levels[cell_types + 1], levels = cell_types_levels),
        cell_types_code = cell_types + 1
    )
    
    return(list(
        data_matrix = data_matrix,
        cellData = cellData,
        sample_levels = sample_levels,
        cell_types_levels = cell_types_levels
    ))
}

#' Create train-test split for hypergate analysis
#' 
#' @param data_matrix Expression data matrix
#' @param cellData Cell metadata
#' @param sample Current sample being processed
#' @param cell_type Current cell type being processed
#' @param subsample_size Total number of cells in training set
#' @param train_ratio Fraction of target cells to use for training (e.g., 0.7 for 70% train, 30% test)
#' @param ratio_nontarget_target Ratio of non-target to target cells in training set
#' @return List containing training data, test data, cell names, and metadata
create_train_test_split <- function(data_matrix, cellData, sample, cell_type, 
                                  subsample_size, train_ratio = 0.7, 
                                  ratio_nontarget_target = 1) {
    
    # Validate parameters
    if (ratio_nontarget_target <= 0) {
        warning("Ratio of nontarget to target cells must be positive. Using default ratio of 1:1")
        ratio_nontarget_target <- 1
    }
    
    if (train_ratio <= 0 || train_ratio >= 1) {
        warning("Train ratio must be between 0 and 1. Using default ratio of 0.7")
        train_ratio <- 0.7
    }
    
    # Select cell type subset
    sample_id <- cellData$sample == sample
    cell_type_id <- cellData$cell_types == cell_type
    
    # Obtain target and non-target events
    sample_target <- sample_id & cell_type_id
    sample_nontarget <- sample_id & !cell_type_id
    
    # Subset data matrix
    data_sample <- data_matrix[sample_id, ]
    data_target <- data_matrix[sample_target, ]
    data_nontarget <- data_matrix[sample_nontarget, ]
    
    # Subset cellData
    cellData_target <- cellData[sample_target, ]
    cellData_nontarget <- cellData[sample_nontarget, ]
    
    # Step 1: Calculate training set composition based on subsample_size and ratio
    # subsample_size = n_target_train + n_nontarget_train
    # ratio_nontarget_target = n_nontarget_train / n_target_train
    # Therefore: n_target_train = subsample_size / (1 + ratio_nontarget_target)
    n_target_train <- floor(subsample_size / (1 + ratio_nontarget_target))
    n_nontarget_train <- subsample_size - n_target_train
    
    # Step 2: Calculate total target cells needed (train + test)
    # n_target_train = n_target_total * train_ratio
    # Therefore: n_target_total = n_target_train / train_ratio
    n_target_total_needed <- ceiling(n_target_train / train_ratio)
    n_target_test <- n_target_total_needed - n_target_train
    
    # Check availability
    n_available_target <- nrow(data_target)
    n_available_nontarget <- nrow(data_nontarget)
    
    if (n_target_total_needed > n_available_target) {
        warning(paste("Not enough target cells available. Needed:", n_target_total_needed, 
                     "Available:", n_available_target))
        n_target_total_needed <- n_available_target
        n_target_train <- floor(n_target_total_needed * train_ratio)
        n_target_test <- n_target_total_needed - n_target_train
        # Recalculate training set size
        n_nontarget_train <- min(floor(n_target_train * ratio_nontarget_target), n_available_nontarget)
        actual_subsample_size <- n_target_train + n_nontarget_train
        warning(paste("Adjusted training set size to:", actual_subsample_size))
    }
    
    if (n_nontarget_train > n_available_nontarget) {
        warning(paste("Not enough non-target cells available for training. Needed:", n_nontarget_train,
                     "Available:", n_available_nontarget))
        n_nontarget_train <- n_available_nontarget
    }
    
    # Generate random indices for reproducible splits
    set.seed(42)  # For reproducibility
    
    # Step 3: Subsample target cells
    target_subsample_indices <- sample(1:n_available_target, n_target_total_needed, replace = FALSE)
    data_target_sub <- data_target[target_subsample_indices, ]
    cellData_target_sub <- cellData_target[target_subsample_indices, ]
    
    # Step 4: Split target cells into train/test
    target_train_indices <- sample(1:n_target_total_needed, n_target_train, replace = FALSE)
    target_test_indices <- setdiff(1:n_target_total_needed, target_train_indices)
    
    # Extract train/test target data
    data_target_train <- data_target_sub[target_train_indices, ]
    data_target_test <- data_target_sub[target_test_indices, ]
    cellData_target_train <- cellData_target_sub[target_train_indices, ]
    cellData_target_test <- cellData_target_sub[target_test_indices, ]
    
    # Step 5: Sample non-target cells for training
    nontarget_train_indices <- sample(1:n_available_nontarget, n_nontarget_train, replace = FALSE)
    data_nontarget_train <- data_nontarget[nontarget_train_indices, ]
    cellData_nontarget_train <- cellData_nontarget[nontarget_train_indices, ]
    
    # Step 6: Use remaining non-target cells for test set
    nontarget_test_indices <- setdiff(1:n_available_nontarget, nontarget_train_indices)
    data_nontarget_test <- data_nontarget[nontarget_test_indices, ]
    cellData_nontarget_test <- cellData_nontarget[nontarget_test_indices, ]
    
    # Step 7: Combine train and test datasets
    data_train <- rbind(data_target_train, data_nontarget_train)
    cellData_train <- rbind(cellData_target_train, cellData_nontarget_train)
    
    data_test <- rbind(data_target_test, data_nontarget_test)
    cellData_test <- rbind(cellData_target_test, cellData_nontarget_test)
    
    # Step 8: Extract cell names for Python compatibility
    train_target_cell_names <- rownames(data_target_train)
    test_target_cell_names <- rownames(data_target_test)
    train_nontarget_cell_names <- rownames(data_nontarget_train)
    test_nontarget_cell_names <- rownames(data_nontarget_test)
    
    all_train_cell_names <- rownames(data_train)
    all_test_cell_names <- rownames(data_test)
    sample_cell_names <- rownames(data_sample)
    
    return(list(
        # Training data
        data_train = data_train,
        cellData_train = cellData_train,
        
        # Test data  
        data_test = data_test,
        cellData_test = cellData_test,
        
        # Full sample data for evaluation
        data_sample = data_sample,
        cell_type_id_sample = cell_type_id[sample_id],
        
        # Cell names for Python compatibility (instead of indices)
        train_target_cell_names = train_target_cell_names,
        test_target_cell_names = test_target_cell_names,
        train_nontarget_cell_names = train_nontarget_cell_names,
        test_nontarget_cell_names = test_nontarget_cell_names,
        all_train_cell_names = all_train_cell_names,
        all_test_cell_names = all_test_cell_names,
        sample_cell_names = sample_cell_names,
        
        # Summary statistics
        n_target_total = n_available_target,
        n_target_subsample = n_target_total_needed,
        n_target_train = n_target_train,
        n_target_test = n_target_test,
        n_nontarget_total = n_available_nontarget,
        n_nontarget_train = n_nontarget_train,
        n_nontarget_test = length(nontarget_test_indices),
        actual_training_size = nrow(data_train),
        requested_training_size = subsample_size,
        train_ratio = train_ratio,
        ratio_nontarget_target = ratio_nontarget_target
    ))
}

#' Run hypergate analysis on training data
#' 
#' @param split_data Output from create_train_test_split function
#' @param cell_type Current cell type being analyzed
#' @return List containing hypergate results and timing information
run_hypergate_analysis <- function(split_data, cell_type) {
    
    # Get cell type ID
    ct_ID <- which(levels(split_data$cellData_train$cell_types) == cell_type)
    
    # Run hypergate on training data
    start_time <- Sys.time()
    hg_output <- hypergate(
        xp = split_data$data_train,
        gate_vector = split_data$cellData_train$cell_types_code,
        level = ct_ID,
        verbose = FALSE
    )
    
    gating_predicted <- subset_matrix_hg(hg_output, split_data$data_train)
    training_time <- Sys.time() - start_time
    
    # Get human readable output
    hg_out_info <- hgate_info(
        hg_output,
        xp = split_data$data_train,
        gate_vector = split_data$cellData_train$cell_types_code,
        level = ct_ID
    )
    
    # Apply gate to all sample events
    start_time <- Sys.time()
    bm <- boolmat(gate = hg_output, xp = split_data$data_sample)
    application_time <- Sys.time() - start_time
    
    # Determine which cells are in the gate
    if (length(bm) > nrow(split_data$data_sample)) {
        in_gate <- rowSums(bm) == ncol(bm)  # All gates need to be true
    } else {
        in_gate <- bm
    }
    
    total_time <- training_time + application_time
    
    return(list(
        hg_output = hg_output,
        hg_out_info = hg_out_info,
        in_gate = in_gate,
        training_time = training_time,
        application_time = application_time,
        total_time = total_time
    ))
}

#' Compute performance metrics (precision, recall, F1)
#' 
#' @param true_labels True cell type labels
#' @param predicted_labels Predicted gate results
#' @return List containing precision, recall, and F1 score
compute_performance_metrics <- function(true_labels, predicted_labels) {
    
    # Create confusion matrix
    conf_matrix <- table(true_labels, predicted_labels)
    
    # Handle edge cases where confusion matrix might not be 2x2
    if (nrow(conf_matrix) < 2 || ncol(conf_matrix) < 2) {
        return(list(precision = NA, recall = NA, f1 = NA))
    }
    
    # Calculate metrics
    tp <- conf_matrix[2, 2]  # True positives
    fp <- sum(conf_matrix[1, 2])  # False positives
    fn <- sum(conf_matrix[2, 1])  # False negatives
    
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1 <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
    
    return(list(
        precision = precision,
        recall = recall,
        f1 = f1,
        confusion_matrix = conf_matrix
    ))
}

#' Save analysis results to CSV files
#' 
#' @param file_path Original HDF5 file path 
#' @param cell_type_key Cell type key used in analysis
#' @param sample Current sample
#' @param cell_type Current cell type
#' @param output_table Detailed results table
#' @param subsample_size Subsample size used
#' @param save_individual Whether to save individual sample/cell type results
#' @param save_path Path to save the results to. If not used 
save_results <- function(file_path, cell_type_key, sample, cell_type, 
                        output_table, subsample_size, save_individual = TRUE, save_path = NULL) {
    
    if (save_individual) {
        if (length(save_path) == 0) {
            save_path <- dirname(file_path)
        }
        data_name <- gsub('.h5ad', '', basename(file_path))
        
        filename <- paste0(
            save_path, '/',
            data_name, '_',
            sample, '_',
            'hypergate', '_',
            cell_type_key, '_',
            cell_type, '_',
            as.character(subsample_size), '.csv'
        )
        
        write.csv(x = output_table, file = filename, row.names = FALSE)
    }
}

#' Main pipeline function to run complete hypergate analysis
#' 
#' @param file_path Path to HDF5 file
#' @param sample_key Key for sample information
#' @param cell_type_key Key for cell type information
#' @param subsample_size Total number of cells in training set
#' @param train_ratio Fraction of target cells to use for training
#' @param ratio_nontarget_target Ratio of non-target to target cells in training set
#' @param save Whether to save results to files
#' @param save_path Path to save the data, defaults to file_path
#' @return Data frame with summary results
run_hypergate_pipeline <- function(file_path, sample_key, cell_type_key, 
                                 subsample_size, train_ratio = 0.7,
                                 ratio_nontarget_target = 1, save = FALSE, save_path=NULL) {
    
    # Load and preprocess data
    cat("Loading and preprocessing data...\n")
    processed_data <- load_and_preprocess_data(file_path, sample_key, cell_type_key)
    
    # Initialize results table
    res_table <- data.frame(
        sample = character(),
        cell_type = character(),
        set_size = character(),
        ratio = numeric(),
        train_ratio = numeric(),
        time = numeric(),
        score = character(),
        value = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Store train/test indices for Python compatibility
    train_test_indices <- list()
    
    # Iterate over all samples and cell types
    for (sample in processed_data$sample_levels) {
        for (cell_type in processed_data$cell_types_levels) {
            
            # Skip not annotated cell types
            if (cell_type == 'not annotated') {
                next
            }
            
            cat(paste("Processing sample:", sample, "cell type:", cell_type, "\n"))
            
            # Create train-test split
            split_data <- create_train_test_split(
                processed_data$data_matrix,
                processed_data$cellData,
                sample,
                cell_type,
                subsample_size,
                train_ratio,
                ratio_nontarget_target
            )
            
            # Store cell names for Python compatibility
            split_key <- paste(sample, cell_type, sep = "_")
            train_test_indices[[split_key]] <- list(
                train_target_cell_names = split_data$train_target_cell_names,
                train_nontarget_cell_names = split_data$train_nontarget_cell_names,
                test_target_cell_names = split_data$test_target_cell_names,
                test_nontarget_cell_names = split_data$test_nontarget_cell_names,
                all_train_cell_names = split_data$all_train_cell_names,
                all_test_cell_names = split_data$all_test_cell_names,
                sample_cell_names = split_data$sample_cell_names,
                train_ratio = train_ratio,
                ratio_nontarget_target = ratio_nontarget_target,
                n_target_train = split_data$n_target_train,
                n_target_test = split_data$n_target_test,
                n_nontarget_train = split_data$n_nontarget_train,
                n_nontarget_test = split_data$n_nontarget_test,
                actual_training_size = split_data$actual_training_size,
                requested_training_size = subsample_size
            )
            
            # Run hypergate analysis
            hg_results <- run_hypergate_analysis(split_data, cell_type)
            
            # Compute performance on full sample
            full_metrics <- compute_performance_metrics(
                split_data$cell_type_id_sample,
                hg_results$in_gate
            )
            
            # Compute performance on training set
            train_mask <- rownames(split_data$data_sample) %in% rownames(split_data$data_train)
            train_metrics <- compute_performance_metrics(
                split_data$cell_type_id_sample[train_mask],
                hg_results$in_gate[train_mask]
            )
            
            # Compute performance on test set
            test_mask <- rownames(split_data$data_sample) %in% rownames(split_data$data_test)
            test_metrics <- compute_performance_metrics(
                split_data$cell_type_id_sample[test_mask],
                hg_results$in_gate[test_mask]
            )
            
            # Create detailed output table
            cellID <- rownames(split_data$data_sample)
            train_info <- train_mask
            test_info <- test_mask
            output_table <- data.frame(
                cellID = cellID,
                in_gate = hg_results$in_gate,
                true_label = split_data$cell_type_id_sample,
                train_set = train_info,
                test_set = test_info,
                stringsAsFactors = FALSE
            )
            
            # Add results to summary table - Full sample
            full_results <- data.frame(
                sample = rep(sample, 3),
                cell_type = rep(cell_type, 3),
                set_size = rep('full', 3),
                ratio = rep(ratio_nontarget_target, 3),
                train_ratio = rep(train_ratio, 3),
                time = rep(as.numeric(hg_results$total_time), 3),
                score = c('precision', 'recall', 'f1'),
                value = c(full_metrics$precision, full_metrics$recall, full_metrics$f1),
                stringsAsFactors = FALSE
            )
            
            # Add results to summary table - Training set
            train_results <- data.frame(
                sample = rep(sample, 3),
                cell_type = rep(cell_type, 3),
                set_size = rep('train', 3),
                ratio = rep(ratio_nontarget_target, 3),
                train_ratio = rep(train_ratio, 3),
                time = rep(as.numeric(hg_results$training_time), 3),
                score = c('precision', 'recall', 'f1'),
                value = c(train_metrics$precision, train_metrics$recall, train_metrics$f1),
                stringsAsFactors = FALSE
            )
            
            # Add results to summary table - Test set
            test_results <- data.frame(
                sample = rep(sample, 3),
                cell_type = rep(cell_type, 3),
                set_size = rep('test', 3),
                ratio = rep(ratio_nontarget_target, 3),
                train_ratio = rep(train_ratio, 3),
                time = rep(as.numeric(hg_results$application_time), 3),
                score = c('precision', 'recall', 'f1'),
                value = c(test_metrics$precision, test_metrics$recall, test_metrics$f1),
                stringsAsFactors = FALSE
            )
            
            res_table <- rbind(res_table, full_results, train_results, test_results)
            
            # Save individual results if requested
            if (save != FALSE) {
                if (length(save_path) == 0) {
                    save_path <- dirname(file_path)
                    }
                save_results(
                    file_path, cell_type_key, sample, cell_type,
                    output_table, split_data$actual_training_size, save_path=save_path
                )
            }
        }
    }
    
    # Save summary results and indices
    if (save != FALSE) {
        if (length(save_path) == 0) {
            save_path <- dirname(file_path)
        }
        data_name <- gsub('.h5ad', '', basename(file_path))
        
        # Save summary table
        summary_filename <- paste0(save_path, '/', data_name, '_', cell_type_key, '_hypergate.csv')
        write.csv(x = res_table, file = summary_filename, row.names = FALSE)
        
        # Save train/test cell names as JSON for Python compatibility
        cell_names_filename <- paste0(save_path, '/', data_name, '_', cell_type_key, '_train_test_cell_names.json')
        write_json(train_test_indices, cell_names_filename, pretty = TRUE)
        
        cat(paste("Results saved to:", summary_filename, "\n"))
        cat(paste("Train/test cell names saved to:", cell_names_filename, "\n"))
    }
    
    # Return results and cell names
    return(list(
        results = res_table,
        train_test_cell_names = train_test_indices
    ))
}