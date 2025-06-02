# Set the working directory to the script's location
setwd(getSrcDirectory(function(dummy) {dummy}))

# Set a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Required packages
packages <- c("tidyverse", "Synth", "readODS")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))

# Clear the console
cat("\014")

# Sheet type checker functions
is_germans_only <- function(sheet_name) {
  grepl("Germans_only", sheet_name)
}

is_german_language <- function(sheet_name) {
  grepl("German_language", sheet_name)
}

# Function to create directories for each analysis
create_analysis_directories <- function(analysis_name) {
  base_dir <- file.path("Results", paste0(analysis_name))
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "pdf"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "txt"), showWarnings = FALSE, recursive = TRUE)
  return(base_dir)
}

# Function to calculate RMSPE for a specific period
calculate_period_rmspe <- function(actual, synthetic, start_year, end_year, year_index) {
  period_indices <- which(year_index >= start_year & year_index <= end_year)
  differences <- actual[period_indices] - synthetic[period_indices]
  rmspe <- sqrt(mean(differences^2))
  return(rmspe)
}

# Function to save synthetic control details
save_synth_details <- function(synth_out, dataprep_out, author, base_dir, citation_var) {
  output_file <- file.path(base_dir, "txt", paste0(author, ".txt"))
  sink(output_file)
  cat(paste0("Synthetic Control Analysis for ", author, "\n"))
  cat(paste0("Using citation variable: ", citation_var, "\n\n"))
  
  # Write weights
  cat("Synthetic Control Weights:\n")
  solution_w <- synth_out$solution.w
  all_names <- dataprep_out$names.and.numbers$unit.names
  all_numbers <- dataprep_out$names.and.numbers$unit.numbers
  name_map <- setNames(all_names, all_numbers)
  control_numbers <- as.numeric(rownames(solution_w))
  
  weights_df <- data.frame(
    Name = name_map[as.character(control_numbers)],
    Weight = as.numeric(solution_w)
  )
  weights_df <- weights_df[order(-weights_df$Weight), ]
  weights_df$Weight <- sprintf("%.6f", weights_df$Weight)
  significant_weights <- weights_df[as.numeric(weights_df$Weight) > 0.001, ]
  print(significant_weights, row.names = FALSE)
  
  # Write predictor balance with modified special predictors
  cat("\nPredictor Balance:\n")
  dataprep_mod <- dataprep_out
  special_predictors <- grep("special", rownames(dataprep_mod$X0))
  dataprep_mod$X0[special_predictors,] <- dataprep_mod$X0[special_predictors,] * 1000000
  dataprep_mod$X1[special_predictors,] <- dataprep_mod$X1[special_predictors,] * 1000000
  dataprep_mod$Z0[special_predictors,] <- dataprep_mod$Z0[special_predictors,] * 1000000
  dataprep_mod$Z1[special_predictors,] <- dataprep_mod$Z1[special_predictors,] * 1000000
  
  synth.tables <- synth.tab(dataprep.res = dataprep_mod, synth.res = synth_out)
  print(synth.tables$tab.pred)
  
  # Calculate and write RMSPEs
  Y1plot <- dataprep_out$Y1plot
  synthetic <- dataprep_out$Y0plot %*% synth_out$solution.w
  years <- 1878:1932
  
  pre_rmspe <- calculate_period_rmspe(Y1plot[,1], synthetic[,1], 1878, 1916, years)
  post_rmspe <- calculate_period_rmspe(Y1plot[,1], synthetic[,1], 1917, 1932, years)
  rmspe_ratio <- post_rmspe / pre_rmspe
  
  cat("\nRoot Mean Square Prediction Error (RMSPE):")
  cat("\nPre (1878-1916):", pre_rmspe)
  cat("\nPost (1917-1932):", post_rmspe)
  cat("\nPost/Pre Ratio:", rmspe_ratio, "\n")
  
  sink()
  
  return(rmspe_ratio)
}

# Function to create plot
create_synth_plot <- function(synth_out, dataprep_out, citation_var) {
  Y1plot <- dataprep_out$Y1plot
  synthetic <- dataprep_out$Y0plot %*% synth_out$solution.w
  
  par(mai = c(0.6, 1.2, 0.2, 1.2),
      family = "sans",
      cex = 1.2,
      cex.axis = 1.2,
      cex.lab = 1.2,
      tck = 0.01,
      lwd = 0.8,
      las = 1,
      mgp = c(3.5, 0.8, 0))
  
  y_range <- range(c(Y1plot, synthetic))
  y_gap <- diff(y_range) * 0.25
  y_limits <- c(y_range[1], y_range[2] + y_gap)
  
  # Modify y-axis label based on citation variable
  y_label <- if(citation_var == "cite_German") {
    "N\uadgram share in German (%)"
  } else {
    "N\uadgram share in English (%)"
  }
  
  plot(1878:1932, Y1plot[, 1],
       ylim = y_limits,
       xlim = c(1878, 1932),
       type = "n",
       xlab = "",
       ylab = y_label,
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n")
  
  axis(1, at = seq(1880, 1930, by = 10), cex.axis = 1.2)
  
  y_ticks <- pretty(y_limits)
  y_labels <- format(y_ticks, scientific = TRUE)
  y_labels <- gsub("-", "\uad", y_labels)
  axis(2, at = y_ticks, labels = y_labels, las = 0, cex.axis = 1.2)
  
  box(lwd = 0.8)
  
  lines(1878:1932, Y1plot[, 1], lwd = 2)
  lines(1878:1932, synthetic[, 1], lty = 2, lwd = 2)
  
  abline(v = 1917, lty = 2, lwd = 0.8)
  
  legend("topright",
         legend = c("Actual", "Synthetic"),
         lty = c(1, 2),
         lwd = 1,
         bty = "n",
         cex = 1.2)
}

# Function to save plots in PDF format
save_synth_plots <- function(synth_out, dataprep_out, author, citation_var, base_dir) {
  pdf(file = file.path(base_dir, "pdf", paste0(author, ".pdf")),
      width = 9.2, height = 5)
  create_synth_plot(synth_out, dataprep_out, citation_var)
  dev.off()
}

# Function to check if all donors are German writers
check_all_german_writers <- function(data, donors) {
  all_german <- data %>%
    filter(Name %in% donors) %>%
    pull(wrote_German) %>%
    all(. == 1)
  return(all_german)
}

# Predictor selection logic
get_predictors <- function(data, donors, treated_author, sheet_name) {
  # Base predictors (always included)
  base_predictors <- c("YearofPublication",
                       "Socialist", 
                       "Political")
  
  # Add YearofTranslationtoEnglish for non-German language sheets
  if(!is_german_language(sheet_name)) {
    base_predictors <- c(base_predictors, "YearofTranslationtoEnglish")
  }
  
  # Add language indicators only if not Germans_only
  if(!is_germans_only(sheet_name)) {
    base_predictors <- c(base_predictors,
                         "wrote_English", "wrote_French",
                         "wrote_Greek", "wrote_Latin", 
                         "wrote_Italian", "wrote_Spanish",
                         "wrote_German")
  }
  
  return(base_predictors)
}

# Function to perform synthetic control analysis
perform_synth_analysis <- function(data, treated_author, all_potential_controls, base_dir, sheet_name) {
  tryCatch({
    predictors <- get_predictors(data, all_potential_controls, treated_author, sheet_name)
    character_vars <- predictors[grep("wrote_|Socialist|Political", predictors)]
    
    # Remove treated author from donor pool if present
    actual_donors <- setdiff(all_potential_controls, treated_author)
    
    # Set citation variable based on German_language sheet name
    citation_var <- if(is_german_language(sheet_name)) {
      "cite_German"
    } else {
      "cite_English"
    }
    
    # Check if the selected citation variable exists in the data
    if(!citation_var %in% names(data)) {
      warning(paste("Citation variable", citation_var, "not found for", treated_author, "- defaulting to cite_English"))
      citation_var <- "cite_English"
    }
    
    # Filter data for treated author and actual donors and required citation variable
    analysis_data <- data %>%
      filter(Name %in% c(treated_author, actual_donors)) %>%
      filter(!is.na(!!sym(citation_var)))
    
    # Only check YearofTranslationtoEnglish if using English citations and not German language
    if(citation_var == "cite_English" && !is_german_language(sheet_name) && 
       "YearofTranslationtoEnglish" %in% names(analysis_data)) {
      predictor_data <- analysis_data %>%
        group_by(Name) %>%
        filter(!is.na(YearofTranslationtoEnglish)) %>%
        ungroup()
      
      # Check if we still have the treated unit and at least one control after filtering
      if(!(treated_author %in% predictor_data$Name)) {
        stop("Treated unit filtered out due to missing YearofTranslationtoEnglish")
      }
      if(!any(actual_donors %in% predictor_data$Name)) {
        stop("All control units filtered out due to missing YearofTranslationtoEnglish")
      }
      
      analysis_data <- predictor_data
    }
    
    # Verify citations start at 0
    if(min(analysis_data[[citation_var]], na.rm = TRUE) != 0) {
      warning(paste("Citations do not start at 0 for", citation_var, "- adjusting values"))
      min_cite <- min(analysis_data[[citation_var]], na.rm = TRUE)
      analysis_data[[citation_var]] <- analysis_data[[citation_var]] - min_cite
    }
    
    # Prepare data
    analysis_data <- analysis_data %>%
      group_by(Name) %>%
      mutate(Author_ID = cur_group_id()) %>%
      ungroup()
    
    for (var in character_vars) {
      analysis_data[[var]] <- as.numeric(as.factor(analysis_data[[var]])) - 1
    }
    
    special_predictors <- list(
      list(citation_var, 1914:1916, "mean"),
      list(citation_var, 1908:1910, "mean"),
      list(citation_var, 1902:1904, "mean"),
      list(citation_var, 1896:1898, "mean"),
      list(citation_var, 1890:1892, "mean"),
      list(citation_var, 1884:1886, "mean"),
      list(citation_var, 1878:1880, "mean")
    )
    
    treated_id <- analysis_data %>%
      filter(Name == treated_author) %>%
      distinct(Author_ID) %>%
      pull(Author_ID)
    
    control_ids <- analysis_data %>%
      filter(Name != treated_author) %>%
      distinct(Author_ID) %>%
      pull(Author_ID)
    
    dataprep_out <- dataprep(
      foo = as.data.frame(analysis_data),
      predictors = predictors,
      predictors.op = "mean",
      special.predictors = special_predictors,
      dependent = citation_var,
      unit.variable = "Author_ID",
      time.variable = "Year",
      treatment.identifier = treated_id,
      controls.identifier = control_ids,
      time.predictors.prior = 1878:1916,
      time.optimize.ssr = 1878:1916,
      unit.names.variable = "Name",
      time.plot = 1878:1932
    )
    
    # Run synthetic control
    synth_out <- tryCatch({
      synth(data.prep.obj = dataprep_out)
    }, error = function(e) {
      warning(paste("Synth optimization error for", treated_author, ":", conditionMessage(e)))
      return(NULL)
    })
    
    # Handle failed convergence
    if(is.null(synth_out) || is.null(synth_out$solution.w)) {
      warning(paste("Optimization failed to converge for", treated_author, "- skipping"))
      cat(paste0("Failed to converge: ", treated_author, "\n"), 
          file = file.path(base_dir, "txt/failed_convergence.txt"), 
          append = TRUE)
      return(NA)
    }
    
    rmspe_ratio <- save_synth_details(synth_out, dataprep_out, treated_author, base_dir, citation_var)
    save_synth_plots(synth_out, dataprep_out, treated_author, citation_var, base_dir)
    
    return(rmspe_ratio)
    
  }, error = function(e) {
    cat("Error processing", treated_author, ":", conditionMessage(e), "\n")
    return(NA)
  })
}

# Function to calculate both versions of p-values from RMSPE ratios
calculate_p_values <- function(rmspe_ratios) {
  treated_ratio <- rmspe_ratios[1]  # Marx's ratio
  placebo_ratios <- rmspe_ratios[-1]
  
  # Remove NA values for calculation
  valid_placebo_ratios <- placebo_ratios[!is.na(placebo_ratios)]
  total_valid_placebos <- length(valid_placebo_ratios)
  
  # Version 1: Exclude Marx from denominator (original method)
  # P-value = (# placebos >= Marx's ratio) / (total # placebos)
  p_value_exclude_marx <- sum(valid_placebo_ratios >= treated_ratio, na.rm = TRUE) / total_valid_placebos
  
  # Version 2: Include Marx in denominator 
  # P-value = (# placebos >= Marx's ratio + 1) / (total # placebos + 1)
  p_value_include_marx <- (sum(valid_placebo_ratios >= treated_ratio, na.rm = TRUE) + 1) / (total_valid_placebos + 1)
  
  return(list(
    p_value_exclude_marx = p_value_exclude_marx,
    p_value_include_marx = p_value_include_marx,
    n_valid_placebos = total_valid_placebos,
    n_placebos_ge_marx = sum(valid_placebo_ratios >= treated_ratio, na.rm = TRUE)
  ))
}

# Main analysis function
main <- function() {
  # Read data from Data folder
  data <- read.csv("Data/all_authors_with_citations_and_indicators.csv",
                   header = TRUE, stringsAsFactors = FALSE)
  
  # Check for required columns without stopping
  required_cols <- c("cite_English", "cite_German")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if(length(missing_cols) > 0) {
    warning(paste("Missing columns:", paste(missing_cols, collapse=", "), 
                  "\nWill default to using cite_English where available"))
  }
  
  # Read donor pool data from Data folder
  donor_data <- read_ods("Data/SCM_donor_pool_analysis.ods")
  sheet_names <- list_ods_sheets("Data/SCM_donor_pool_analysis.ods")
  
  # Initialize results dataframe for p-values
  p_values_results <- data.frame(
    Analysis = character(),
    P_Value_Exclude_Marx = numeric(),
    P_Value_Include_Marx = numeric(),
    N_Valid_Placebos = numeric(),
    N_Placebos_GE_Marx = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each sheet
  for (sheet_name in sheet_names) {
    base_dir <- create_analysis_directories(sheet_name)
    
    sheet_data <- read_ods("Data/SCM_donor_pool_analysis.ods", sheet = sheet_name)
    donors <- sheet_data$Name[!is.na(sheet_data$Name)]
    
    cat(paste("\nProcessing sheet:", sheet_name))
    cat("\nNumber of donors found:", length(donors))
    cat("\nDonors:", paste(donors, collapse=", "), "\n")
    
    analysis_set <- c("Karl Marx", donors)
    
    # FIXED: Create full potential control pool that includes Marx for placebo tests
    # When analyzing Marx, use only donors as controls
    # When analyzing any donor (placebo test), use Marx + other donors as controls
    full_control_pool <- c("Karl Marx", donors)
    
    # Define required columns based on sheet name
    base_columns <- c("Name", "Year", "cite_English", "cite_German",
                      "YearofPublication",
                      "Socialist", "Political")
    
    # Add YearofTranslationtoEnglish for non-German language sheets
    if(!is_german_language(sheet_name)) {
      base_columns <- c(base_columns, "YearofTranslationtoEnglish")
    }
    
    # Add language variables only if not Germans_only
    if(!is_germans_only(sheet_name)) {
      base_columns <- c(base_columns,
                        "wrote_English", "wrote_French",
                        "wrote_Greek", "wrote_Latin", 
                        "wrote_Italian", "wrote_Spanish",
                        "wrote_German")
    }
    
    # Clean data
    cleaned_data <- data %>%
      select(any_of(base_columns)) %>%
      mutate(
        Year = as.integer(Year),
        cite_English = as.numeric(cite_English),
        cite_German = as.numeric(cite_German),
        YearofPublication = as.integer(YearofPublication)
      ) %>%
      filter(Year >= 1878, Year <= 1932)
    
    # Process all authors and collect RMSPE ratios
    rmspe_ratios <- sapply(analysis_set, function(author) {
      # FIXED: For Marx, use only donors as controls
      # For placebo tests (donors), use full control pool (Marx + other donors)
      if(author == "Karl Marx") {
        control_pool <- donors
      } else {
        control_pool <- full_control_pool
      }
      
      perform_synth_analysis(cleaned_data, author, control_pool, base_dir, sheet_name)
    })
    
    # Calculate both versions of p-values
    p_value_results <- calculate_p_values(rmspe_ratios)
    
    # Add to p-values results
    p_values_results <- rbind(p_values_results, 
                              data.frame(
                                Analysis = sheet_name,
                                P_Value_Exclude_Marx = p_value_results$p_value_exclude_marx,
                                P_Value_Include_Marx = p_value_results$p_value_include_marx,
                                N_Valid_Placebos = p_value_results$n_valid_placebos,
                                N_Placebos_GE_Marx = p_value_results$n_placebos_ge_marx,
                                stringsAsFactors = FALSE
                              ))
    
    # Save individual analysis results
    results_df <- data.frame(
      Name = analysis_set,
      RMSPE_ratio = rmspe_ratios,
      stringsAsFactors = FALSE
    )
    
    write.csv(results_df, 
              file.path(base_dir, paste0(sheet_name, "_RMSPE_ratios.csv")), 
              row.names = FALSE)
    
    # Print results for this analysis
    cat(paste("\nResults for", sheet_name, ":"))
    cat(paste("\nMarx RMSPE ratio:", round(rmspe_ratios[1], 4)))
    cat(paste("\nP-value (exclude Marx):", round(p_value_results$p_value_exclude_marx, 4)))
    cat(paste("\nP-value (include Marx):", round(p_value_results$p_value_include_marx, 4)))
    cat(paste("\nPlacebos >= Marx:", p_value_results$n_placebos_ge_marx, "/", p_value_results$n_valid_placebos))
    cat("\n")
  }
  
  # Save overall p-values results to Results folder
  write.csv(p_values_results, "Results/SCM_p_values.csv", row.names = FALSE)
  
  # Print summary
  cat("\n=== SUMMARY OF ALL ANALYSES ===\n")
  for(i in 1:nrow(p_values_results)) {
    cat(sprintf("\n%s:", p_values_results$Analysis[i]))
    cat(sprintf("\n  P-value (exclude Marx): %.4f", p_values_results$P_Value_Exclude_Marx[i]))
    cat(sprintf("\n  P-value (include Marx): %.4f", p_values_results$P_Value_Include_Marx[i]))
    cat(sprintf("\n  Valid placebos: %d", p_values_results$N_Valid_Placebos[i]))
  }
  cat("\n")
}

# Run the analysis
main()