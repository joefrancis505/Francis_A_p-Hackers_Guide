# Set the working directory to the script's location
setwd(getSrcDirectory(function(dummy) {dummy}))

# Set a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
if (!require(tidyverse, quietly = TRUE)) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require(Synth, quietly = TRUE)) {
  install.packages("Synth")
  library(Synth)
}

# Clear the console
cat("\014")

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
  
  y_label <- "N\uadgram share in English (%)"
  
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

# Function to save plots in PDF format only
save_synth_plots <- function(synth_out, dataprep_out, author, citation_var, base_dir) {
  pdf(file = file.path(base_dir, "pdf", paste0(author, ".pdf")),
      width = 9.2, height = 5)
  create_synth_plot(synth_out, dataprep_out, citation_var)
  dev.off()
}

# Function to perform synthetic control analysis with simplified predictors
perform_synth_analysis <- function(data, treated_author, donors, base_dir, analysis_name) {
  tryCatch({
    # Basic predictors as requested
    predictors <- c("YearofPublication", "wrote_English", "wrote_German", "YearofTranslationtoEnglish")
    
    # Citation variable - using cite_English only as specified
    citation_var <- "cite_English"
    
    # Remove treated author from donor pool if present
    actual_donors <- setdiff(donors, treated_author)
    
    # Filter data for treated author and actual donors
    analysis_data <- data %>%
      filter(Name %in% c(treated_author, actual_donors)) %>%
      filter(!is.na(!!sym(citation_var)))
    
    # Verify we have YearofTranslationtoEnglish
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
    
    # Special predictors for specific time periods as requested
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

# Function to calculate p-values from RMSPE ratios
calculate_p_values <- function(rmspe_ratios) {
  treated_ratio <- rmspe_ratios[1]  # Karl Marx's ratio (always first in the list)
  placebo_ratios <- rmspe_ratios[-1]
  total_placebos <- length(placebo_ratios)
  p_value <- sum(placebo_ratios >= treated_ratio, na.rm = TRUE) / total_placebos
  return(p_value)
}

# Main analysis function
main <- function() {
  # Read data
  data <- read.csv("Data/all_authors_with_citations_and_indicators.csv",
                   header = TRUE, stringsAsFactors = FALSE)
  
  # Check for required column
  if(!"cite_English" %in% names(data)) {
    stop("Required column 'cite_English' not found in the data")
  }
  
  # First filter to only Socialist authors and remove Marx entirely
  socialist_authors <- data %>%
    filter(Socialist == 1, Name != "Marx", Name != "Karl Marx") %>%
    pull(Name) %>%
    unique()
  
  # Define the analyses
  analyses <- list(
    "Socialist" = list(
      donors = socialist_authors
    ),
    "Socialist_excluding_Wilde" = list(
      donors = socialist_authors[socialist_authors != "Oscar Wilde"]
    ),
    "Socialist_excluding_Sombart_and_Wilde" = list(
      donors = socialist_authors[!socialist_authors %in% c("Sombart", "Oscar Wilde")]
    ),
    "Socialist_excluding_Bebel_Sombart_and_Wilde" = list(
      donors = socialist_authors[!socialist_authors %in% c("Sombart", "Oscar Wilde", "August Bebel")]
    ),
    "Socialist_excluding_Bebel_George_Sombart_and_Wilde" = list(
      donors = socialist_authors[!socialist_authors %in% c("Sombart", "Oscar Wilde", "August Bebel", "Henry George")]
    )
  )
  
  # Initialize results dataframe for p-values
  p_values_results <- data.frame(
    Analysis = character(),
    P_Value = numeric(),
    N_Donors = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each analysis
  for (analysis_name in names(analyses)) {
    base_dir <- create_analysis_directories(analysis_name)
    donors <- analyses[[analysis_name]]$donors
    
    cat(paste("\nProcessing analysis:", analysis_name))
    cat("\nNumber of donors found:", length(donors))
    cat("\nDonors:", paste(donors, collapse=", "), "\n")
    
    # Set the treated author as Karl Marx
    treated_author <- "Karl Marx"
    
    # Clean data - select only the columns we need
    cleaned_data <- data %>%
      select(Name, Year, cite_English, YearofPublication, 
             wrote_English, wrote_German, YearofTranslationtoEnglish, Socialist) %>%
      mutate(
        Year = as.integer(Year),
        cite_English = as.numeric(cite_English),
        YearofPublication = as.integer(YearofPublication),
        wrote_English = as.integer(wrote_English),
        wrote_German = as.integer(wrote_German),
        Socialist = as.integer(Socialist)
      ) %>%
      filter(Year >= 1878, Year <= 1932)
    
    # Process Karl Marx first (treated analysis)
    marx_rmspe <- perform_synth_analysis(cleaned_data, "Karl Marx", donors, base_dir, analysis_name)
    
    # Process placebo tests with Marx included in donor pool
    placebo_donors <- c("Karl Marx", donors)
    placebo_rmspes <- sapply(donors, function(author) {
      perform_synth_analysis(cleaned_data, author, placebo_donors, base_dir, analysis_name)
    })
    
    # Combine results
    rmspe_ratios <- c(marx_rmspe, placebo_rmspes)
    names(rmspe_ratios) <- c("Karl Marx", donors)
    
    # Calculate p-value with Karl Marx as the treated unit
    p_value <- calculate_p_values(rmspe_ratios)
    
    # Add to p-values results
    p_values_results <- rbind(p_values_results, 
                              data.frame(
                                Analysis = analysis_name,
                                P_Value = p_value,
                                N_Donors = length(donors),
                                stringsAsFactors = FALSE
                              ))
    
    # Save individual analysis results
    results_df <- data.frame(
      Name = c("Karl Marx", donors),
      RMSPE_ratio = rmspe_ratios,
      stringsAsFactors = FALSE
    )
    
    write.csv(results_df, 
              file.path(base_dir, paste0(analysis_name, "_RMSPE_ratios.csv")), 
              row.names = FALSE)
  }
  
  # Save overall p-values results
  write.csv(p_values_results, "Results/SCM_p_values_socialists.csv", row.names = FALSE)
}

# Run the analysis
main()