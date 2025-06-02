# Set the working directory to the script's location
setwd(getSrcDirectory(function(dummy) {dummy}))

# Load required libraries
if (!require(ngramr, quietly = TRUE)) {
  install.packages("ngramr")
  library(ngramr)
}

cat("\014")

# Create Results directory if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Create a data frame with terms and their biographical information
terms_df <- data.frame(
  name = c(
    "Robert Owen", "Charles Fourier", "Sismondi", 
    "Thomas Carlyle", "Blanqui", "Rodbertus", 
    "Proudhon", "Moses Hess", "Bakunin", 
    "Karl Marx", "John Ruskin", "Ferdinand Lassalle", 
    "Henry George", "August Bebel", "Kropotkin", 
    "Lujo Brentano", "Eduard Bernstein", "Edward Bellamy", 
    "Oscar Wilde", "Sombart"
  ),
  birth_year = c(
    1771, 1772, 1773, 1795, 1805, 1805,
    1809, 1812, 1814, 1818, 1819, 1825,
    1839, 1840, 1842, 1844, 1850, 1850,
    1854, 1863
  ),
  death_year = c(
    1858, 1837, 1842, 1881, 1881, 1875,
    1865, 1875, 1876, 1883, 1900, 1864,
    1897, 1913, 1921, 1931, 1932, 1898,
    1900, 1941
  )
)

# Create a results data frame to store our findings
results <- data.frame(
  Name = terms_df$name,
  Birth_Year = terms_df$birth_year,
  Death_Year = terms_df$death_year,
  Before_Trend = numeric(nrow(terms_df)),
  After_Trend = numeric(nrow(terms_df)),
  Before_Avg = numeric(nrow(terms_df)),
  After_Avg = numeric(nrow(terms_df)),
  Death_Freq = numeric(nrow(terms_df)),
  Trend_Change = numeric(nrow(terms_df)),
  Pct_Change = numeric(nrow(terms_df)),
  Trend_Category = character(nrow(terms_df))
)

# Function to calculate linear trend
calculate_trend <- function(years, frequencies) {
  if(length(years) < 5) return(NA)
  model <- lm(frequencies ~ years)
  return(coef(model)[2])  # Return the slope
}

# Loop through each person
for(i in 1:nrow(terms_df)) {
  name <- terms_df$name[i]
  death_year <- terms_df$death_year[i]
  
  cat(sprintf("Fetching data for %s...\n", name))
  
  # Get data from Google Ngrams
  tryCatch({
    data <- ngram(
      phrases = name,
      corpus = "en-2019",
      year_start = death_year - 55,  # Extra 5 years to ensure we have enough data
      year_end = death_year + 55,    # Extra 5 years to ensure we have enough data
      smoothing = 3                  # Apply some smoothing to reduce noise
    )
    
    # Convert to per million for readability
    data$freq_million <- data$Frequency * 1000000
    
    # Get data for periods before and after death
    before_data <- data[data$Year >= (death_year - 50) & data$Year < death_year, ]
    after_data <- data[data$Year > death_year & data$Year <= (death_year + 50), ]
    
    # Find frequency at death year (or closest)
    death_idx <- which.min(abs(data$Year - death_year))
    results$Death_Freq[i] <- data$freq_million[death_idx]
    
    # Calculate trends
    results$Before_Trend[i] <- calculate_trend(before_data$Year, before_data$freq_million)
    results$After_Trend[i] <- calculate_trend(after_data$Year, after_data$freq_million)
    
    # Calculate averages
    results$Before_Avg[i] <- mean(before_data$freq_million)
    results$After_Avg[i] <- mean(after_data$freq_million)
    
    # Calculate trend change
    results$Trend_Change[i] <- results$After_Trend[i] - results$Before_Trend[i]
    
    # Calculate percentage change in trend
    if(!is.na(results$Before_Trend[i]) && abs(results$Before_Trend[i]) > 0.0000001) {
      results$Pct_Change[i] <- (results$Trend_Change[i] / abs(results$Before_Trend[i])) * 100
    }
    
    # Categorize trend changes
    before <- results$Before_Trend[i]
    after <- results$After_Trend[i]
    
    if(is.na(before) || is.na(after)) {
      results$Trend_Category[i] <- "Insufficient data"
    } else if(before > 0 && after > 0) {
      if(after > before) results$Trend_Category[i] <- "Accelerated growth"
      else results$Trend_Category[i] <- "Slowed growth"
    } else if(before < 0 && after < 0) {
      if(after < before) results$Trend_Category[i] <- "Accelerated decline"
      else results$Trend_Category[i] <- "Slowed decline"
    } else if(before <= 0 && after > 0) {
      results$Trend_Category[i] <- "Reversal (decline to growth)"
    } else if(before >= 0 && after < 0) {
      results$Trend_Category[i] <- "Reversal (growth to decline)"
    } else {
      results$Trend_Category[i] <- "No significant change"
    }
    
  }, error = function(e) {
    cat(sprintf("Error fetching data for %s: %s\n", name, e$message))
    # Leave all values as NA for this person
  })
}

# Function to format results for better readability - only format displayed columns
format_table <- function(results) {
  formatted <- results
  
  # Only format columns that will be displayed
  display_cols <- c("Before_Trend", "After_Trend", "Trend_Change")
  
  for(col in display_cols) {
    if(col %in% colnames(formatted)) {
      # For trend values (slopes), use 5 decimal places without scientific notation
      formatted[[col]] <- sprintf("%.5f", formatted[[col]])
    }
  }
  
  return(formatted)
}

# Sort results by death year for chronological viewing
results_by_death_year <- results[order(results$Death_Year),]
formatted_results <- format_table(results_by_death_year)

# Print interpretation guide in a more logical position
cat("\nInterpretation Guide:\n")
cat("- Before_Trend: Rate of change per year in 50 years before death\n")
cat("- After_Trend: Rate of change per year in 50 years after death\n")
cat("- Trend_Change: Difference between After_Trend and Before_Trend\n")
cat("- Positive Trend_Change: Improvement in trend after death\n")
cat("- Negative Trend_Change: Worsening trend after death\n")
cat("- Trend categories provide a qualitative description of the change\n")

# Print just the main results table sorted by death year, now including Birth_Year
cat("\nTrend Analysis by Death Year (50-year window):\n")
print(formatted_results[, c("Name", "Birth_Year", "Death_Year", "Before_Trend", "After_Trend", 
                            "Trend_Change", "Trend_Category")], row.names = FALSE)

# Save results to CSV in the Results directory
write.csv(results, "Results/historical_figures_trend_analysis.csv", row.names = FALSE)