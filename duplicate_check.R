args <- commandArgs(trailingOnly = TRUE)
location <- args[1]
filename <- args[2]

  # Read the data from the CSV file
  df <- read.table(paste0(location, filename, "_hwe.fam"), stringsAsFactors = FALSE)
  
  # Create a vector to keep track of seen values in the second column
  seen_values <- c()
  
  # Loop through each row of the data frame
  for (i in 1:nrow(df)) {
    current_value <- df[i, 2]  # Access the value in the second column
    
    # Check if the value has been seen before
    if (current_value %in% seen_values) {
      # If it has been seen, append "_repeat"
      df[i, 2] <- paste0(current_value, "_repeat_", i)
    } else {
      # If it's the first time seeing this value, add it to seen_values
      seen_values <- c(seen_values, current_value)
    }
  }
  
  # Save the updated data to a new CSV file
  write.table(df, paste0(location, filename, "_hwe.fam"), row.names = FALSE, col.names = FALSE, quote = FALSE)
