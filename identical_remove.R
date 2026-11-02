args <- commandArgs(trailingOnly = TRUE)
location <- args[1]
filename <- args[2]

# Define file paths
genome_file <- paste0(location, "pihat_min_iden.genome")
imiss_file <- paste0(location, filename, ".imiss")
output_file <- paste0(location, "identical_remove.txt")


# Read the .genome file
genome_data <- read.table(genome_file, header = TRUE, stringsAsFactors = FALSE)

# Read the .imiss file
imiss_data <- read.table(imiss_file, header = TRUE, stringsAsFactors = FALSE)

# Extract FID and IID from both individuals in the .genome file
individuals <- genome_data[, c("FID1", "IID1", "FID2", "IID2")]

# Get the call rate (F_MISS) for the individuals in the .imiss file
# Create a function to extract call rates based on FID and IID
get_call_rate <- function(fid, iid, imiss_data) {
  idx <- which(imiss_data$FID == fid & imiss_data$IID == iid)
  if (length(idx) == 1) {
    return(1 - imiss_data$F_MISS[idx])  # Call rate = 1 - F_MISS
  } else {
    return(NA)  # If no match found, return NA
  }
}

# Initialize vectors to store the Family ID and Individual ID of the individual with the lower call rate
lower_call_fid <- character(0)
lower_call_iid <- character(0)

# Loop through each pair of individuals and compare their call rates
for (i in 1:nrow(individuals)) {
  call_rate1 <- get_call_rate(individuals$FID1[i], individuals$IID1[i], imiss_data)
  call_rate2 <- get_call_rate(individuals$FID2[i], individuals$IID2[i], imiss_data)
  
  # Compare the call rates and store the one with the lower call rate
  if (!is.na(call_rate1) & !is.na(call_rate2)) {
    if (call_rate1 < call_rate2) {
      lower_call_fid <- c(lower_call_fid, individuals$FID1[i])
      lower_call_iid <- c(lower_call_iid, individuals$IID1[i])
    } else {
      lower_call_fid <- c(lower_call_fid, individuals$FID2[i])
      lower_call_iid <- c(lower_call_iid, individuals$IID2[i])
    }
  }
}

# Combine the results into a data frame
results <- data.frame(FID = lower_call_fid, IID = lower_call_iid, stringsAsFactors = FALSE)

# Save the results to the output file without headers
write.table(results, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
