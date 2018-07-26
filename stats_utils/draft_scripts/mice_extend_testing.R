# Get data:
# https://stefvanbuuren.name/fimd/sec-knowledge.html
library(mice)
head(boys)
data <- boys[, c("age", "hgt", "wgt", "hc", "reg")]
imp <- mice(data, print = FALSE, seed = 71712)
# imp40 <- mice.mids(imp, maxit = 2, print = F)
# imp40$iteration
long <- mice::complete(imp, "long", include = TRUE)
long$whr <- with(long, 100 * wgt / hgt)
head(long)
dim(long)
imp.itt <- as.mids(long)
# Check it works:
imp.itt_2 <- mice.mids(imp.itt,
											 maxit = 20, # max iterations per imputation
											 print = F # omit printing of the iteration cycle
)
imp$iteration
imp.itt_2$iteration
plot(imp.itt_2)

# Write to file:
fwrite(long, 'imp.itt_long.tsv',
			 sep = '\t',
			 na = 'NA',
			 col.names = TRUE,
			 row.names = FALSE, # don't keep, IDs are in .id
			 quote = FALSE,
			 nThread = num_cores # use all cores available to write out
)

# read back in:
long <- 'imp.itt_long.tsv'
input_data <- fread(long,
										sep = '\t',
										header = TRUE,
										stringsAsFactors = FALSE)

input_data <- as.mids(as.data.frame(input_data))
head(input_data)
class(input_data)

# Try again:
long_from_file <- mice.mids(input_data,
											 maxit = 20, # max iterations per imputation
											 print = F # omit printing of the iteration cycle
)
plot(long_from_file)
