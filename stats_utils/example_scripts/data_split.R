#############
# Setup
library(dplyr)
#############

#############
# Create a fake screen_data
# Continuous variables:
sample_size <- 200
x <- rnorm(n = sample_size, mean = 20, sd = 2)
y <- rnorm(n = sample_size, mean = 100, sd = 30)
z <- rnorm(n = sample_size, mean = 50, sd = 10)

# Binomial (for eg gender):
gender <- rbinom(sample_size, size = 1, prob = 0.50)
# Good to check, should be very close to 50%:
sum(length(which(gender == 0))) / sample_size

# Create the fake data frame:
screen_data <- data.frame(x, y, z, gender)
# Explore it:
class(screen_data)
str(screen_data)
dim(screen_data)
head(screen_data)
tail(screen_data)
#############

#############
# Create a training data set with random sampling:
data_train <- sample_frac(screen_data, 0.6)
# Explore:
class(data_train)
str(data_train)
dim(data_train)
head(data_train)
tail(data_train)

# Create a test set with the remaining rows:
train_index <- as.numeric(rownames(data_train))
train_index
data_test <- screen_data[ -train_index, ]
# Explore:
class(data_test)
str(data_test)
dim(data_test)
head(data_test)
tail(data_test)
# Sanity checks to make sure everything is OK:
which(rownames(data_test) %in% rownames(data_train)) # Should be zero (all FALSE)
length(which(rownames(data_test) %in% rownames(screen_data))) # Should be equal to fraction remaining (ie 40% here)
dim(screen_data)

# All good except we have subset based on rownames, which are not indices!
# Reset them if they are of no importance:
# (This is a quick fix but not the most appropriate)
rownames(data_test)
rownames(data_test) <- NULL
rownames(data_test)
#############

#############
# Divide the test set in two, one for model validation, second for actual model test:
model_validating <- sample_frac(data_test, 0.5)
# Explore:
class(model_validating)
str(model_validating)
dim(model_validating)
head(model_validating)
tail(model_validating)
# Sanity checks to make sure everything is OK:
length(which(rownames(model_validating) %in% rownames(data_test))) # Should be equal to fraction sampled
dim(data_test)
#############

#############
# Get the indices and split as above:
valid_index <- as.numeric(rownames(model_validating))
valid_index
length(valid_index)
# Use the same code as above, use "-" not "!".
# ! is for logicals, while - gives negatives indices
# https://stackoverflow.com/questions/12328056/how-do-i-delete-rows-in-a-data-frame
# data_test2 <- data_test[!(data_test$substudy_part_id %in% valid_index), ] 
data_test2 <- data_test[ -valid_index, ]
# Explore:
class(data_test2)
str(data_test2)
dim(data_test2)
head(data_test2)
tail(data_test2)
# Sanity checks to make sure everything is OK:
length(which(rownames(data_test2) %in% rownames(model_validating))) # Should be zero
length(which(rownames(data_test2) %in% rownames(screen_data))) # Should be equal to fraction remaining of sampled
dim(model_validating)
#############

#############
# Summary and last check:
dim(data_train)
dim(model_validating)
dim(data_test2)

my_splits <- c(data_train, model_validating, data_test2)
sapply(my_splits, summary)
boxplot(data_train$x, model_validating$x, data_test2$x)
cor(model_validating$x, data_test2$x)
plot(model_validating$x, data_test2$x)
# values are random though, in a real data set all these should match
# (and other consistency checks)
#############