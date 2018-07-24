#!/usr/bin/env Rscript

######################
# R script to run with docopt for command line options:
'
run_mice_impute.R
=====================

Author: Antonio J Berlanga-Taylor
Release: |version|
Date: |today|

Purpose
========

This is a wrapper for the R package mice for multiple imputation.
See:
https://cran.r-project.org/web/packages/mice/index.html
https://stefvanbuuren.name/fimd/


Usage and options
=================


To run, type:
  Rscript run_mice_impute.R -I <INPUT_FILE> [options]

Usage: run_mice_impute.R (-I <INPUT_FILE>)
			 run_mice_impute.R [options]
			 run_mice_impute.R [-h | --help]

Options:
-I <INPUT_FILE>                 Input file name
-O <OUTPUT_FILE>                Output file name
--session=<R_SESSION_NAME>      R session name if to be saved. [default: FALSE]
-h --help                       Show this screen
--ID_col <int>                  Column with ID values for rows. [default: 1]
--vars_interest <FILE_VARS>     Variables of interest for sanity checking
--num_cores <int>               Number of cores to use for parallelising.
--derived_vars <FILE_EXCLUDE>   Derived variables to exclude from imputation
--missingness_cut <int>         Percent of missigness to exclude. [default: 30]
--mincor <num>                  Minimum correlation for quickpred. [default: 0.3]
--seed <int>                    Random seed number for reproducible analysis. [default: 123456]
-m <int>                        Number of imputed datasets. [default: 5]
--maxit <int>                   Number of iterations per dataset to impute. [default: 30]
--pred <FILENAME>               Predictor matrix to determines relationship between variables


Input:

Input file: A tab separated file with headers. This is usually a phenotype file with individuals (samples)
as rows and phenotype information (age, gender, etc.) as columns.
The file is read with data.table and stringsAsFactors = FALSE
A per column and per row cleaned data frame is expected with the first column as the IDs.

Variables of interest: A tab separated file with no header containing one column with variable names per row is
required. Variable names must match the column names in the input file.
The first variable will be used as the main variable of interest. Several plots
and diagnostic tests are performed using these.

Derived variables: (optional) a file in the same format as --vars_interest with
variable names that are derived, not needed for imputation and which should be excluded
from this step (an re-calculated manually afterwards). See below for more information.

--missingness_cut is used for information purposes only, no action is taken. You should
clean per row and per column first. The value passed is used for both rows and columns.


Output:

  A tab separated file with headers with the first imputed dataset.
  A tab separated file with headers in the long format with all imputed datasets. This
can then be read back into R and converted to a mids object with as.mids().
  With --session, an Rdata object that contains all the imputed (mids object) datasets,
this can then be re-loaded into R and ran in pooled analysis (slow). The session only saves
the imputed object.
  Several diagnostic plots (convergence, strip plots, bwplots, density plots).

Notes:

# Parallelising

The script tries to use multiple cores if available. Using all cores can slow
down the computer considerably so the number of cores is passed as
max(1, detectCores() - 1). If num_cores is provided there is no checking.

num_cores = 1 will still use the parallel setup code so expect it to be
slower than code without. If you want to use everything you have simply pass
num_cores as the number of cores available.

The cluster type is set to FORK, on Windows this probably will not run.

Parallelising will create num_cores * m imputed datasets.

# Predictor matrix

You can change the predictor matrix which specifies the variables that are used to
impute each incomplete variable with
pred

Running:
# ini <- mice(input_data, maxit = 0, print = F)
# ini$pred

will output the predictor matrix.
A value 1 indicates that the column variable was used to impute the row variable.
A row with all zeros indicates no missing values in that variable
You can change it to suit the specific relationships you need

By default, this script uses the option quickpred()
for a quick selection procedure of predictors based on
variable correlations of at least p = 0.30 (set using the option mincor).

# Variables of interest for basic sanity checking

You need to supply a variable of interest and some covariates to the
vars_interest
option.
The first variable will be taken as the variable of interest.

# Changing the method used for imputation

Check the options for defaultMethod in mice. Currently you cannot change these
with this script (make a copy of the script and modify separately if needed).

# Imputing derived variables

This is not straightforward. For an overview see:
https://stefvanbuuren.name/fimd/sec-knowledge.html

In particular, note that if derived variables are not needed for imputation,
exclude them from imputation, impute the originals, and then
calculate the derived variables. This is the approach used here, known as
impute, then transform (Von Hippel 2009). Provide the --derived_vars for this,
these will be excluded from imputation.

You will then need to recalculate them for the missing data.
eg if BMI is present and not needed for imputation (ie it provides no additional
information useful to impute other variables), pass it in the derived_vars option,
this will be excluded for the imputation procedure.
Once this is done, use the output from this script and run:

(this is from https://stefvanbuuren.name/fimd/sec-knowledge.html):
data <- boys[, c("age", "hgt", "wgt", "hc", "reg")]
imp <- mice(data, print = FALSE, seed = 71712)
long <- mice::complete(imp, "long", include = TRUE)
long$whr <- with(long, 100 * wgt / hgt)
imp.itt <- as.mids(long)

The first three steps are carried out here, the last need to be done manually
according to your needs.

If you are planning to carry out a model which includes interactions,
non-linear effects of variables or Cox proportional hazards model (which is
non-linear), this script and the MICE approach used here will not be appropriate
as the imputation model and the subsequent analysis model will be incompatible.
In this case you may need to use the SMC-FCS approach. See:
https://cran.r-project.org/web/packages/smcfcs/vignettes/smcfcs-vignette.html

# Saving the imputed R object

You can save the imputed object (the result of mice()) with --session
This may be slow though but required in order to run pooled analyses
using the multiple imputed datasets.


Documentation
=============

For more information see:

|url|
' -> doc

# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)
# See:
# https://cran.r-project.org/web/packages/docopt/docopt.pdf
# docopt(doc, args = commandArgs(TRUE), name = NULL, help = TRUE,
# version = NULL, strict = FALSE, strip_names = !strict,
# quoted_args = !strict)

# Print to screen:
str(args)
######################

######################
##########
# References

# Textbook online:
# Flexible Imputation of Missing Data. Second Edition.
# https://stefvanbuuren.name/fimd/

# Review (medicine/epidemiology):
# Multiple imputation for missing data in epidemiological and clinical research:
# potential and pitfalls
# https://www.bmj.com/content/338/bmj.b2393

# An example from epidemiology and pitfalls of imputation:
# Derivation and validation of QRISK, a new cardiovascular disease risk score for the United Kingdom: prospective open cohort study | The BMJ
# https://www.bmj.com/content/335/7611/136?ijkey=08580dfe6855e093441ea71ad418c11325310acd&keytype2=tf_ipsecsha
# Doubts about QRISK score: total / HDL cholesterol should be important. | The BMJ
# https://www.bmj.com/rapid-response/2011/11/01/doubts-about-qrisk-score-total-hdl-cholesterol-should-be-important
# QRISK - authors response | The BMJ
# https://www.bmj.com/rapid-response/2011/11/01/qrisk-authors-response


# Multiple imputation software
# http://www.stefvanbuuren.nl/mi/Software.html
#	Multivariate Imputation by Chained Equations
# http://stefvanbuuren.github.io/mice/
# Multivariate Imputation by Chained Equations
# https://github.com/stefvanbuuren/mice
# Ad hoc methods and mice
# https://gerkovink.github.io/miceVignettes/Ad_hoc_and_mice/Ad_hoc_methods.html

# Tutorials, etc. (use mice vignettes first):
#	Imputing Missing Data with R; MICE package | DataScience+
# https://datascienceplus.com/imputing-missing-data-with-r-mice-package/
# MICE: Multivariate Imputation by Chained Equations in R - MICE in R - Draft.pdf
#	http://www.stefvanbuuren.nl/publications/MICE%20in%20R%20-%20Draft.pdf
# Multiple Imputation for General Missing Data Patterns in the Presence of High-dimensional Data | Scientific Reports
# https://www.nature.com/articles/srep21689

# Easy tutorial using library mice:
# http://web.maths.unsw.edu.au/~dwarton/missingDataLab.html

# https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/
# https://www.jstatsoft.org/article/view/v045i03

# Parallel computation:
# https://stefvanbuuren.name/fimd/parallel-computation.html
# https://stackoverflow.com/Questions/24040280/Parallel-Computation-of-Multiple-Imputation-by-Using-Mice-R-Package

# Suggestions on how to carry out imputation:
# https://stefvanbuuren.name/fimd/conclusion-5.html
# https://stats.stackexchange.com/questions/219013/how-do-the-number-of-imputations-the-maximum-iterations-affect-accuracy-in-mul

# Some notes:
# Roughly one imputation per percent of incomplete data (White et al.,2011),
# but the more the better, 100 can easily be run on small datasets on a laptop
# Roughly 20-30 iterations should be enough, use plot() to check convergence:
# http://stats.stackexchange.com/questions/219013/how-do-the-number-of-imputations-the-maximum-iterations-affect-accuracy-in-mul

# Some studies imputed up to 70% of missing data, seems dangerous but choice
# is between biases introduced from listwise, pairwise deletions vs imputation
# See QRISK problems where some key variables were not used for imputation
# What to do with large cohorts which have sub-sampling
# (usually because of financing)?

# Consider library(Hmisc) with aregImpute() using additive regression,
# bootstrapping, and predictive mean matching (pmm)
# It also adapts the method based on variable type automatically
# PMM (with mice or others) for numerical variables
# For categorical in mice use
# polyreg(Bayesian polytomous regression) for factor variables with >= 2 levels
# proportional odds model for ordered variables with >= 2 levels

# Transform and normalise before imputation?
# https://www.theanalysisfactor.com/multiple-imputation-5-recent-findings-that-change-how-to-use-it/
# Don't transform
# don't round off binary vars after imputing
# iterate as many as missing percentage
# Create multiplicative terms before imputing if the model will contain them
# (for eg interaction term or a quadratic), otherwise imputation may be biased (von Hippel, 2009)
# Also probably best not to remove outliers at all
##########

##########
# TO DO:
# Legends and methods output
# Report (box 2 in https://www.bmj.com/content/338/bmj.b2393):
# Amount of missing data:
# number of missing values for each variable
# number of cases with complete data for each important component of the analysis.
# Give reasons for missing values if possible
# indicate how many individuals were excluded because of missing data when
# reporting the flow of participants through the study.
# If possible, describe reasons for missing data in terms of other variables
# (rather than just reporting a universal reason such as treatment failure)

# Comparison of distribution of key variables in individuals
# with and without missing data
# Provide imputed and non-imputed datasets
# Plausibility of the missing at random assumption (discussion and sensitivity analysis)
# For sensitivity analysis example see:
# https://gerkovink.github.io/miceVignettes/Sensitivity_analysis/Sensitivity_analysis.html

# full report of imputation method including:
# assumptions of method
# software used and of key settings
# number of imputed datasets created
# variables were included in the imputation procedure
# treatment of non-normally distributed and binary/categorical variables
# Were statistical interactions also included in imputation models?
# Provide results from analyses restricted to complete cases
# Consider that complete cases analysis may be biased and multiple imputation
# could correct biases (if missing at random).
# investigate the robustness of key inferences to possible departures
# from the missing at random assumption by assuming a range of missing
# not at random mechanisms in sensitivity analyses. See the ampute package for this.
##########
######################

######################
# args:

num_cores <- as.integer(args[['--num_cores']])

# set rownames and exclude from imputation calcs:
ID_col <- as.integer(args[['--ID_col']])
# ID_col <- 1

# pass a few variables of interest for sanity checking:
vars_interest <- as.character(args[['--vars_interest']])
# vars_interest <- c('age', 'hyp')
print(sprintf('Variables of interest are: %s', vars_interest))
# The first var will be the var of interest:
main_var <- vars_interest[[1]]
# main_var <- 'bmi'
print(sprintf('The main variable of interest is: %s', main_var))

# variables such as BMI which should be excluded from the imputation
# procedure if not needed for imputation and recalculated afterwards:
derived_vars <- as.character(args[['--derived_vars']])
# derived_vars <- c('bmi')
print(sprintf('Derived variables to exclude from imputation are: %s', derived_vars))

# Percent of missigness to report (for information only, no action is done):
# Run the per col and per row cleaning script first
missingness_cut <- as.integer(args[['--missingness_cut']])
# missingness_cut <- 30

# for mice():
mincor <- as.numeric(args[['--mincor']]) # minimum correlation for quickpred
# mincor <- 0.3

# Random number seed for reproducible analysis:
seed <- as.integer(args[['--seed']])
# seed <- 123456

# Number of imputed datasets, 5 is default, roughly 1 per percent missing
m <- as.integer(args[['-m']])
# m <- 5

# number of iterations per dataset to impute
maxit <- as.integer(args[['--maxit']])
# maxit <- 30

# Provide a predictor matrix that determines the relationship between variables
# that will be used for imputation:
pred <- as.character(args[['--pred']])
# initial imputation
# to show predictor matrix and methods chosen:
# ini <- mice(input_data, maxit = 0, print = F)
# ini$pred
# to remove a predictor do eg
# pred[ ,"hyp"] <- 0
# object is matrix nested in list

# Provide a methods list to use if not the default:
# TO DO: leave for later and rely on default choices:
# meth <- as.character(args[['--meth']])
# Methods:
# overview of the methods in mice can be found by
# methods(mice)
# ini <- mice(nhanes2, maxit = 0)
# ini$meth
# Change as eg:
# meth["bmi"] <- "norm"
# object is a list
######################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
	this.file = NULL
	# This file may be 'sourced'
	for (i in -(1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}

	if (!is.null(this.file)) return(dirname(this.file))

	# But it may also be called from the command line
	cmd.args = commandArgs(trailingOnly = FALSE)
	cmd.args.trailing = commandArgs(trailingOnly = TRUE)
	cmd.args = cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
	res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

	# If multiple --file arguments are given, R uses the last one
	res = tail(res[res != ""], 1)
	if (0 < length(res)) return(dirname(res))

	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}
Rscripts_dir <- LocationOfThisScript()
print('Location where this script lives:')
Rscripts_dir
# R scripts sourced with source() have to be in the same directory as this one
# (or the path constructed appropriately with file.path) eg:
#source(file.path(Rscripts_dir, 'moveme.R')) #, chdir = TRUE)
######################

######################
# Import libraries
library(mice) # multiple imputation
library(VIM) # visualise missingness
# library(miceadds) # still needed after mice update?
library(parallel) # run imputations with multiple cores
library(data.table)
library(svglite) # prefer over base R svg()
library(lattice) # for density plots

# source functions from a different R script:
source(file.path(Rscripts_dir, 'ggtheme.R')) #, chdir = TRUE)
######################

######################
##########
# Read files, this is with data.table:
if (is.null(args[['-I']]) == FALSE) {
	input_name <- as.character(args[['-I']])
	# For tests:
	# input_name <- 'XXX'
	# setwd('~/xxxx/')
	input_data <- fread(input_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
} else {
	# Stop if arguments not given:
	print('You need to provide an input file. This has to be tab separated with headers.')
	stopifnot(!is.null(args[['-I']]) == TRUE)
}

print('File being used: ')
print(input_name)
##########

##########
# Set output file names:
suffix <- 'imputed'
if (is.null(args[['-O']])) {
	stopifnot(!is.null(args[['-I']]))
	print('Output file name not given. Using:')
	# Split infile name at the last '.':
	input_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
	output_file_name <- sprintf('%s_%s.tsv', input_name, suffix)
	print('Outfile with first complete imputed dataset will be named:')
	print(output_file_name)
	# Save a name also for the long format:
	output_file_name_long <- sprintf('%s_%s_long.tsv', input_name, suffix)
	print('Outfile with all imputed datasets in long format will be named:')
	print(output_file_name_long)
} else {
	output_file_name <- as.character(args[['-O']])
	# output_file_name <- 'testing'
	output_file_name <- sprintf('%s.%s', output_file_name, suffix)
	print(sprintf('Output file names: %s', output_file_name))
}
##########
######################

######################
# Set a seed for reproducible analysis:
set.seed(seed = seed)

# Set-up multiple cores if needed
# Using all cores can slow down the computer, leave one free
# TO DO: pass num_cores
# if num_cores then num_cores = num_cores
# else
num_cores <- max(1, detectCores() - 1)
# Get multiple copies of R running and communicating:
# TO DO: check if FORK is OK for most situations:
cl <- makeCluster(num_cores, type = "FORK") # PSOCK is default but needs env vars passed
# Pass a seed:
clusterSetRNGStream(cl, iseed = seed)
# Use the following if PSOCK is needed:
# Export variables and libraries to so that they are available to all cores:
# clusterExport(cl, input_data) # export all objects needed for function
# clusterEvalQ(cl, library(mice)) # export all libraries needed
# At the end run stopCluster(cl)
# run gc() and rm() if needed # only gc() for garbage collection
######################

######################
# Set method and predictor matrix if provided:
# TO DO:
# if meth given
# meth = meth
# else meth = NULL

# if pred given
# pred = pred
# else
pred <- quickpred(input_data, mincor = mincor)

# if vars_interest
# read in vars_interest
# else use columns 2-5 assuming column 1 has IDs or pass
# and the first var = main_var
# with all being = vars_interest

# if derived_vars
# read in derived_vars
# exclude from data to impute
# include after imputation
# calculate derived var from imputed data
# else pass
# derived_vars = variables such as BMI which should be excluded from the imputation
# procedure but then recalculated afterwards.
######################

######################
# Load the data and inspect
# For testing:
# input_name <- 'nhanes' # use nhanes and nhanes2 as examples
# data(list = input_name)
# input_data <- nhanes

input_data <- as.data.frame(input_data)
# summary(input_data)
# dim(input_data)

# Inspect the missing data pattern
svg('missingness_pattern.svg')
missingness <- md.pattern(input_data, plot = TRUE)
missingness
dev.off()
# TO DO: print out legend
# TO DO: save as table with caption
# each row corresponds to a missing data pattern (1=observed, 0=missing).
# Rows and columns are sorted in increasing amounts of missing information.
# The last column and row contain row and column counts, respectively.
# lower right corner = total missing data points
# last row shows total missing values for each variable
# last column (no header) shows number of missing variables
# first column (no header) x observations with y vars missing

# Check proportion of missing data:
prop_NA <- function(x) {sum(is.na(x)) / length(x) * 100}
# Individuals with more than X% of missing variables:
rows_missing <- apply(input_data, 1, prop_NA) # by rows
rows_above_cut <- nrow(input_data[which(rows_missing > missingness_cut), ])
print(sprintf('Number of rows with >%s%% missing data: %s', missingness_cut, rows_above_cut))
# By columns:
cols_missing <- apply(input_data, 2, prop_NA)
cols_above_cut <- ncol(input_data[, which(cols_missing > missingness_cut)])
print(sprintf('Number of columns with >%s%% missing data: %s', missingness_cut, cols_above_cut))

# See pattern using VIM and mice libraries
# Plot missing values:
svg('missingness_vars_interest_VIM.svg')
# TO DO:
# if vars_interest plot vars_interest, else plot all
# save legend
aggr_plot <- aggr(input_data[, vars_interest],
									only.miss = TRUE, # Plot only missing variables
									col = c('lightgrey', 'red'), # 1 colour for missing data, 2 observed, 3 imputed
									numbers = T,
									sortVars = T,
									labels = names(input_data[, vars_interest]),
									cex.axis = 0.4,
									gap = 2,
									ylab = c('Proportion of missing data', 'Pattern')
									)
dev.off()
######################


######################
# Impute missing data with multiple cores
# Keep samples IDs but exclude them from imputation (non missing and
# would be used to estimate imputation if left as column).
# Check this works as expected with output:
rownames(input_data) <- input_data[[ID_col]]
input_data <- input_data[, -ID_col]
head(input_data)

# Run imputation
# The following will yield num_cores * m imputed datasets
# which will be contained in imp_pars as a list object
# Each list within, eg imp_pars[[1]] will correspond to the structure of
# a mids object, where imp_pars[[1]][1] is data,
# imp_pars[[1]][2] contains the imputed data for each variable, etc.
# mice::ibind merges and attributes it as class mids below
imp_pars <-
	parLapply(cl = cl,
						X = 1:num_cores,
						fun = function(no) {
							mice(input_data,
									 m = m, # Number of imputed datasets, 5 is default
									 maxit = maxit, # max iterations per imputation
									 # quickpred = set the minimum correlation for variable
									 # selection in the predictor matrix:
									 pred = pred,
									 print = F, # omit printing of the iteration cycle
									 diagnostics = TRUE,
									 # meth = 'pmm', # predictive mean matching, leave empty for
				                           # auto selection depending on variable type
				           seed = seed
									 )
							}
						)

# Merge the datasets and create a mids object:
imp_merged <- imp_pars[[1]]
for (n in 2:length(imp_pars)) {
	imp_merged <- mice::ibind(imp_merged,
											imp_pars[[n]])
}

# Free up the cores taken:
stopCluster(cl)
gc(verbose = TRUE) # Prob not necessary but ensure R returns memory to the OS
######################

######################
# Explore the imputed data

##########
# Explore attributes of imputed object which contains the multiply imputed
# data set (class mids) and all information from the procedure including:
# original data, imputed values, number of missing values,
# number of iterations, etc.
attributes(imp_merged)
# imp_merged$data # original data
# imp_merged$imp # imputed data
imp_merged$meth # imputation method used, "" empty string means no NAs
##########

##########
# Explore by visualising the main variables of interest
# Plot vars of interest original data:
# TO DO: save legend
out <- vector(mode = 'list', length = length(vars_interest))
names(out) <- vars_interest
out
for (i in vars_interest) {
	xlab <- sprintf('%s %s, observed values', input_name, i)
	out[[i]] <- densityplot(input_data[[i]],
													xlab = xlab)
}
# Save to disk, one plot per file:
for (i in names(out)) {
	plot_name <- sprintf('densityplots_%s_%s.svg', input_name, i)
	svg(plot_name)
# cols_plot <- max(1, 2 %% length(vars_interest))
# par(mfrow = c(length(vars_interest), cols_plot))
  print(out[[i]])
  dev.off()
}

# Plot all numerical variables with 2 or more missing values:
# densityplot(imp_merged, ~ bmi) # will give only one var, has to be unquoted
# densityplot(imp_merged, ~ bmi | .imp) # will plot each imputed dataset separately
# TO DO: save legend
svg('densityplots_imputed.svg')
lattice::densityplot(imp_merged)
dev.off()
# blue is observed, magenta imputed
##########

##########
# Inspect the convergence of the algorithm
# mice() implements an iterative MCMC type of algorithm.
# Trace lines generated by the algorithm to study convergence:
# Plot convergence of imputed data, only plots the last 3 variables:
svg('convergence_plot_imputations.svg')
plot(imp_merged)
dev.off()
# TO DO: save legend
# The plot shows the mean (left) and standard deviation (right)
# of the imputed values only.
# In general, we would like the streams to intermingle
# and be free of any trends at the later iterations.
##########

##########
# Check the imputed data vs the variables of interest:
# TO DO: require unquoted variables to be passed
# TO DO: save legend
# svg('missing_data_scatterplots_vars_interest.svg')
# xyplot(imp_merged,
# 			 bmi ~ age,
# 			 # pch = 1, cex = 1, strip = T,
# 			 ylab = sprintf('%s as predicted by %s', main_var, vars_interest[-1]),
# 			 xlab = '',
# 			 strip = strip.custom(factor.levels = labels),
# 			 type = c('p')) # Magenta are imputed, blue observed
# dev.off()


# Further exploratory plots
# # TO DO: save legends
# Diagnostics for plausible values. compare imputed vs observed values
# Assuming data are missing completely at random (MCAR)
# imputations should have the same distribution as the observed data.
# Distributions may differ because missing data are
# missing at random (MAR) or non-random (MNAR)
# Very large discrepancies should not exist though, check with:
svg('bwplots_imputations.svg')
bwplot(imp_merged)
dev.off()

# Stripplots might look better, check the first 10 imputed datasets:
svg('stripplots_imputations.svg')
stripplot(imp_merged,
					subset = (.imp == 1 | .imp == 2 | .imp == 3 | .imp == 4 | .imp == 5 |
											.imp == 6 | .imp == 7 | .imp == 8 | .imp == 9 | .imp == 10),
					# col = mdc(1:2), #col = mdc(1:2), pch=20, cex=1.5,
					pch = 1, cex = 0.7,
					strip = strip.custom(par.strip.text = list(cex = 0.7)))
# Magenta are imputed, blue observed
dev.off()

# Can also run for variables of interest only:
# TO DO: needs unquoted vars
# stripplot(imp_merged, vars_interest~.imp, pch = 20, cex = 2)

# TO DO: save legend
# Legend: Strip plot of observed (blue) and imputed (red) values for each variable imputed.
# The figure shows whether distributions are similar.
# If using the PMM method, imputed values have the same gaps as the observed data,
# and are always within the range of the observed data.

# Under MCAR, univariate distributions of the observed and imputed data
# are expected to be identical.
# Under MAR, they can be different, both in location and spread,
# but their multivariate distribution is assumed to be identical.

# The only way to test the amount of missing data that is to test
# imputed vs complete cases and see if the results change
##########

######################
# TO DO: needed?
# Extend the number of iterations
# five iterations works well in practice
# can extend with mice.mids():
# imp_extended <- mice.mids(imp_merged, maxit = 35, print = F)
# plot(imp_extended)
######################

######################
##########
# Extract the completed data:
imp_merged_comp <- complete(imp_merged, 1) # first imputed data set
# complete() provides more options for exploring with include and action
# for outputting both observed and imputed and long, broad, etc. formatted datasets
# Sanity check the number of missing values, will error if complete though:
# md.pattern(imp_merged_comp)

# TO DO: save all datasets for regression with pooled data in eg --session.
##########

##########
# TO DO: this would be manual:
# # Sanity check observed and imputed data
# # Variables used as predictors for imputation of each incomplete variable:
# imp_merged$pred
# # Implausible results for specific variables:
# which(imp_merged$imp$vitd12 <= 1 | imp_merged$imp$vitd12 >= 250)
# which(imp_merged$imp$calendar_age_ra <= 60 | imp_merged$imp$calendar_age_ra >= 100)
##########
######################


######################
# Save files
# First imputed dataset in broad format:
fwrite(imp_merged_comp, output_file_name,
			 sep = '\t', na = 'NA',
			 col.names = TRUE, row.names = FALSE,
			 quote = FALSE,
			 nThread = num_cores # use all cores available to write out
)
# All imputed datasets in long format:
fwrite(imp_merged, output_file_name_long,
			 sep = '\t', na = 'NA',
			 col.names = TRUE, row.names = FALSE,
			 quote = FALSE,
			 nThread = num_cores # use all cores available to write out
)
######################

######################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
objects_to_save <- c('imp_merged') # needs to be a character vector for save()

# Filename to save current R session, data and objects at the end:
if (is.null(args[['--session']]) == FALSE) {
	save_session <- as.character(args[['--session']]) #args $ `--session`
	R_session_saved_image <- sprintf('%s.RData', save_session)
	print(sprintf('Saving an R session image as: %s', R_session_saved_image))
	save(list = objects_to_save, file = R_session_saved_image, compress = 'gzip')
	# TO DO: check time to compress and if gzip is equivalent and faster:
	save(list = objects_to_save, file = 'imputation_test.RData')
} else {
	print('Not saving an R session image, this is the default. Specify the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf
# by default, it doesn't look like there is an easy way to suppress it, so deleting here:
print('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
print('Finished successfully.')
sessionInfo()
q()

# Next: run the script for xxx
######################