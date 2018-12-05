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

Usage: run_mice_impute.R (-I <INPUT_FILE>) (--vars-interest <FILE_VARS>) [options]
       run_mice_impute.R (-I <INPUT_FILE>) (--extend) (--maxit <int>) [options]
       run_mice_impute.R (-I <INPUT_FILE>) (--dry-run) [options]
       run_mice_impute.R [-h | --help]

Options:
-I <INPUT_FILE>                 Input file name
-O <OUTPUT_PREFIX>              Output file prefix, _imputed.tsv or imputed_long.tsv will be added
--session                       R session if to be saved, only saves the mids imputed object
-h --help                       Show this screen
--ID-col <int>                  Column with ID values for rows [default: 1]
--vars-interest <FILE_VARS>     Variables of interest for sanity checking
--num-cores <int>               Number of cores to use for parallelising
--derived-vars <FILE_EXCLUDE>   Derived variables to exclude from imputation
--missingness-cut <int>         Percent of missigness to exclude [default: 30]
--mincor <num>                  Minimum correlation for quickpred [default: 0.3]
--seed <int>                    Random seed number for reproducible analysis [default: 123456]
-m <int>                        Number of imputed datasets to generate [default: 5]
--maxit <int>                   Number of iterations per dataset to impute [default: 30]
--pred <FILENAME>               Predictor matrix to determine relationship between variables
--meth <FILENAME>               Methods to use for imputation for each variable
--dry-run                       Run a dry imputation (maxit = 0) to get the predictor matrix and methods
--extend                        Run further iterations (not implemented, use session instead)


Input
=====

Input file: A tab separated file with headers. This is usually a phenotype file
with individuals (samples) as rows and phenotype information (age, gender, etc.)
as columns. The file is read with data.table and stringsAsFactors = FALSE
A per column and per row cleaned data frame is expected with the first column
as the IDs.

Variables of interest: A tab separated file with no header containing one column with
variable names per row is required. Variable names must match the column names in the
input file (and have NAs). The first variable will be used as the main variable of
interest. Several plots and diagnostic tests are performed using these.

Derived variables: (optional) a file in the same format as vars-interest with
variable names that are derived, not needed for imputation and which should be excluded
from this step (an re-calculated manually afterwards). See below for more information.

The missingness-cut option is used for information purposes only, no action is taken. You should
clean per row and per column first. The value passed is used for both rows and columns.

The dry-run option will save the methods and predictors from mice. You need to provide the input
file and optionally the derived_vars file if excluding any variables.

Use the extend option if you need to increase the number of iterations to check
if convergence is reached. Pass the long format (as output by mice::complete())
as input plus an increased maxit.

Output
======

A tab separated file with headers with the first imputed dataset.
A tab separated file with headers in the long format with all imputed datasets.
This can then be read back into R and converted to a mids object with as.mids().
With session (often slow), an Rdata object that contains all the imputed
(mids object) datasets, this can then be re-loaded into R and ran in pooled
analysis. The session only saves the imputed object.
Several diagnostic plots (convergence, strip plots, bwplots, density plots).

Notes
=====

# Parallelising

The script tries to use multiple cores if available. Using all cores can slow
down the computer considerably so the number of cores is passed as
max(1, detectCores() - 1). If num-cores is provided there is no checking.

num-cores = 1 will still use the parallel setup code so expect it to be
slower than code without. If you want to use everything you have simply pass
num-cores as the number of cores available.

The cluster type is set to FORK, on Windows this probably will not run.

Parallelising will create num-cores * m imputed datasets.

# Changing the methods and predictor matrix used for imputation

You can change the predictor matrix which specifies the variables that are used to
impute each incomplete variable with pred. Create a dry run first (remember to
also provide the derive-vars file if excluding any variables):

  > Rscript run_mice_impute.R #-I my_data.tsv #--dry-run
    (ignore #)

Modify the predictor matrix and methods as required, then run the imputation with:

  > Rscript run_mice_impute.R #-I my_data.tsv #--pred predictor_matrix_my_data.tsv #--meth methods_my_data.tsv
    (ignore #)

In the output, a value of 1 indicates that the column variable was used to
impute the row variable. A row with all zeros indicates no missing values in that
variable. You can change it to suit the specific relationships you need.

By default, this script uses the option quickpred() for a quick selection
procedure of predictors based on variable correlations of at least p = 0.30
(set using the option mincor).

# Variables of interest for sanity checking and diagnostics

You need to supply a variable of interest and some covariates to the
vars-interest option as a single column file without header.
These need to have NAs to be used. The first variable will be taken as the
variable of interest.

# Imputing derived variables

If derived variables are not needed for imputation, exclude them from imputation,
impute the originals, and then calculate the derived variables.
This is the approach used here, known as impute, then transform (Von Hippel 2009).

Provide the derived-vars option as a single column file without header for this.
These variables will be excluded from imputation. You will then need to
recalculate them for the missing data. eg if BMI is present and not needed
for imputation (it provides no additional information useful to impute other
variables), pass it in the derived-vars option, this will be excluded for the
imputation procedure. Once this is done, use the output from this script
and run eg:

(this is from https://stefvanbuuren.name/fimd/sec-knowledge.html):
data <- boys[, c("age", "hgt", "wgt", "hc", "reg")]
imp <- mice(data, print = FALSE, seed = 71712)
long <- mice::complete(imp, "long", include = TRUE)
long$whr <- with(long, 100 * wgt / hgt)
imp.itt <- as.mids(long)

The first three steps are carried out here, the last need to be done manually
according to your needs.

You could also consider passive imputation for this type of variables,
as explained in:
https://gerkovink.github.io/miceVignettes/Passive_Post_processing/Passive_imputation_post_processing.html

For an overview see:
https://stefvanbuuren.name/fimd/sec-knowledge.html

If you are planning to carry out a model which includes interactions,
non-linear effects of variables or Cox proportional hazards model (which is
non-linear), this script and the MICE approach used will not be appropriate
as the imputation model and the subsequent analysis model will not be compatible.
In this case you may need to use the SMC-FCS approach (as suggested elsewhere).
See:
https://cran.r-project.org/web/packages/smcfcs/vignettes/smcfcs-vignette.html

' -> doc
# Load docopt:
library(docopt, quietly = TRUE)
# Retrieve the command-line arguments:
args <- docopt(doc)
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
# Derivation and validation of QRISK, a new cardiovascular disease risk score for
# the United Kingdom: prospective open cohort study | The BMJ
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
# Multiple Imputation for General Missing Data Patterns in the Presence of
# High-dimensional Data | Scientific Reports
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

# An example of how to create a predictor matrix and modify it:
# https://rpubs.com/kaz_yos/mice-exclude

# TO DO: could add more visualisations eg:
# https://datascienceplus.com/graphical-presentation-of-missing-data-vim-package/
# https://cran.r-project.org/web/packages/VIMGUI/vignettes/VIM-Imputation.pdf

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
# (for eg interaction term or a quadratic), otherwise imputation
# may be biased (von Hippel, 2009)
# Also probably best not to remove outliers at all
##########
######################

######################
# This function allows other R scripts to obtain the path to a script directory
# (ie where this script lives). Useful when using source('some_script.R')
# without having to pre-specify the location of where that script is.
# This is taken directly from:
# How to source another_file.R from within your R script molgenis/molgenis-pipelines Wiki
# https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
# Couldn't find a licence at the time (12 June 2018)
LocationOfThisScript = function() # Function LocationOfThisScript returns the
	                                # location of this .R script
	                                # (may be needed to source other files in same dir)
{
	this.file = NULL
	# This file may be 'sourced'
	for (i in -(1:sys.nframe())) {
		if (identical(sys.function(i),
									base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}

	if (!is.null(this.file)) return(dirname(this.file))

	# But it may also be called from the command line
	cmd.args = commandArgs(trailingOnly = FALSE)
	cmd.args.trailing = commandArgs(trailingOnly = TRUE)
	cmd.args = cmd.args[seq.int(from = 1,
															length.out = length(cmd.args) - length(cmd.args.trailing))]
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
library(parallel) # run imputations with multiple cores,
                  # this is core package only needs loading, not installing
library(data.table) # data manipulation and fast reading and writing
library(svglite) # prefer over base R's svg()
library(lattice) # for density plots

# source functions from a different R script:
# source(file.path(Rscripts_dir, 'ggtheme.R')) #, chdir = TRUE)
######################

######################
# Get arguments

##########
# minimum correlation for quickpred:
mincor <- as.numeric(args[['--mincor']])
# TO DO: add minpuc for minimum proportion of usable cases, default is 0
# https://www.rdocumentation.org/packages/mice/versions/3.3.0/topics/quickpred
# TO DO: mice can be very slow with large datasets (eg Airwave ~35000 rows, ~230 columns
# with <30% missing rows and columns, takes >3 days runnining in parallel
# profile script
# TO DO: check logged events for problems, eg:
# ini <- mice(data, maxit = 0)   # recommended
# head(ini$loggedEvents, 2)
# TO DO: check/plot the global influx-outflux pattern with mice::flux(), see
# https://stefvanbuuren.name/fimd/missing-data-pattern.html#sec:flux
# see below for fluxplot()

num_cores <- as.integer(args[['--num-cores']])

# set rownames and exclude from imputation calcs:
ID_col <- as.integer(args[['--ID-col']])
# ID_col <- 1

# Percent of missigness to report (for information only, no action is done):
# Run the per col and per row cleaning script first
missingness_cut <- as.integer(args[['--missingness-cut']])
# missingness_cut <- 30

# for mice():
# Random number seed for reproducible analysis:
seed <- as.integer(args[['--seed']])
# seed <- 123456

# Number of imputed datasets, 5 is default, roughly 1 per percent missing
m <- as.integer(args[['-m']])
# m <- 5

# number of iterations per dataset to impute
maxit <- as.integer(args[['--maxit']])
# maxit <- 30
##########

##########
# Read files, this is with data.table:
if (!is.null(args[['-I']])) {
	input_name <- as.character(args[['-I']])
	# input_name <- 'nhanes.tsv'
	input_data <- fread(input_name,
											sep = '\t',
											header = TRUE,
											stringsAsFactors = FALSE)
} else {
	# Stop if arguments not given:
	print('You need to provide an input file. This has to be tab separated with headers.')
	stopifnot(!is.null(args[['-I']]))
}

print('File being used: ')
print(input_name)
##########

##########
# Keep sample IDs but exclude them from imputation (non missing and
# would be used to estimate imputation if left as column).
# Only run if -I given but without --extend
if (!is.null(args[['-I']]) &  # arg is NULL
		args[['--extend']] == FALSE) {  # arg is boolean
	# When writing out, the long format has .id, as character,
	# corresponding to the row names
	# row names are also saved to the dataframes with fwrite at the end.
	input_data <- as.data.frame(input_data) # needed for plot functions and rownames
	rownames(input_data) <- input_data[[ID_col]]
	input_data <- input_data[, -c(ID_col)]
	print('First few rows of data:')
	head(input_data)
}
##########

##########
# Set output file names:
suffix <- 'imputed'
if (is.null(args[['-O']])) {
	stopifnot(!is.null(args[['-I']]))
	print('Output file name not given.')
	# Split infile name at the last '.':
	output_name <- strsplit(input_name, "[.]\\s*(?=[^.]+$)", perl = TRUE)[[1]][1]
	output_file_name <- sprintf('%s_%s.tsv', output_name, suffix)
	# Save a name also for the long format:
	output_file_name_long <- sprintf('%s_%s_long.tsv', output_name, suffix)
} else {
	output_name <- as.character(args[['-O']])
	output_file_name <- sprintf('%s_%s.tsv', output_name, suffix)
	# Save a name also for the long format:
	output_file_name_long <- sprintf('%s_%s_long.tsv', output_name, suffix)
}
print('Outfile with first complete imputed dataset will be named:')
print(output_file_name)
print('Outfile with all imputed datasets in long format will be named:')
print(output_file_name_long)
##########

##########
# Pass a few variables of interest for sanity checking:
if (!is.null(args[['--vars-interest']]) &
		args[['--dry-run']] == FALSE) { # dry-run doesn't need vars_interest
	vars_interest <- as.character(args[['--vars-interest']])
	# vars_interest <- 'vars_interest.tsv'
	vars_interest <- fread(vars_interest,
												 sep = '\t',
												 header = FALSE,
												 stringsAsFactors = FALSE)
	vars_interest <- as.character(vars_interest[, V1])
	print('Variables of interest:')
	print(vars_interest)
	# The first var will be the var of interest:
	main_var <- vars_interest[[1]]
	# main_var <- 'bmi'
	print(sprintf('The main variable of interest is: %s', main_var))
	# TO DO:
	# also calculate intraclass correlation (ICC) to identify cluster structure
	# in our dataset, as in:
	# https://gerkovink.github.io/miceVignettes/Multi_level/Multi_level_data.html
	# library(pan); library(multilevel)
	# step 7 for the vars_interest (which must have NAs), as eg
	# icc(aov(main_var ~ input_data[, vars_interest[[2]]], data = input_data))
	# and for each pairwise correlation
	# If there is multilevel structure present, account for it in imputation
	# use --dry-run, modify pred and meth using a fixed effects approach and impute
	# Compare ICC for observed, imputed 1, and imputed fixed effects for vars_interest
}
# else {
# 	# Stop if arguments not given:
# 	print('You need to provide a file with some variables of interest to test.')
# 	stopifnot(!is.null(args[['--vars-interest']]) &
# 						args[['--dry-run']] == FALSE)
# }
##########

##########
# variables such as BMI which should be excluded from the imputation
# procedure if not needed for imputation and recalculated afterwards:
if (!is.null(args[['--derived-vars']])) {
	derived_vars <- as.character(args[['--derived-vars']])
	# derived_vars <- 'derived_vars.tsv'
	derived_vars <- fread(derived_vars,
												 sep = '\t',
												 header = FALSE,
												 stringsAsFactors = FALSE)
	derived_vars <- as.character(derived_vars[, V1])
	print(sprintf('Derived variables to exclude from imputation: %s', derived_vars))
	# Exclude from dataset provided:
	input_data <- input_data[, !(names(input_data) %in% derived_vars)]
	print('First few rows of data after excluding derived vars:')
	head(input_data)
} else {
	print('Derived variables not provided, using all columns for imputation.')
}
##########

##########
# Create a run dry in order to print out correlation matrix and methods chosen:
if (args[['--dry-run']] == TRUE) { # arg is boolean
	print('Running a dry imputation to get methods and predictor matrix.')
	dry_mice <- mice(input_data, maxit = 0, print = F)
	# Save predictor matrix:
	# pred[ ,"hyp"] <- 0
	fwrite(as.data.frame(dry_mice$pred),
				 sprintf('predictor_matrix_%s', output_file_name),
				 sep = '\t',
				 na = 'NA',
				 col.names = TRUE,
				 row.names = TRUE,
				 quote = FALSE
				 )
	# Save methods:
	# overview of the methods in mice can be found by
	# methods(mice)
	# dry_mice$meth
	# Change as eg:
	# meth["bmi"] <- "norm"
	fwrite(as.list(dry_mice$meth),
				 sprintf('methods_%s', output_file_name),
				 sep = '\t',
				 na = 'NA',
				 col.names = TRUE,
				 row.names = FALSE
	)
	# And quit after this:
	print('Dry run finished, exiting.')
	sessionInfo()
	q()
	}
##########

##########
# TO DO:
# Provide a post imputation option, see step 4 in:
# https://gerkovink.github.io/miceVignettes/Passive_Post_processing/Passive_imputation_post_processing.html
# where eg:
# do a dry-run, then provide --meth --pred and --post, then run actual imputation
# post allows processing of particular variables to constrain them (eg limit them
# to positive values, given range, etc.)
##########

##########
# Set method and predictor matrix if provided
# Provide a predictor matrix that determines the relationship between variables
# that will be used for imputation:
if (!is.null(args[['--pred']])) { # arg is NULL
	pred_name <- as.character(args[['--pred']])
	# pred_name <- 'output_dry_run/predictor_matrix_nhanes_imputed.tsv'
	pred <- as.data.frame(fread(pred_name,
													    sep = '\t',
													    header = TRUE,
													    stringsAsFactors = FALSE,
													)
										)
	# Needs to be passed as a matrix, first column are variable names, move to rownames:
  rownames(pred) <- pred[[1]]
  pred <- pred[, -1]
	pred <- as.matrix(pred)
  print(sprintf('Predictor matrix provided, using this instead of mice pred: %s',
								pred_name)
				)
	print('First few rows of predictor matrix provided: ')
	head(pred)
} else {
	pred <- quickpred(input_data, mincor = mincor)
	print('Predictor matrix not provided, ')
  print('using mice defaults with quickpred and ')
  print(sprintf('%s for minimum correlation between variables.', mincor))
	}

# Provide a methods list to use if not the default:
if (!is.null(args[['--meth']])) { # arg is NULL
	meth_name <- as.character(args[['--meth']])
	# meth_name <- 'methods_nhanes_imputed.tsv'
	print(sprintf('Methods provided, using this instead of mice defaults: %s',
								meth_name)
	)
	meth <- fread(meth_name,
								sep = '\t',
								header = TRUE,
								stringsAsFactors = FALSE
								# na.strings = ""
								)
	meth_names <- names(meth)
	meth <- as.character(meth[1, ])
	names(meth) <- meth_names
	# Convert NA to literal empty strings:
	meth <- lapply(meth, function(x){replace(x, x == "NA", as.character(""))})
	print('First few elements of methods provided: ')
	head(meth)
} else {
	meth <- NULL
	print('Methods not provided, using mice defaults')
	}
##########

##########
# Set a seed for reproducible analysis:
set.seed(seed = seed)

# Set-up multiple cores if needed
if (!is.null(args[['--num-cores']])) { # arg is NULL
	num_cores <- as.integer(args[['--num-cores']])
	print(sprintf('Number of cores provided, using: %s', num_cores))
} else {
	# Using all cores can slow down the computer, leave one free:
	num_cores <- max(1, detectCores() - 1)
	print(sprintf('Detected cores, using: %s', num_cores))
}

# Setup the cluster
# FORK runs only in Unix like, PSOCK is default but needs env vars passed to each core
cl <- makeCluster(num_cores, type = "FORK")
# Pass a seed:
clusterSetRNGStream(cl, iseed = seed)
# Use the following if PSOCK is needed:
# Export variables and libraries to so that they are available to all cores:
# clusterExport(cl, input_data) # export all objects needed for function
# clusterEvalQ(cl, library(mice)) # export all libraries needed
# At the end run stopCluster(cl)
# run gc() and rm() if needed # only gc() for garbage collection
##########
######################

######################
# Inspect data

# Only run if -I given but without --extend:
if (!is.null(args[['-I']]) &  # arg is NULL
		args[['--extend']] == FALSE) {  # arg is boolean
	# summary(input_data)
	# dim(input_data)

	# Inspect the missing data pattern:
	# TO DO: print out legend
	svg(sprintf('missingness_pattern_%s.svg', output_name))
	missingness <- md.pattern(input_data, plot = TRUE)
	dev.off()

	# Save missingness pattern:
	fwrite(as.data.frame(missingness),
				 sprintf('missingness_pattern_%s.tsv', output_name),
				 sep = '\t',
				 na = 'NA',
				 col.names = TRUE,
				 row.names = TRUE, # keep as these are the ID column
				 quote = FALSE
	)
	# TO DO: save as table with caption
	# each row corresponds to a missing data pattern (1=observed, 0=missing).
	# Rows and columns are sorted in increasing amounts of missing information.
	# number of rows is equal to the number of patterns identified
	# first column (no header) x observations with y vars missing
	# The last column and row contain row and column counts, respectively.
	# lower right corner = total missing data points
	# last row shows total missing values for each variable
	# last column (no header) shows number of missing variables


	# TO DO:
	# Does the missing data of var x depend on var y?
	# Plot histograms conditional on missingness for vars_interest eg:
	# https://gerkovink.github.io/miceVignettes/Missingness_inspection/Missingness_inspection.html
	# (in step 8)
	# R <- is.na(boys$gen)
	# histogram(~age|R, data=boys)

	# TO DO:
	# Add a fluxplot to identify powerful predictors
	# https://gerkovink.github.io/miceVignettes/Sensitivity_analysis/Sensitivity_analysis.html
	# fx <- fluxplot(input_data)
	#

	# Check proportion of missing data:
	prop_NA <- function(x) {sum(is.na(x)) / length(x) * 100}
	# Individuals with more than X% of missing variables:
	rows_missing <- apply(input_data, 1, prop_NA) # by rows
	rows_above_cut <- nrow(input_data[which(rows_missing > missingness_cut), ])
	print(sprintf('Number of rows with >%s%% missing data: %s',
								missingness_cut, rows_above_cut))
	# By columns:
	cols_missing <- apply(input_data, 2, prop_NA)
	cols_above_cut <- ncol(input_data[, which(cols_missing > missingness_cut)])
	print(sprintf('Number of columns with >%s%% missing data: %s',
								missingness_cut, cols_above_cut))

	# See pattern using VIM and mice libraries
	# Plot aggregations for missing/imputed values:
	# TO DO: could move this to post imputation to make full use
	# input_data gets rewritten below though, and is used after input from both
	# long format for --extension and imputation
	svg(sprintf('missingness_vars_interest_VIM_%s.svg', output_name))
	# TO DO: save legend
	aggr_plot <- aggr(input_data,#[, vars_interest],
										combined = FALSE, # plot bar and pattern separately
										only.miss = FALSE, # Plot combinations only for missing variables
										numbers = TRUE,
										sortVars = TRUE,
										labels = names(input_data),#[, vars_interest]),
										cex.axis = 0.4,
										gap = 2,
										ylab = c('Proportion of missing data', 'Pattern')
	)
	dev.off()
	summary(aggr_plot)
	# Barplot (left) shows the proportion of missing or imputed values in each variable.
	# Aggregation plot (middle) shows all existing combinations of of  missing  (red),
	# imputed (orange) and observed (blue) values.
	# Barplot (right) shows the frequencies of different variable combinations
	}
######################

######################
# Extend the number of iterations with option --extend
# Get the appropriate long file and convert with as.mids(), then impute
# Run if both arguments are given:
if (!is.null(args[['-I']]) &  # arg is NULL
		args[['--extend']] == TRUE) { # arg is boolean
	input_name <- as.character(args[['-I']])
	# input_name <- 'nhanes_imputed_long.tsv'
	input_data <- fread(input_name,
											sep = '\t',
											header = TRUE,
											stringsAsFactors = FALSE)
	print('File being used to extend imputations: ')
	print(input_name)
	input_data <- as.mids(as.data.frame(input_data))
	# head(input_data)
	# class(input_data)
} else { # if not extend
  print('--extend not in options.')
}
######################


######################
# Impute missing data with multiple cores

##########
# Run imputation
# The following will yield num_cores * m imputed datasets
# which will be contained in imp_pars as a list object
# Each list within, eg imp_pars[[1]] will correspond to the structure of
# a mids object, where imp_pars[[1]][1] is data,
# imp_pars[[1]][2] contains the imputed data for each variable, etc.
# mice::ibind merges and attributes it as class mids below

# TO DO: check adding extension works OK, when parallelising and with ibind()
# Only run if -I given but without --extend
if (!is.null(args[['-I']]) &  # arg is NULL
		args[['--extend']] == FALSE) {  # arg is boolean
	print('Starting imputations.')
	print(sprintf('Total number of imputed datasets to complete: %s', num_cores * m))
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
										 meth = meth,
					           seed = seed
										 )
								}
							)
} else if (!is.null(args[['-I']]) & # if both arguments are given run
					 args[['--extend']] == TRUE) {
	# imp_pars <- mice.mids(input_data, maxit = 35, print = F)
  print('Extending iterations.')
	imp_pars <-
		parLapply(cl = cl,
							X = 1:num_cores,
							fun = function(no) {
								mice.mids(input_data,
													maxit = maxit, # max iterations per imputation
													print = F # omit printing of the iteration cycle
								)
								}
							)
	# plot(imp_pars)
	}

# Merge the datasets and create a mids object:
imp_merged <- imp_pars[[1]]
for (n in 2:length(imp_pars)) {
	imp_merged <- mice::ibind(imp_merged,
											imp_pars[[n]])
}
##########

##########
# Free up the cores taken:
stopCluster(cl)
gc(verbose = TRUE) # Prob not necessary but ensure R returns memory to the OS
# TO DO: convert all to functions...
##########
######################

######################
# Explore the imputed data

##########
# Explore attributes of imputed object which contains the multiply imputed
# data set (class mids) and all information from the procedure including:
# original data, imputed values, number of missing values,
# number of iterations, etc.
# attributes(imp_merged)

# Save predictor matrix:
fwrite(as.data.frame(imp_merged$pred),
			 sprintf('predictor_matrix_%s', output_file_name),
			 sep = '\t',
			 na = 'NA',
			 col.names = TRUE,
			 row.names = TRUE,
			 quote = FALSE
)
# Save methods:
# imputation method used, "" empty string means no NAs
fwrite(as.list(imp_merged$meth),
			 sprintf('methods_%s', output_file_name),
			 sep = '\t',
			 na = 'NA',
			 col.names = TRUE,
			 row.names = FALSE
)

# TO DO: save to file as methods description:
# Report (box 2 in https://www.bmj.com/content/338/bmj.b2393):
# Amount of missing data:
# number of missing values for each variable
# number of cases with complete data for each important component of the analysis.
# Give reasons for missing values if possible
# indicate how many individuals were excluded because of missing data when
# reporting the flow of participants through the study.
# If possible, describe reasons for missing data in terms of other variables
# (rather than just reporting a universal reason such as treatment failure)
# Use extended STROBE guidelines for reporting of multiple imputation analyses
# https://www.strobe-statement.org/index.php?id=available-checklists

# Comparison of distribution of key variables in individuals
# with and without missing data
# Provide imputed and non-imputed datasets
# Plausibility of the missing at random assumption (discussion and sensitivity analysis)
# For sensitivity analysis example see:
# https://gerkovink.github.io/miceVignettes/Sensitivity_analysis/Sensitivity_analysis.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5358992/pdf/clep-9-157.pdf
# Also see how to evaluate the effect of missing data on statistical
# inferences with mice::ampute:
# https://rianneschouten.github.io/mice_ampute/vignette/ampute.html

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
# Also:
# After adding back derived variables, run scatter plots of eg BMI observed vs
# BMI calculated post imputation; plus any other relationships that should match,
# both to check concordance and check imputation values are plausible
##########

##########
# Explore by visualising the main variables of interest
# Plot vars of interest original data:
# TO DO: save legend
out <- vector(mode = 'list', length = length(vars_interest))
names(out) <- vars_interest
for (i in vars_interest) {
	xlab <- sprintf('%s %s, observed values', input_name, i)
	out[[i]] <- densityplot(input_data[[i]],
													xlab = xlab)
}
# Save to disk, one plot per file:
for (i in names(out)) {
	plot_name <- sprintf('densityplots_%s_%s.svg', output_name, i)
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
svg(sprintf('densityplots_imputation_%s.svg', output_name))
lattice::densityplot(imp_merged)
dev.off()
# blue is observed, magenta imputed
##########

##########
# Inspect the convergence of the algorithm
# mice() implements an iterative MCMC type of algorithm.
# Trace lines generated by the algorithm to study convergence:
# Plot convergence of imputed data, only plots the last 3 variables:
svg(sprintf('convergence_plots_imputation_%s.svg', output_name))
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
# TO DO: pass as unquoted variables
# TO DO: save legend
# svg(sprintf('missing_data_scatterplots_vars_interest_%s.svg', output_name))
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
svg(sprintf('bwplots_imputation_%s.svg', output_name))
bwplot(imp_merged,
			 subset = (.imp == 0 | # get the original data
			 					 .imp == 1 | .imp == 2 | .imp == 3 | .imp == 4 | .imp == 5 |
			 					 .imp == 6 | .imp == 7 | .imp == 8 | .imp == 9 | .imp == 10),
			 # col = mdc(1:2), #col = mdc(1:2), pch=20, cex=1.5,
			 pch = 1, cex = 0.7,
			 strip = strip.custom(par.strip.text = list(cex = 0.7))
			 )
dev.off()

# Stripplots might look better, check the first 10 imputed datasets:
svg(sprintf('stripplots_imputation_%s.svg', output_name))
stripplot(imp_merged,
					subset = (.imp == 0 | # get the original data
										.imp == 1 | .imp == 2 | .imp == 3 | .imp == 4 | .imp == 5 |
										.imp == 6 | .imp == 7 | .imp == 8 | .imp == 9 | .imp == 10),
					# col = mdc(1:2), #col = mdc(1:2), pch=20, cex=1.5,
					pch = 1, cex = 0.7,
					strip = strip.custom(par.strip.text = list(cex = 0.7))
					)
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

######################
# TO DO: this would be manual:
# # Sanity check observed and imputed data
# # Variables used as predictors for imputation of each incomplete variable:
# imp_merged$pred
# # Implausible results for specific variables:
# which(imp_merged$imp$vitd12 <= 1 | imp_merged$imp$vitd12 >= 250)
# which(imp_merged$imp$calendar_age_ra <= 60 | imp_merged$imp$calendar_age_ra >= 100)
######################

######################
# Extract the completed data:
imp_merged_comp <- complete(imp_merged, 1) # first imputed data set
# complete() provides more options for exploring with include and action
# for outputting both observed and imputed and long, broad, etc. formatted datasets
# Sanity check the number of missing values, will error if complete though:
# md.pattern(imp_merged_comp)

# Save the long format if needed for extensions, rownames are taken as .id:
imp_merged_long <- complete(imp_merged, action = 'long', include = TRUE)
# inlcude = T, saves also the original data, as is suggested
# This format can then be passed back to increase the number of imputations with:
# as.mids(imp_merged_long)
######################

######################
# Save files
# First imputed dataset in broad format:
fwrite(imp_merged_comp, output_file_name,
			 sep = '\t',
			 na = 'NA',
			 col.names = TRUE,
			 row.names = TRUE, # keep as these are the ID column
			 quote = FALSE,
			 nThread = num_cores # use all cores available to write out
)
# All imputed datasets in long format:
fwrite(imp_merged_long, output_file_name_long,
			 sep = '\t',
			 na = 'NA',
			 col.names = TRUE,
			 row.names = FALSE, # don't keep, IDs are in .id
			 quote = FALSE,
			 nThread = num_cores # use all cores available to write out
)

gc()
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
if (!is.null(args[['--session']])) { # arg is NULL
	save_session <- sprintf('%s_%s.RData', output_name, suffix)
	print(sprintf('Saving an R session image as: %s', save_session))
	save(list = objects_to_save, file = save_session, compress = 'gzip')
} else {
	print('Not saving an R session image, this is the default. Specify
         the --session option otherwise')
}

# If using Rscript and creating plots, Rscript will create the file Rplots.pdf
# by default, it doesn't look like there is an easy way to suppress it,
# so deleting here:
print('Deleting the file Rplots.pdf...')
system('rm -f Rplots.pdf')
print('Finished successfully.')
sessionInfo()
q()

# Next: run pooled analyses
######################
