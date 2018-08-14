#!/usr/bin/env bash

###########################
# Remove non-printable characters from files
# July 2018
# Antonio
###########################

###########################
# Set bash script options

# exit when a command fails:
set -o errexit

# exit if any pipe commands fail:
set -o pipefail

# exit if there are undeclared variables:
set -o nounset

# trace what gets executed:
set -o xtrace
set -o errtrace
###########################

###########################
# Variables to substitute:
INFILE=$1
OUTFILE=$2
###########################

###########################
# Files have non-printable characters:
# <AE> is £ in 'Auto_025_290216.csv'
# ? is (R) copyright symbol in 'trt_025_290216.csv'
# Change LANG as my default is UTF-8 and these are latin1?
LC_ALL=C LC_TYPE=C LANG=C sed $'s/[^[:print:]\t]//g' ${INFILE} > ${OUTFILE}
###########################