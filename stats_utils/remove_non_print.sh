#!/usr/bin/env bash

###########################
# Remove non-printable characters from files
# July 2018
# Antonio
###########################

###########################
# Usage eg:
#remove_non_print.sh infile outfile

# Several more notes below for identifying non-ascii characters and other commands that may be helpful.
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
# Remove all special characters in Linux text:
#https://stackoverflow.com/questions/43108359/how-to-remove-all-special-characters-in-linux-text

# Also: change LANG (if eg system is UTF-8) to avoid error: "error: illegal byte sequence"
#https://stackoverflow.com/questions/19242275/re-error-illegal-byte-sequence-on-mac-os-x#19770395

# eg:
# <AE> is £
# ? is (R) copyright symbol
# Change LANG if default is UTF-8:
LC_ALL=C LC_TYPE=C LANG=C sed $'s/[^[:print:]\t]//g' ${INFILE} > ${OUTFILE}
###########################

##########################
# A few notes that may be helpful

#####
#For carriage returns (eg opening a terminal generated file in Mac and then saving will
# introduce carriage returns which Unix sees as ^M):
#From a Mac to Unix try:
#https://leemendelowitz.github.io/blog/remove-carriage-return-control-character.html
# tr '\r' '\n' < file_in > file_out
#####

#####
# Search for non-ascii in VIM
# https://stackoverflow.com/questions/16987362/how-to-get-vim-to-highlight-non-ascii-characters
#/[^\x00-\xFF]

#In VIM try:
#metacharacter ^M replacement with vim:
#https://unix.stackexchange.com/questions/32001/what-is-m-and-how-do-i-get-rid-of-it
#:s/^M$//
#:%s/^M/\r/g # this will insert a carriage return instead
#(Press Ctrl+V Ctrl+M to insert ^M.)
#####

#####
#Other options, try:
#http://stackoverflow.com/questions/3337936/remove-non-ascii-characters-from-csv

#cat -A FILENAME (to print all characters) # -v in Mac

#Check codes to character
#http://en.wikipedia.org/wiki/ASCII
#http://www.asciitable.com/

#sed -i 's/[\d128-\d255]//g' FILENAME
#or
#cat -A FILE | sed 's/[\d128-\d255]//g'
#####

#####
#http://alvinalexander.com/blog/post/linux-unix/how-remove-non-printable-ascii-characters-file-unix
#tr -cd '\11\12\15\40-\176' < file-with-binary-chars > clean-file
#####

#####
# Grep non-ascii:
# https://stackoverflow.com/questions/3001177/how-do-i-grep-for-all-non-ascii-characters
#pcregrep --color='auto' -n '[^\x00-\x7F]' folder/*
#grep --color='auto' -P -n "[\x80-\xFF]" file.xml # Linux
#####
##########################
