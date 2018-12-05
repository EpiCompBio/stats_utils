##################################################
# Dockerfile for
# https://github.com/EpiCompBio/stats_utils
##################################################


############
# Base image
############

FROM continuumio/miniconda3
# It runs on Debian GNU/Linux 8; use e.g. uname -a ; cat /etc/issue.net
# https://hub.docker.com/r/continuumio/miniconda/
# Or simply run:
# docker run --rm -ti continuumio/miniconda3
# docker run --rm -ti ubuntu


#########
# Contact
#########
MAINTAINER Antonio Berlanga-Taylor <a.berlanga@imperial.ac.uk>


#########################
# Update/install packages
#########################

# Install system dependencies
# If running on Debian and anaconda/miniconda image, use apt-get:
RUN apt-get update && apt-get upgrade -qy apt-utils

RUN apt-get install -qy gcc \
    g++ \
    tzdata \
    wget \
    bzip2 \
    unzip \
    sudo \
    bash \
    fixincludes

# Get plotting libraries:
RUN apt-get install -qy \
            inkscape \
            graphviz

#########################
# Install conda
#########################

# Miniconda:
#RUN cd /usr/bin \
#    && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda

#RUN export PATH="/usr/local/miniconda/bin:$PATH"

# Add conda channels, last to be added take priority
# Don't mix conda-forge and/or bioconda with defaults channel in R as packages
# will conflict with other and fail
# channels
RUN conda config --add channels bioconda
    && conda config --add channels conda-forge
    && conda config --remove channels defaults
    && conda config --remove channels r

# Update conda:
RUN conda update -y conda

#########################
# Install dependencies
#########################

##########
# Install all packages needed
# Major packages:
RUN conda install python=3.6
    && conda install -y r=3.4

# Install python packages:
RUN pip install --upgrade pip numpy
    && pip install cython
    && pip install pandas

# Install project specific packages:
#RUN conda create -y -n stats_utils r=3.4
RUN conda install -y r-mice r-docopt r-ggplot2 r-dplyr r-data.table r-caret r-lattice r-VIM r-ggthemes
RUN R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("cowplot", repos = "http://cran.us.r-project.org") ; library("cowplot")' 
RUN R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("svglite", repos = "http://cran.us.r-project.org") ; library("svglite")' 

# Note parallel is a core R package, not on CRAN, and only needs to be loaded, as already installed with base R

# Install rpy2 with conda as pip version causes conflicts:
#RUN conda install -y rpy2
##########

##########
# Install packages that aren't on pip or conda:

# eg:
# FlashPCA2 and dependencies:

# Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
# Eigen has a conda recipe though: conda install -c menpo eigen
#RUN cd ; \
#    mkdir flashpca2 ; \
#    wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz ; \
#    tar xvfz 3.3.4.tar.gz \
#    mv eigen-eigen-5a0156e40feb eigen ; \
#    cd eigen ;
##########


##############################
# Install package of interest
##############################

RUN pip install git+git://github.com/EpiCompBio/stats_utils.git

############################
# Default action to start in
############################
# Only one CMD is read (if several only the last one is executed)
#ENTRYPOINT ['/xxx']
#CMD echo "Hello world"
#CMD project_quickstart.py
#CMD ["/bin/bash"]
CMD ["/bin/bash"]

# To build run as:
#docker build --no-cache=true -t antoniojbt/pipe_tests_alpine .

# To run e.g.:
# docker run --rm -ti antoniojbt/pipe_tests

# If mounting a volume do e.g.:
# docker run -v /host/directory:/container/directory --rm -ti antoniojbt/pipe_tests
# docker run -v ~/Documents/github.dir/docker_tests.dir:/home/ --rm -ti antoniojbt/pipe_tests_alpine

# Create a shared folder between docker container and host
#VOLUME ["/shared/data"]
