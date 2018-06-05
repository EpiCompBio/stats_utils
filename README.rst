.. .. include:: substitution_vars.rst

.. GitHub doe not render rst substitutions

.. copy across your travis "build..." logo so that it appears in your Github page

.. .. image:: https://travis-ci.org/EpiCompBio/stats_utils.svg?branch=master
    :target: https://travis-ci.org/EpiCompBio/stats_utils

.. do the same for ReadtheDocs image:

.. note that if your project is called project_Super readthedocs will convert
.. it to project-super

.. .. image:: https://readthedocs.org/projects/stats_utils/badge/?version=latest
    :target: http://stats_utils.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

 .. Edit manually:

.. .. Zenodo gives a number instead, this needs to be put in manually here:
   .. image:: https://zenodo.org/badge/#######.svg
      :target: https://zenodo.org/badge/latestdoi/#####

################################################
stats_utils
################################################


.. The following is a modified template from RTD
    http://www.writethedocs.org/guide/writing/beginners-guide-to-docs/#id1

.. For a discussion/approach see 
    http://tom.preston-werner.com/2010/08/23/readme-driven-development.html

A collection of scripts for common procedures (e.g. PCA)


Requirements
------------

Various and this will probably get outdated quickly. Please see the individual script requirements.

See also requirements files and Dockerfile for more information.

Most of the scripts are R or Python though so at the least you'll need:

* R >= 3.2
* Python >= 3.5
* r-docopt
* r-data.table
* r-ggplot2


Installation
------------

.. code-block:: bash
   
    pip install git+git://github.com/EpiCompBio/stats_utils.git


To use
------

.. code-block:: bash

    # Create a folder or a whole data science project, e.g. project_quickstart -n my_project
    cd my_project/results
    mkdir tests
    cd tests
    # Paths aren't set so you'll need to add full paths to each script and have
    the necessary dependencies installed
    # Simulate some data:
    python simulate_cont_var.py --createDF --sample-size=10000 --var-size=2000 -O cont_var_sim_data
    # Run principal components on it:
    Rscript run_PCA.R -h
    Rscript run_PCA.R -I test_file.tsv -O my_PCA
    # Check the outputs: 
    head my_PCA.tsv
    open my_PCA.svg

Contribute
----------

- Issue Tracker: github.com/EpiCompBio/stats_utils/issues
- Source Code: github.com/EpiCompBio/stats_utils
- Pull requests welcome!


Support
-------

If you have any issues, pull requests, etc. please report them in the issue tracker. 


