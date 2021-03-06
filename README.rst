

.. copy across your travis "build..." logo so that it appears in your Github page

.. image:: https://travis-ci.org/EpiCompBio/stats_utils.svg?branch=master
    :target: https://travis-ci.org/EpiCompBio/stats_utils

.. do the same for ReadtheDocs image:

.. note that if your project is called project_Super readthedocs will convert
.. it to project-super

.. image:: https://readthedocs.org/projects/stats_utils/badge/?version=latest
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
    
    # You may want to create a specific environment with conda first, then run:
    pip install git+git://github.com/EpiCompBio/stats_utils.git


To use
------

.. code-block:: bash

    # Create a folder or a whole data science project, e.g. project_quickstart -n my_project
    cd my_project/results
    mkdir tests
    cd tests
    # You may need to install missing dependencies, e.g.:
    conda install r-docopt r-data.table r-ggplot2 r-cowplot r-ggthemes
    # Simulate some data:
    simulate_cont_var.py -h
    simulate_cont_var.py --createDF --sample-size=1000 --var-size=50 -O cont_var_sim_data
    # The file will have rows as features/variables and columns as
    # samples/individuals. Transpose it for prcomp in run_PCA:
    transpose.R -I cont_var_sim_data.tsv
    # Run principal components:
    run_PCA.R -h
    run_PCA.R -I cont_var_sim_data.transposed.tsv
    # Check the outputs: 
    head cont_var_sim_data* | cut -f1-5
    open top_10_PCs_cont_var_sim_data.transposed.pca.svg


Contribute
----------

- `Issue Tracker`_

.. _`Issue Tracker`: github.com/EpiCompBio/stats_utils/issues

Pull requests welcome!


Support
-------

If you have any issues, pull requests, etc. please report them in the issue tracker. 
