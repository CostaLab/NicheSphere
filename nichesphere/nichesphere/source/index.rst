.. Nichesphere documentation master file, created by
   sphinx-quickstart on Mon Feb  3 18:13:32 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Nichesphere's documentation!
=======================================


Introduction
============

**Nichesphere** is an sc-verse compatible Python library which allows the user to find differential co-localization domains / niches based on cell type pair co-localization probabilities in different conditions. Cell type pair co-localization probabilities can be obtained in different ways, for example, through deconvolution of spatial transcriptomics / PIC-seq data (getting the probabilities of finding each cell type in each spot / multiplet) ; or counting cell boundaries overlaps for each cell type pair in single cell spatial data (MERFISH , CODEX ...).

It also offers the possibility to look at biological process based differential communication among differential co-localization domains based on Ligand-Receptor pairs expression data, such as results from CrossTalkeR [ref.].


Installation
============

In your terminal window run::

    conda create --name test python=3.10
    conda activate test
    conda install pip
    pip install jupyterlab
    conda install --channel conda-forge pygraphviz
    conda install conda-forge::git
    pip install git+https://github.com/CostaLab/Nichesphere#subdirectory=nichesphere


Tutorials
=========

In our first example we will use data from the Myocardial Infarction atlas from Kuppe, C. et. Al., 2022 to find differential co-localization domains related to ischemia. The data you'll need to run the tutorials can be found here: https://zenodo.org/records/15790389?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjU3MTkzYmVlLWI4ZWItNGNlZS1hZGJhLWRiNjAxMWU0MjNmNyIsImRhdGEiOnt9LCJyYW5kb20iOiI1MTkzMGYzZDA1YmIyODg4MjZmODc2NjNlMTFiODcxZiJ9.RDwYH0UnKsr8Kw03P5WRuL1Gbr2GL10IzUHW_i8LU3GNFYSlP2LmOwRRApLP2hnZiO-yrg88F12WJjKAU_nn7A 

.. nbgallery::
    :caption: Notebooks:
    :glob:

    notebooks/Nichesphere_tutorial_MIvisium_coloc.ipynb
    notebooks/Nichesphere_tutorial_MIvisium_comm.ipynb


API
===

:doc:`coloc`

:doc:`comm`

:doc:`tl`


