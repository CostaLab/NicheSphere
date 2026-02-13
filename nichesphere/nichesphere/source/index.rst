.. Nichesphere documentation master file, created by
   sphinx-quickstart on Mon Feb  3 18:13:32 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NicheSphere's documentation!
=======================================

**NicheSphere** is an sc-verse compatible Python library which allows the user to find differential co-localization domains / niches and associated biological processes based on cell type pair co-localization probabilities in different conditions. 


.. toctree::
   :maxdepth: 2
   :caption: About NicheSphere

   about


.. toctree::
   :maxdepth: 2
   :caption: Installation

   installation


Tutorials
=========

In our first example we will use data from the Myocardial Infarction atlas from Kuppe, C. et. Al., 2022 to find differential co-localization domains related to ischemia. The data you'll need to run the tutorials can be found here: https://doi.org/10.5281/zenodo.15790389  

.. nbgallery::
    :caption: Notebooks:
    :glob:

    notebooks/Nichesphere_tutorial_MIvisium_coloc.ipynb
    notebooks/Nichesphere_tutorial_MIvisium_comm.ipynb
    notebooks/Nichesphere_tutorial_MIvisium_preprocessing.ipynb


Docker image and summarized analysis tutorial
=============================================

We provide access to a Docker image, available at: https://gitlab.com/sysbiobig/ismb-eccb-2025-tutorial-vt3/container_registry. The Docker image comes preconfigured with all necessary libraries, tools, and software required to follow the hands-on exercises. Additionally, the repository at https://gitlab.com/sysbiobig/ismb-eccb-2025-tutorial-vt3 contains a summarized Nichesphere co-localization + communication analysis tutorial.


API
===

:doc:`coloc`

:doc:`comm`

:doc:`tl`


