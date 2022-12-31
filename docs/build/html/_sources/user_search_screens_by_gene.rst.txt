Gene Search Functions (user_search_screens_by_gene)
==============================================================================

Functions
=========
.. automodule:: user_search_screens_by_gene
.. autofunction:: compile_dfs

=====
compile_dfs()
=====

    compile_dfs() creates a list of dataFrames corresponding to each virus in [viruses]
    
    **Attributes**
    
    virus : str
        the end of the filepath corresponding to where each virus is in the filepath
    data_path : str
        the user specific datapath given when running the script
    dfs : list
        list of dataFrames, one for each virus

    **Returns**
    
    dfs : list
        list of dataFrames, one for each virus
    dictionary_of_viruses_dfs : dict
        dict of dataFrames, one for each virus
        will never be empty; by default Ebola virus will be selected  

.. automodule:: user_search_screens_by_gene
.. autofunction:: dict_gene_searcxh

=====
dict_gene_search()
=====

    dict_gene_search() iterates through the dataFrames for each virus in {viruses} 
    and finds genes across all screens who's selected metric meets a threshold criteria
    then returns those genes in a single dataFrame for ease of observation
    
    **Attributes**
    
    data : list
        list of rows containing the gene and its selected threshold that meets 
        the user inputted threshold/condition
    dictionary_of_viruses_dfs : dict
        dict of dataFrames, one for each virus
        will never be empty; by default Ebola virus will be selected  
    gene : str
        user inputted gene to be searched throughout each screen
    metric: str
        user inputted metric
    metric_range : list
        range of selected metric aka the threshold
    gene_metric : series
        element within [data] that contains value of metric for specific gene
    list_of_dfs : list
        list of all concatonated dfs
    
    **Returns**

    final_df : df
        df of concatonated series -> dfs

Usage
=====

Dependencies
============

**channels:**
  - conda-forge
  - defaults

**dependencies:**
  - blas=1.0
  - bottleneck=1.3.4
  - brotlipy=0.7.0
  - bzip2=1.0.8
  - ca-certificates=2022.5.18.1
  - certifi=2022.5.18.1
  - cffi=1.15.0
  - click=8.0.4
  - dash=2.6.2
  - dash-core-components=2.6.2
  - dash-html-components=2.6.2
  - dataclasses=0.8
  - flask=2.0.3
  - flask-compress=1.5.0
  - icu=68.1
  - intel-openmp=2021.4.0
  - itsdangerous=2.0.1
  - jinja2=3.0.3
  - libcxx=12.0.0
  - libffi=3.3
  - libuv=1.40.0
  - markupsafe=2.0.1
  - mkl=2021.4.0
  - mkl-service=2.4.0
  - mkl_fft=1.3.1
  - mkl_random=1.2.2
  - ncurses=6.3
  - nodeenv=1.6.0
  - nodejs=16.13.1
  - numexpr=2.8.1
  - numpy=1.22.3
  - numpy-base=1.22.3
  - openssl=1.1.1o
  - packaging=21.3
  - pandas=1.4.2
  - pandas-stubs=1.2.0.58
  - pip=21.2.4
  - plotly=5.9.0
  - pycparser=2.21
  - pyparsing=3.0.4
  - pyright=1.1.248
  - python=3.10.8
  - python-dateutil=2.8.2
  - python_abi=3.10
  - pytz=2021.3
  - readline=8.1.2
  - setuptools=61.2.0
  - scipy=1.9.0
  - six=1.16.0
  - sqlite=3.38.3
  - tenacity=8.0.1
  - tk=8.6.11
  - tzdata=2022a
  - werkzeug=2.0.3
  - wheel=0.37.1
  - xz=5.2.5
  - zlib=1.2.12