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
The script can be simply run by the command:
    $ python user_search_screens_by_gene.py
-----

