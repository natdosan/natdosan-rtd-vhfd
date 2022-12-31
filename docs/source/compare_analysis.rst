Dual Screen Analysis Plots (compare_analysis)
==============================================================================

Functions
=========
.. currentmodule:: compare_analysis
.. autofunction:: file

=====
file()
=====

    file() returns the path to the gene_summary.txt file for a specific virus
    input the virus acronym and the path to where all the data is held. 
    When running locally, need to provide data path to dash app if not using data/mageck_nextflow_out/
    
    **Attributes**
    
    virus : str
        the virus selected by the user
    datapath : str
        the file path to the gene_summary file for a specific virus

    **Returns**
    
    filepath : str
        the file path to the gene_summary file for a specific virus

.. currentmodule:: single_analysis
.. autofunction:: get_gRNA

=====
get_gRNA()
=====

    get_gRNA() returns the path to the .srgrna_summary.txt file for a specific virus
    input the virus acronym and the path to where all the data is held. 
    When running locally, need to provide data path to dash app if not using data/mageck_nextflow_out/

    **Attributes**

    virus : str
        the virus selected by the user
    datapath : str
        the file path to the gene_summary file for a specific virus

    **Returns**

    filepath : str
        the file path to the gene_summary file for a specific virus

.. currentmodule:: single_analysis
.. autofunction:: data_clean

=====
data_clean()
=====

    data_clean() cleans the data inside the dataTable columns and returns the merged df

    **Attributes**
    
    x_param : str 
        virus 1 parameter (user selected)
    y_param : str 
        virus 2 parameter (user selected)

    **Returns**
    
    df_merge : dataFrame
        merged and cleaned dataFrame of all the parameters

.. currentmodule:: single_analysis
.. autofunction:: compare_plot

=====
compare_plot()
=====

    compare_plot is a method to create a screen comparison scatter plot

    **Attributes**
    
    datapath : str
        filepath to the user selected virus
    input_genes : list
        list of all user selected genes
    parameter : str 
        user selected parameter to compare across screens
    virus1 : str
        1st selected virus
    virus 2 : str
        2nd selected virus
    scale : list
        list of all colorscales for the gradient plot
    background_color : str
        string from a list of strings; each string is the name of a plotly background color 
    
    **Returns**
    
    fig : plotly.express.scatter 
        scatter plot between two metrics

Usage
=====
The app can be simply run by the command:
    $ python app.py

    the default port should be http://localhost:8080/

To view the compare screens page, click on compare-app/, then
    - select the Virus Dataset for the X axis
    - select the Virus Dataset for the Y axis
    - select the Analysis Metric
    - input any gene names into the "Search For Genes" search bar
-----