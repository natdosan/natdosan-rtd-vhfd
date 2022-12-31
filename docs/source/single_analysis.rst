Single Screen Analysis Plots (single_analysis)
==============================================================================

Functions
=========
.. currentmodule:: single_analysis
.. autofunction:: file

======
file()
======
    file() returns the path to the gene_summary.txt file for a specific virus
    input the virus acronym and the path to where all the data is held. 
    When running locally, need to provide data path to dash app if not using data/mageck_nextflow_out

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
.. autofunction:: top_hits

=====
top_hits()
=====

    top_hits() is a method used to return the top n inputted hits in the screen

    **Attributes**

    filepath : str
        the file path to the gene_summary file for a specific virus
    num : int
        the number of significant genes, as selected by the user
    metric : str
        metric of significance (positive score, rank, etc.), as selected by the user

    **Returns**

    df : dataFrame
        has all the genes except for the significant ones
    df_max : dataFrame
        has all the significant genes

.. currentmodule:: single_analysis
.. autofunction:: sig_alpha

=====
sig_alpha()
=====

    sig_alpha() is a method used to sort dataframes based on significance 
    and assign color values to each gene for when the screen is plotted
    
    **Attributes**
    
    gen_df and df_max : dataFrame
        two dfs returned from top_hits() 
    gene inputs : list
        list of genes the user wants to highlight in the graph (from callback)
    Genes_Alpha : list
        gives the x-coordinate of each point, genes are ordered alphabetically on the x-axis

    **Returns**
    
    grey_df : dataFrame
        df of all non-significant genes (without top n hit genes)
    red_df : dataFrame
        df of only significant genes
    final_df : dataFrame
        df of both significant and non significant genes

.. currentmodule:: single_analysis
.. autofunction:: sig_scatter_plot

=====
sig_scatter_plot()
=====

    Produce a plot of significance vs. alphabetical order

    **Attributes**
    
    data_path : str
        the file path to the gene_summary file for a specific virus
    num : int
        number of user inputted top hits 
    metric : str
        string corresponding to the metric chosen for analysis
    gene_inputs : list
        list of selected genes to be highlighted
    virus : str
        string corresponding to the virus in the virus dictionary in the main single_complete-app file
    backgroundColor : str
        string referencing which Plotly background color to be selected

    **Returns**
    
    fig : Plotly graph
        scatter plot of a numerical vs categorical variable

.. currentmodule:: single_analysis
.. autofunction:: sig_negative_log_rank_plot

=====
sig_negative_log_rank_plot()
=====

    Plot of one metric against another, initially will show singificance vs. rank
    
    **Attributes**
    
    metric_x : str 
        the metric that will be plotted on the x-axis
    metric_y : str 
        the metric that will be plotted on the y-axis
    data_path : str
        the file path to the gene_summary file for a specific virus
    num : int
        number of user inputted top hits 
    metric : str
        string corresponding to the metric chosen for analysis
    gene_inputs : list
        list of selected genes to be highlighted
    virus : str
        string corresponding to the virus in the virus dictionary in the main single_complete-app file
    backgroundColor : str
        string referencing which Plotly background color to be selected

    **Returns**
    
    fig : plotly.express.scatter
        scatter plot of a numerical vs numericals variable

.. currentmodule:: single_analysis
.. autofunction:: filter_df

=====
filter_df()
=====

    filter_df finds the limits of the effect size and 
    uses that to calculate and return a dataframe with 
    cols listed in the last lines of the functions

    **Attributes**
    
    effect_size : list
        list containing a min/max for a specified range; can be altered
    pv : float
        the p-value cutoff

    **Returns**
    
    df : dataFrame
        dataFrame for a selected virus from the master dropdown

.. currentmodule:: single_analysis
.. autofunction:: get_sig_df

=====
get_sig_df()
=====

    Create table of significant and enriched genes 

    **Returns**
    
    df : dataFrame
        original df where the Color attribute ends in a specific string

.. currentmodule:: single_analysis
.. autofunction:: positive_volcano_plot

=====
positive_volcano_plot()
=====

    Creates a Positive Volcano Plot

    **Attributes**
    
    df : dataFrame
        initial df to be used for plotting
    backgroundColor : str
        string referencing which Plotly background color to be selected

    **Returns**

    fig : plotly.express.scatter
        Positive volcano plot

.. currentmodule:: single_analysis
.. autofunction:: negative_volcano_plot

=====
negative_volcano_plot()
=====

    Creates a Negative Volcano Plot

    **Attributes**

    df : dataFrame
        initial df to be used for plotting
    lfc_range: list
        list containing a min/max for a specified range; can be altered
    p_value : float
        p_value specified
    backgroundColor : str
        string referencing which Plotly background color to be selected

    **Returns**
    
    fig : plotly.express.scatter
        Negative volcano plot

.. currentmodule:: single_analysis
.. autofunction:: pie_plot

=====
pie_plot()
=====

    creates a pie chart (plotly) 

    **Attributes**
    
    df : dataFrame
        dataFrame of top n hit genes

    **Returns**
    
    fig : Plotly figure
        pie chart

.. currentmodule:: single_analysis
.. autofunction:: sunburst_plot

=====
sunburst_plot()
=====

    creates a sunburst (plotly) figure

    **Attributes**
    
    df : dataFrame
        dataFrame of top n hit genes

    **Returns**
    
    fig : plotly.express.sunburst
        sunburst chart
    
.. currentmodule:: single_analysis
.. autofunction:: swarmplot

=====
swarmplot()
=====

    Strip Plot showing gRNA counts and ranks for the top n significant genes
    
    **Attributes**
    
    df_top_n_hits : dataFrame 
        dataFrame containing n rows corresponding with the top n significant user selected genes
    df_gRNA : dataFrame 
        dataFrame containing the sgrna names for all the genes
    background_color : str
        name of the background color 
    gRNA_subsets : list
        list of top n gRNA dfs 
    df_gRNA_subset : dataFrame
        subset df for each Gene and its corresponding sgrna(s)
    merged_df : dataFrame
        concatonated df of all gRNA dfs

    **Returns**
    
    fig : plotly.express.strip
        strip plot showing gRNA counts and ranks

.. currentmodule:: single_analysis
.. autofunction:: box_plot

=====
box_plot()
=====

    Box Plot showing gRNA counts and ranks for the top n significant genes
    
    **Attributes**
    
    df_top_n_hits : dataFrame 
        dataFrame containing n rows corresponding with the top n significant user selected genes
    df_gRNA : dataFrame 
        dataFrame containing the sgrna names for all the genes
    background_color : str
        name of the background color 
    gRNA_subsets : list
        list of top n gRNA dfs 
    merged_df : dataFrame
        concatonated df of all gRNA dfs

    **Returns**
    
    fig : plotly.express.box
        Box plot showing gRNA counts, quartiles, and ranks

.. currentmodule:: single_analysis
.. autofunction:: scatter_3d

=====
scatter_3d()
=====

    scatter_3d() plots 3 metrics for a simple visualization
    This plot is used to generate images included in figures/
    It is not essential to running the app

    **Attributes:**
    
    df: dataFrame
        df for the user inputted virus

    **Returns**
    
    fig : plotly.express.scatter
        plot showing 3 different metrics plotted

Usage
=====
The app can be simply run by the command:
    $ python app.py

    the default port should be http://localhost:8080/

To view the single screens page, click on single-complete-app/, then

    - select the Virus Dataset
    - select the Analysis Metric
    - input any gene names into the "Search For Genes" search bar
-----