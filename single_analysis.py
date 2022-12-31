import plotly.express as px
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import zscore
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from functools import reduce
import os

def file(virus, data_path):
    """
    file() returns the path to the gene_summary.txt file for a specific virus
    input the virus acronym and the path to where all the data is held. 
    When running locally, need to provide data path to dash app if not using data/mageck_nextflow_out/
    
    Attributes
    ----------
    virus : str
        the virus selected by the user
    datapath : str
        the file path to the gene_summary file for a specific virus

    Returns
    -------
    filepath: str
        the file path to the gene_summary file for a specific virus
    """
 
    for subdir, dirs, files in os.walk(data_path):
        for filename in files:
            filepath = subdir + os.sep + filename
            if filepath.endswith("gene_summary.txt") and virus in filepath: 
                return filepath

def get_gRNA(virus, data_path):
    """
    get_gRNA() returns the path to the .srgrna_summary.txt file for a specific virus
    input the virus acronym and the path to where all the data is held. 
    When running locally, need to provide data path to dash app if not using data/mageck_nextflow_out/
    
    Attributes
    ----------
    virus : str
        the virus selected by the user
    datapath : str
        the file path to the gene_summary file for a specific virus

    Returns
    -------
    filepath: str
        the file path to the gene_summary file for a specific virus
    """
 
    for subdir, dirs, files in os.walk(data_path):
        for filename in files:
            filepath = subdir + os.sep + filename
            if filepath.endswith(".sgrna_summary.txt") and virus in filepath: 
                return filepath  

def top_hits(filepath, num, metric):
    """
    top_hits() is a method used to return the top n inputted hits in the screen

    Attributes
    ----------
    filepath : str
        the file path to the gene_summary file for a specific virus
    num : int
        the number of significant genes, as selected by the user
    metric : str
        metric of significance (positive score, rank, etc.), as selected by the user

    Returns
    -------
    df : dataFrame
        has all the genes except for the significant ones
    df_max : dataFrame
        has all the significant genes
    """

    df = pd.read_csv(filepath, sep = '\t')
    df = df.rename(columns={'id':'Genes'})

    # compute negative logs into new columns
    df['-log10(pos|score)'] = -np.log10(df['pos|score'])
    df['-log10(neg|score)'] = -np.log10(df['neg|score'])
    df['pos|p-value'] = -np.log(df['pos|p-value'])
    df['neg|p-value'] = -np.log(df['neg|p-value'])
    df['Rank'] = df['pos|rank'].rank(axis = 0, ascending=True)

    # create df of top n significant genes
    df_max = df.nlargest(num, metric)

    df.drop(df_max.index, axis=0,inplace=True)
    df = df.sort_values(by='Genes', key=lambda col: col.str.lower())
    df_max = df_max.sort_values(by='Genes', key=lambda col: col.str.lower())
    df = df.reset_index(drop=True)
    df_max = df_max.reset_index(drop=True)

    return df, df_max

def sig_alpha(gen_df, df_max, gene_inputs, metric, virus):
    """
    sig_alpha() is a method used to sort dataframes based on significance 
    and assign color values to each gene for when the screen is plotted
    
    Attributes
    ----------
    gen_df and df_max : dataFrame
        two dfs returned from top_hits() 
    gene inputs : list
        list of genes the user wants to highlight in the graph (from callback)
    Genes_Alpha : list
        gives the x-coordinate of each point, genes are ordered alphabetically on the x-axis

    Returns
    -------
    grey_df : dataFrame
        df of all non-significant genes (without top n hit genes)
    red_df : dataFrame
        df of only significant genes
    final_df : dataFrame
        df of both significant and non significant genes
    """

    # grey_df is the data frame of non-significant genes
    grey_x = list()
    for i in gen_df.index:
        grey_x.append(i/len(gen_df['Genes']))
    
    grey_df = gen_df
    grey_df['Genes'] = gen_df['Genes']
    grey_df['Genes_Alpha'] = grey_x
    grey_df['Color'] = 'Not Enriched'

    # red_df is the data frame of significant genes (that will be in red on the graph)
    red_x = list()
    for i in df_max.index:
        red_x.append(i/len(df_max['Genes']))
    red_df = df_max
    red_df['Genes'] = df_max['Genes']
    red_df['Genes_Alpha'] = red_x
    red_df['Color'] = 'Enriched'
    
    # final_df is the concatination of grey and red df
    final_df = pd.concat([grey_df, red_df])
    final_df['Virus'] = [virus]*len(final_df['Genes'])
    
    # selected is a list of all genes the user selected
    for gene in gene_inputs:
        final_df.loc[(final_df["Genes"]==gene),"Color"] = "Selected"

    # returns the three concatonated data frames
    return red_df, grey_df, final_df

def sig_scatter_plot(df_gRNA, data_path, num, metric, gene_inputs, virus, backgroundColor):
    """
    Produce a plot of significance vs. alphabetical order

    Attributes
    ----------
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

    Returns 
    -------
    fig : Plotly graph
        scatter plot of a numerical vs categorical variable
    """

    # initialize dfs 
    filepath = file(virus, data_path)
    df, df_max = top_hits(filepath, num, metric)
    red_df, grey_df, final_df = sig_alpha(df, df_max, gene_inputs, metric, virus)

    # create figure
    fig = px.scatter(
        final_df, 
        x="Genes_Alpha", 
        y=metric, 
        color = 'Color',
        color_discrete_sequence=["grey", "red", "blue"], 
        labels={'Genes_Alpha':'Genes (Sorted Alphabetically)'}, 
        title="List of Selected Genes: " + str(gene_inputs).replace("'", "").replace("]", "").replace("[", ""),
        hover_data=
            {
                "Genes_Alpha": False,
                "Color":False,
                metric:False,
                "Genes": True,
                "pos|rank": True,
                "-log10(pos|score)": True,
                "pos|p-value": True,
                "pos|fdr": True,
                "pos|goodsgrna": True,
                "neg|rank": True,
                "neg|score": True,
                "neg|p-value": True,
                "neg|fdr": True,
                "neg|goodsgrna": True
            },
        template=backgroundColor
    )

    # add trace for hovering
    fig.add_trace(
        go.Scatter(
            x= red_df['Genes_Alpha'],
            y= red_df[metric],
            mode="text",
            name="Toggle Gene Names",
            text=red_df['Genes'],
            textposition="top center",
            hovertemplate = 
            f"{metric}:" +  "%{y}"
        )
    )

    # update CSS elements
    fig.update_layout(
        font_family = 'Inter',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        dragmode = 'select'
    )

    fig.update_xaxes(showticklabels=False)

    return fig

def sig_negative_log_rank_plot(df_gRNA, data_path, num, metric_y, gene_inputs, virus, backgroundColor):
    """
    Plot of one metric against another, initially will show singificance vs. rank
    
    Attributes
    ----------
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

    Returns 
    -------
    fig : Plotly graph
        scatter plot of a numerical vs numericals variable
    """

    # initialize df
    filepath = file(virus, data_path)
    df, df_max = top_hits(filepath, num, metric_y)
    df['Color'] = 'Not Enriched'
    df_max['Color'] = 'Enriched'
    final_df = pd.concat([df, df_max])
    
    # assigning pos / neg ranks
    metric_x = 'Rank'
    if 'neg' in str(metric_y):
        metric_x = 'neg|rank'

    # selected is a list of all genes the user selected
    for gene in gene_inputs:
        final_df.loc[(final_df["Genes"]==gene),"Color"] = "Selected"

    # create initial plot
    fig = px.scatter(
        final_df, 
        x=metric_x, 
        y=metric_y, 
        color = 'Color',
        color_discrete_sequence=["grey", "red", "blue"],
        title = 'Significance Rank Plot',
        hover_data={
            metric_x: False,
            "Color":False,
            "Genes": True,
            "pos|rank": True,
            "pos|score": True,
            "pos|p-value": True,
            "pos|fdr": True,
            "neg|rank": True,
            "neg|score": True,
            "neg|p-value": True,
            "neg|fdr": True
        },
        template=backgroundColor

    )

    fig.update_layout(
        font_family = 'Inter',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    return fig

"""Volcano and Table functions"""
def filter_df(df, df_gRNA, effect_size=[-5,1], pv=2.9957):
    """
    filter_df finds the limits of the effect size and 
    uses that to calculate and return a dataframe with 
    cols listed in the last lines of the functions

    Attributes
    ----------
    effect_size : list
        list containing a min/max for a specified range; can be altered
    pv : float
        the p-value cutoff

    Returns
    -------
    df : dataFrame
        dataFrame for a selected virus from the master dropdown
    """
    
    # create gRNA columns
    df['Positive gRNA'] = df['pos|goodsgrna']
    df['Negative gRNA'] = df['neg|goodsgrna']
    df['Positive gRNA z-score'] = round(zscore(df['pos|goodsgrna']), 4)
    df['Negative gRNA z-score'] = round(zscore(df['neg|goodsgrna']), 4)

    # filtering data frame
    df.lower_limit = effect_size[0]
    df.upper_limit = effect_size[1]
    df["color"] = np.where((df["pos|lfc"] > df.lower_limit) & (-np.log(df["pos|p-value"]) > pv), "Positive hit: Significant & Enriched", "Neither")
    df = df[["id", "pos|rank", "neg|p-value", "pos|p-value", 'pos|score' , 'pos|goodsgrna', 'neg|goodsgrna', 'Positive gRNA z-score', 'Negative gRNA z-score', "pos|lfc", "color"]]
    df = df.set_axis(['Gene id', 'Rank', "Negative P-Value",'Positive P-Value', 'Positive Score', 'Positive gRNA', 'Negative gRNA', 'Positive gRNA z-score', 'Negative gRNA z-score', 'LFC', 'Color'], axis=1, inplace=False)        
    df['Positive P-Value'] = round(-np.log(df['Positive P-Value']), 4)
    df['Negative P-Value'] = round(-np.log(df['Negative P-Value']), 4)
    df['Negative P-Value'] = df['Negative P-Value'].replace([0], 0.00000001)
    return df 

def get_sig_df(df):
    """
    Create table of significant and enriched genes 

    Returns
    -------
    df : dataFrame
        original df where the Color attribute ends in a specific string
    """
    return df[df.Color.str.endswith("Significant & Enriched")]

def positive_volcano_plot(df, backgroundColor):
    """
    Creates a Positive Volcano Plot

    Attributes
    ----------
    df : dataFrame
        initial df to be used for plotting
    backgroundColor : str
        string referencing which Plotly background color to be selected

    Returns
    -------
    fig : Plotly figure
        Positive volcano plot
    """
    df.drop(['Negative P-Value', 'Negative gRNA', 'Negative gRNA z-score'], axis = 1,inplace = True)

    fig = px.scatter(
        df,
        x= "LFC",
        y= 'Positive P-Value',
        color="Color",
        color_discrete_sequence=["red", "gray"],
        hover_data=["Gene id", "Positive P-Value", "LFC"],
        template=backgroundColor

    )

    fig.update_layout(
        xaxis_title="Log Fold Change",
        yaxis_title="-Log(Positive P-Value)",
        legend_title="Color Legend",
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    return fig

def negative_volcano_plot(df, lfc_range, p_value, backgroundColor):
    """
    Creates a Negative Volcano Plot

    Attributes
    ----------
    df : dataFrame
        initial df to be used for plotting
    lfc_range: list
        list containing a min/max for a specified range; can be altered
    p_value : float
        p_value specified
    backgroundColor : str
        string referencing which Plotly background color to be selected

    Returns
    -------
    fig : Plotly figure
        Negative volcano plot
    """

    df["Color"] = np.where((df["LFC"] < lfc_range[0]) & ((df["Negative P-Value"]) > p_value), "Negative hit: Significant & Enriched",  "Neither")
    df.drop(['Positive P-Value', 'Positive gRNA', 'Positive gRNA z-score'], axis = 1,inplace = True)

    fig = px.scatter(
        df,
        x= "LFC",
        y= 'Negative P-Value',
        color="Color",
        color_discrete_sequence=["gray", "blue"],
        hover_data=["Gene id", "Negative P-Value", "LFC"],
        template=backgroundColor
    )

    fig.update_layout(
        xaxis_title="Log Fold Change",
        yaxis_title="-Log(Negative P-Value)",
        legend_title="Color Legend",
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    return fig

def sunburst_plot(df):
    """
    creates a sunburst (plotly) figure

    Attributes
    ----------
    df : dataFrame
        dataFrame of top n hit genes

    Returns
    -------
    fig : Plotly figure
        sunburst chart
    """

    df['Screen'] = 'Screen'

    fig = px.sunburst(
        df,
        path=["Screen", 'Gene id', 'Positive gRNA'],
        title='*Likely to be replaced by Box Whisker Plot* Sunburst Chart of Gene Metric Significance'
    )
    fig.update_layout(
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    return fig

def pie_plot(df):
    """
    creates a pie chart (plotly) 

    Attributes
    ----------
    df : dataFrame
        dataFrame of top n hit genes

    Returns
    -------
    fig : Plotly figure
        pie chart
    """    
    
    fig = px.pie(
        df,
        values='Positive P-Value',
        names='Gene id',
        title='Proportion of Top n Genes'
    )
    return fig

def swarmplot(df_top_n_hits, df_gRNA, gRNA_metric, background_color):
    """
    Strip Plot showing gRNA counts and ranks for the top n significant genes
    
    Attributes
    ----------
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

    Returns
    -------
    fig : Plotly figure
        strip plot showing gRNA counts and ranks
    """

    # Drop unneeded columns
    df_top_n_hits = df_top_n_hits[['Gene id', 'Rank', "Negative P-Value",'Positive P-Value', 'Positive Score', 'Positive gRNA', 'Negative gRNA', 'Positive gRNA z-score', 'Negative gRNA z-score', 'LFC', 'Color']]     
    scatter_3d(df_top_n_hits)
    gRNA_subsets = list()

    for gene in df_top_n_hits['Gene id']:
        df_gRNA_subset = df_gRNA[df_gRNA['Gene'] == gene]
        df_gRNA_subset = df_gRNA_subset[['sgrna', 'Gene', 'score', 'LFC', 'p.twosided']]
        gRNA_subsets.append(df_gRNA_subset)

    # Merge n dfs in gRNA_subsets
    merged_df = pd.concat(gRNA_subsets)
    merged_df['score'] = -np.log(merged_df['score'])
    merged_df = merged_df.loc[::-1]
    
    fig = px.strip(
        y = merged_df['Gene'],
        x = merged_df[gRNA_metric],
        data_frame = merged_df,
        template=background_color,
        title=f'Strip Chart of gRNA {gRNA_metric}',
        hover_data= {
            'sgrna': True
        }
    )

    fig.update_layout(
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        xaxis_title=f"log(gRNA {gRNA_metric})",
        yaxis_title=f"Top n Significant Genes for {gRNA_metric}",
    )

    return fig

def box_plot(df_top_n_hits, df_gRNA, gRNA_metric, background_color):
    """
    Box Plot showing gRNA counts and ranks for the top n significant genes
    
    Attributes
    ----------
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

    Returns
    -------
    fig : Plotly figure
        Box plot showing gRNA counts, quartiles, and ranks
    """

    # Drop unneeded columns
    df_top_n_hits = df_top_n_hits[['Gene id', 'Rank', "Negative P-Value",'Positive P-Value', 'Positive Score', 'Positive gRNA', 'Negative gRNA', 'Positive gRNA z-score', 'Negative gRNA z-score', 'LFC', 'Color']]     

    gRNA_subsets = list()

    for gene in df_top_n_hits['Gene id']:
        df_gRNA_subset = df_gRNA[df_gRNA['Gene'] == gene]
        df_gRNA_subset = df_gRNA_subset[['sgrna', 'Gene', 'score', 'LFC', 'p.twosided']]
        gRNA_subsets.append(df_gRNA_subset)

    # Merge n dfs in gRNA_subsets
    merged_df = pd.concat(gRNA_subsets)
    merged_df['score'] = -np.log(merged_df['score'])
    merged_df = merged_df.loc[::-1]

    fig = px.box(
        y = merged_df['Gene'],
        x = merged_df['score'],
        data_frame = merged_df,
        template=background_color,
        title=f'Box & Whisker Plot of gRNA {gRNA_metric}',
        hover_data= {
            'sgrna': True
        },
        points=False
        )

    fig.update_layout(
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        xaxis_title=f"log(gRNA {gRNA_metric})",
        yaxis_title=f"Top n Significant Genes for {gRNA_metric}",
    )

    return fig

def scatter_3d(df):
    """
    scatter_3d() plots 3 metrics for a simple visualization
    This plot is used to generate images included in figures/
    It is not essential to running the app

    Attributes: 
    -----------
    df: dataFrame
        df for the user inputted virus

    Returns
    -------
    fig : Plotly figure
        plot showing 3 different metrics plotted
    """

    fig = px.scatter_3d(df, x = 'Positive Score', y = 'Positive gRNA z-score', z = 'Rank', color='Gene id')

    fig.update_layout(
        xaxis_title="Positive Score",
        yaxis_title="Positive gRNA z-score",
        legend_title="Rank",
        font_family = "Inter",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    return fig