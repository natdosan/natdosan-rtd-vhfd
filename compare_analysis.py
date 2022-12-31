import plotly.express as px
import pandas as pd
import numpy as np
from scipy.stats import zscore
import plotly.graph_objects as go
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
    filepath : str
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
    filepath : str
        the file path to the gene_summary file for a specific virus
    """
 
    for subdir, dirs, files in os.walk(data_path):
        for filename in files:
            filepath = subdir + os.sep + filename
            if filepath.endswith(".sgrna_summary.txt") and virus in filepath: 
                print(filepath, virus)  
                return filepath 

def data_clean(df_merge, x_param, y_param):
    """
    data_clean() cleans the data inside the dataTable columns and returns the merged df

    Attributes
    ----------
    x_param : str 
        virus 1 parameter (user selected)
    y_param : str 
        virus 2 parameter (user selected)

    Returns
    -------
    df_merge : dataFrame
        merged and cleaned dataFrame of all the parameters
    """
    
    # drop all cols except for the one of the user inputted parameter
    df_merge = df_merge[['id','pos|lfc_x', 'pos|lfc_y', 'neg|goodsgrna_x', 'pos|goodsgrna_x', 'neg|goodsgrna_y', 'pos|goodsgrna_y', x_param, y_param]]

    # create ratio column for coloring legend
    # ratio is the x param / y param; values further from 1 will be more saturated, and closer to 1 will be more opaque
    df_merge['ratio'] = round(df_merge.apply( lambda row: row[x_param] / row[y_param] if row[y_param] != 0 else row[x_param] / .000000001, axis = 1), 4)

    # round decimal to 4 places
    df_merge[x_param] = round(df_merge[x_param], 4)
    df_merge[y_param] = round(df_merge[y_param], 4)

    # z-score normalization of gRNA
    df_merge['Positive gRNA z-score Virus 1'] = round(zscore(df_merge['pos|goodsgrna_x']), 4)
    df_merge['Positive gRNA z-score Virus 2'] = round(zscore(df_merge['pos|goodsgrna_y']), 4)
    df_merge['Negative gRNA z-score Virus 1'] = round(zscore(df_merge['neg|goodsgrna_x']), 4)
    df_merge['Negative gRNA z-score Virus 2'] = round(zscore(df_merge['neg|goodsgrna_y']), 4)
    print(df_merge.head(), df_merge.columns)

    return df_merge

# plot one virus against another based on user inputted parameter
def compare_plot(data_path, input_genes, parameter, virus1, virus2, scale, background_color):
    """
    compare_plot is a method to create a screen comparison scatter plot

    Attributes
    ----------
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
    
    Returns
    -------
    fig : Plotly figure
        scatter plot between two metrics
    """

    # takes in datapath -> uses file() to retrieve filepath -> filepath used to extract virus into dataFrame
    parameter = str(parameter).replace('[', '').replace(']', '').replace("'", '')
    df1 = pd.read_csv(file(virus1, data_path), sep="\t")
    df2 = pd.read_csv(file(virus2, data_path), sep="\t")
    df1['-log10 (pos|score)'] = -np.log10(df1['pos|score'])
    df2['-log10 (pos|score)'] = -np.log10(df2['pos|score'])
    df1['-log10 (neg|score)'] = -np.log10(df1['neg|score'])
    df2['-log10 (neg|score)'] = -np.log10(df2['neg|score'])

    # merge dataFrames into one with all shared genes 
    df_merge = df1.merge(df2, left_on='id', right_on='id', how='inner')
    df_merge.rename(columns={'id': 'Gene'}, inplace=True)
    x_param = parameter + "_x"
    y_param = parameter + "_y"
    df_merge = df_merge.loc[:, ['Gene', x_param, y_param]]
    df_merge['ratio'] = df_merge.apply( lambda row: row[x_param] / row[y_param] if row[y_param] != 0 else row[x_param] / .000000001, axis = 1)  

    # selected genes
    df_merge['Color'] = 'Unselected'
    selected = list()
    for gene in df_merge['Gene']:
        if gene in input_genes:
            selected.append(gene)
        else:
            selected.append("")
    df_merge['User_Selected'] = selected
    print("List of Selected Genes: " + str(input_genes).replace("'", "").replace("]", "").replace("[", ""))
    
    # main scatter plot
    fig = px.scatter(
        df_merge, 
        x = df_merge[x_param], 
        y = df_merge[y_param], 
        width=700, height=650, 
        text=df_merge['User_Selected'],
        color='ratio',
        color_discrete_sequence=["rgb(180, 151, 231)", "rgb(158, 185, 243)"],
        color_continuous_scale=scale,
        range_color=[0, 2],
        hover_data=[df_merge['Gene'], df_merge[x_param], df_merge[y_param], df_merge['ratio']],
        template=background_color
        )

    #fig.add_trace(
    #    go.Scatter(
    #        x = df_merge[x_param],
    #        y = df_merge[y_param],
    #        mode = "text",
    #        name = "Toggle Selected Genes",
    #        text = df_merge['User_Selected'],
    #        textposition = "top center",
    #        hoverinfo=True
    #    )
    #)

    fig.update_xaxes(
        title_text = str(x_param + " of Virus 1")
    )

    fig.update_yaxes(
        title_text = str(y_param + " of Virus 2")
    )

    fig.update_layout(
        font_family = 'Inter',
        font_color = 'black',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
    )
    
    return fig


