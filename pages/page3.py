
from dash import Dash, html, dcc, Output, Input
import dash, dash_table
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

app = Dash(__name__)

## Read in data
ddg_info = pd.read_csv("/Users/clementineelwes/research/Folding_Dashboard/ddg_info.csv")
gene_pdbs = pd.read_csv("/Users/clementineelwes/research/Folding_Dashboard/gene_pdbs.csv")



### ----------------------
# Layout
layout = html.Div(children=[
    html.Br(),
    html.H1(children='Folding Energies for Genes and Variants'),

    html.Div(children='''
        Here you can explore the (mis)folding energies.
    '''),


    html.Div([

        # Graph container
        html.Div([
            "Gene: ",
            ## Dropdown for genes
            dcc.Dropdown(options=[{'label': gene, 'value': gene} for gene in gene_pdbs['name_of_gene'].unique() if gene != 'gene_name_value'],
                         id="gene_selected",
                         searchable=True,
                         placeholder="Select a gene...",
                         clearable=True),
            dcc.Store(id='pdb_values'), 
            ## Gene graph
            dcc.Graph(id="gene_ddg"),
        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            "Residue number: ",
            dcc.Dropdown(options=[{'label': str(variant), 'value': variant} for variant in ddg_info['pdb_residual'].dropna().unique() if variant != 'pdb_residual_value'],
                         id="residual_selected",
                         searchable=True,
                         placeholder="Select residue number...",
                         clearable=True),
            "Mutation from: ",
            ## Dropdown for variants
            dcc.Dropdown(options=[{'label': variant, 'value': variant} for variant in ddg_info['mut_from'].unique() if variant != 'mut_from_value'],
                         id="mutfrom_selected",
                         searchable=True,
                         placeholder="Select mutation from...",
                         clearable=True),
            "Mutation to: ",
            ## Dropdown for variants
            dcc.Dropdown(options=[{'label': variant, 'value': variant} for variant in ddg_info['mut_to'].unique() if variant != 'mut_to_value'],
                         id="mutto_selected",
                         searchable=True,
                         placeholder="Select mutation to...",
                         clearable=True),
            
            ## Variant graph
            dcc.Graph(id="variant_ddg"),
        ], style={'width': '52%', 'display': 'inline-block', 'padding': 10})
    ], style={'display': 'flex'}),
])

#-------------------------------------------
## Callbacks

## Callback for Dropdown for gene
@app.callback(
    Output(component_id="gene_selected", component_property="children"),
    [Input(component_id="genes", component_property="value")]
)

@app.callback(
    Output(component_id="residual_selected", component_property="children"),
    Input(component_id="variants", component_property="value")
)

@app.callback(
    Output(component_id="mutfrom_selected", component_property="children"),
    Input(component_id="variants", component_property="value")
)

@app.callback(
    Output(component_id="mutto_selected", component_property="children"),
    Input(component_id="variants", component_property="value")
)


## Callback for ddg_for_gene
@app.callback(
    Output(component_id = "gene_ddg", component_property = "figure"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)

def get_pdb_values(gene_pdbs, gene_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    return pdb_values

def ddg_for_gene_plot(gene_selected, pdb_values, residual_selected, mutfrom_selected, mutto_selected):
  
    filtered_ddg_info = ddg_info[ddg_info['pdb'].isin(pdb_values)]
    
    # Calculate median of the variant histogram
    if mutfrom_selected and mutto_selected:
        filtered_ddg_info_var = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
        median_ddg = filtered_ddg_info_var['ddg'].median()
    else:
        median_ddg = None

    figure = px.histogram(filtered_ddg_info, x='ddg', nbins=1000, title='Histogram of ddg values for selected gene')

    # Add a vertical line at the median value
    if median_ddg:
        figure.add_shape(
            go.layout.Shape(
                type="line",
                x0=median_ddg,
                x1=median_ddg,
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                line=dict(
                    color="Red",
                    width=2
                )
            )
        )
        figure.add_annotation(
            go.layout.Annotation(
                x=median_ddg,
                y=max(filtered_ddg_info['ddg'].value_counts()) / 2,
                text=f'Median: {median_ddg:.2f}',
                showarrow=True,
                arrowhead=2
            )
        )

    return figure

## Callback for ddg_for_variant
@app.callback(
    Output(component_id = "variant_ddg", component_property = "figure"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)

def get_pdb_values(gene_pdbs, gene_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    return pdb_values

def ddg_for_variant_plot(ddg_info, pdb_values, residual_selected=None, mutfrom_selected=None, mutto_selected=None):
    # Apply the initial filter based on pdb values
    filtered_ddg_info_var = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
   
        # Create the histogram
    figure = px.histogram(filtered_ddg_info_var, x='ddg', nbins=20, title='Histogram of ddg values for selected variant')
    return figure




if __name__ == '__main__':
    app.run_server(debug=True, port=8052)
