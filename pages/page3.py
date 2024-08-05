
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
ddg_info = pd.read_csv("ddg_infoNOTCH1.csv")
gene_pdbs = pd.read_csv("gene_pdbsNOTCH1.csv")


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
            html.Div([
                "Gene: ",
                ## Dropdown for genes
                dcc.Dropdown(options=[{'label': gene, 'value': gene} for gene in gene_pdbs['name_of_gene'].unique() if gene != 'gene_name_value'],
                             id="gene_selected",
                             value="NOTCH1",
                             searchable=True,
                             placeholder="Select a gene...",
                             clearable=True),
            ], style={'padding': 5}),
            dcc.Store(id='pdb_values'),
            ## Gene graph
            dcc.Graph(id="gene_ddg"),
            ## Markdown below the gene graph
            dcc.Markdown(
                id='gene_ddg_markdown', 
                style={'width': '100%', 'height': 100, 'white-space': 'pre-line'},
            ),
        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),


        html.Div([
            html.Div([
                html.Div([
                    "Residue number: ",
                    dcc.Dropdown(id="residual_selected",
                                 value = 293,
                                 searchable=True,
                                 placeholder="Select residue number...",
                                 clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
                html.Div([
                    "Mutation from: ",
                    dcc.Dropdown(id="mutfrom_selected",
                                 searchable=True,
                                 placeholder="Select mutation from...",
                                 clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
                html.Div([
                    "Mutation to: ",
                    dcc.Dropdown(id="mutto_selected",
                                 searchable=True,
                                 placeholder="Select mutation to...",
                                 clearable=True),
                ], style={'width': '32%', 'display': 'inline-block', 'padding': 5}),
            ], style={'display': 'flex', 'justify-content': 'space-between', 'width': '100%'}),
            dcc.Store(id='filtered_ddg_info'),
            dcc.Store(id='median_ddg'),
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
    Input(component_id="gene_selected", component_property="value")
)

def set_dropdown_options_page3_2a(gene_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info = ddg_info[(ddg_info['pdb'].isin(pdb_values))]
    return [{'label': str(variant), 'value': variant} for variant in filtered_ddg_info['pdb_residual'].unique() if variant != 'pdb_residual_value']

@app.callback(
    Output(component_id = "mutfrom_selected", component_property = "children"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value")]
)

def set_dropdown_options_page3_2b(gene_selected, residual_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info_var1 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected)]
    return [{'label': variant, 'value': variant} for variant in filtered_ddg_info_var1['mut_from'].unique() if variant != 'mut_from_value']

@app.callback(
    Output(component_id = "mutto_selected", component_property = "children"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value")]
)

def set_dropdown_options_page3_2c(gene_selected, residual_selected, mutfrom_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info_var2 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected)]
    return [{'label': variant, 'value': variant} for variant in filtered_ddg_info_var2['mut_to'].unique() if variant != 'mut_to_value']
 


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

# Calculate median of the variant histogram

def calculate_median (pdb_values,residual_selected, mutfrom_selected, mutto_selected):
    if mutfrom_selected and mutto_selected:
        filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
        median_ddg = filtered_ddg_info_var3['ddg'].median()
    else:
        median_ddg = None
    return median_ddg

def ddg_for_gene_plot(gene_selected, pdb_values, median_ddg):
  
    filtered_ddg_info = ddg_info[ddg_info['pdb'].isin(pdb_values)]

    figure = px.histogram(filtered_ddg_info, x='ddg', range_x=[-10, 100], nbins=1000, title=f'Delta Delta G values for {gene_selected}', labels={'ddg': 'ΔΔG (kcal/mol)', 'count': 'Frequency'})

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
                x=median_ddg,y=1,
                xref="x", yref="paper",
                text=f'Variant median: {median_ddg:.2f}',
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

def ddg_for_variant_plot(ddg_info, pdb_values, residual_selected=None, mutfrom_selected=None, mutto_selected=None):
    # Apply the initial filter based on pdb values
    filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
   
        # Create the histogram
    figure = px.histogram(filtered_ddg_info_var3, x='ddg', range_x=[-10, 100], nbins=20, title='Histogram of ddg values for selected variant', labels={'ddg': 'ΔΔG (kcal/mol)', 'count': 'Frequency'})
    return figure

##Callback for markdown text
@app.callback(
    Output(component_id = "gene_ddg_markdown", component_property = "children"),
    [Input(component_id = "gene_selected", component_property = "value"),
     Input(component_id = "pdb_values", component_property = "value"),
     Input(component_id = "residual_selected", component_property = "value"),
     Input(component_id = "mutfrom_selected", component_property = "value"),
     Input(component_id = "mutto_selected", component_property = "value"),]
)
def calculate_percentile(gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected):
    filtered_gene_pdbs = gene_pdbs[gene_pdbs['name_of_gene'] == gene_selected]
    pdb_values = filtered_gene_pdbs['pdb'].unique().tolist()
    filtered_ddg_info = ddg_info[(ddg_info['pdb'].isin(pdb_values))]
    values = filtered_ddg_info['ddg'].values
    filtered_ddg_info_var3 = ddg_info[(ddg_info['pdb'].isin(pdb_values)) & (ddg_info['pdb_residual'] == residual_selected) & (ddg_info['mut_from'] == mutfrom_selected) & (ddg_info['mut_to'] == mutto_selected)]
    median_ddg = filtered_ddg_info_var3['ddg'].median()
    percentile = np.sum(values < median_ddg) / len(values) * 100
    return percentile

def gene_ddg_markdown_text(gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected, median_ddg):
    percentile = calculate_percentile(gene_pdbs, gene_selected, ddg_info, residual_selected, mutfrom_selected, mutto_selected)
    if median_ddg > 2.5:
        return f'Median ΔΔG for selected variant is in the {percentile:.2f} percentile. It is greater than the deleterious value and destabilizing.'
    elif median_ddg > 0.5:
        return f'Median ΔΔG for selected variant is in the {percentile:.2f} percentile. It is greater than the Serrano value of 0.5.'
    else:
        return f'Median ΔΔG for selected variant is in the {percentile:.2f} percentile. It is not destabilizing.'


if __name__ == '__main__':
    app.run_server(debug=True, port=8052)
